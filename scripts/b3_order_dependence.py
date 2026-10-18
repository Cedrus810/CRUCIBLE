"""B3 顺序依赖检查（2026-10-15 计划 §4 B3；用户 2026-10-18 要求，正式采用顺序重拟合之前）。

只读一个 G2 运行目录（IR + scans/），不加载 MACE、不调 xtb、不重扫，不改 provider 默认值。

锁定规格（运行前定，结果出来后不改）：
- 起点与更新规则与 ``refit_torsions_sequential.py`` 相同：第 0 轮 = G2 单遍结果；每次只用本键自己的扫描，
  其他键取当前项扣 MM，fit_fourier(n_max=6)，项放在扫描四元组上，更新后立即替换；跑满 3 轮。
- 每个分子枚举其 mlp 键的**全部顺序**（分子之间互不耦合，分别比较）。
- 每轮、每个分子报告：
  (a) 各键 best-shift 剖面 RMSE 在各顺序间的极差；
  (b) 各键拟合项能量曲线（在本键扫描上、扣常数）相对"正序"（键号升序）结果的最大 RMS 差；
  (c) 分子合并残差（各扫描扣自身常数后合并 RMS）在各顺序间的极差。
- 不设"顺序无关"的判定阈值，只报数。nma 1-3 照常计算，标注扫描质量问题（滞后 57 kJ/mol）。

    python scripts/b3_order_dependence.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_b3_order_20261018 [--rounds 3]
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_torsion_profiles import _bestshift_rmse, _load_fixture, _openmm_profile, _sha256  # noqa: E402
from refit_torsions_sequential import _bond_term_profile, _refit_bond  # noqa: E402

SCAN_QUALITY_FLAG = {("nma", "1-3")}


def _label(b) -> str:
    return f"{b[0]}-{b[1]}"


def _round_state(ir, bonds, scans) -> dict:
    """一轮结束时：各键 best-shift RMSE、各键项曲线（本键扫描上）、分子合并残差。"""
    best, terms, pooled = {}, {}, []
    for b in bonds:
        scan = scans[b]
        e_ref = np.asarray(scan.energy, dtype=np.float64)
        d = _openmm_profile(ir, scan.positions) - e_ref
        pooled.append(d - d.mean())
        best[b] = float(np.sqrt(np.mean((d - d.mean()) ** 2)))
        terms[b] = _bond_term_profile(ir, b, scan.positions)
    return {"best": best, "terms": terms,
            "pooled": float(np.sqrt(np.mean(np.concatenate(pooled) ** 2)))}


def main(argv=None) -> int:
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g2-run", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--rounds", type=int, default=3)
    args = ap.parse_args(argv)
    t0 = time.time()
    g2 = args.g2_run
    names = sorted(p.stem for p in (g2 / "ir").glob("*.json") if (g2 / "scans" / p.stem).is_dir())
    args.outdir.mkdir(parents=True, exist_ok=True)

    order_rows, summary_rows = [], []
    for name in names:
        base, pos, mol = _load_fixture(name)
        provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                      scan_cache=g2 / "scans" / name, model_path=default_model_path(),
                                      cache_only=True)
        ir0, report = provider.parameterize(base, pos, mol)
        bonds = [tuple(r["bond"]) for r in report.rotors if r["status"] == "mlp"]
        scans = report.scans
        forward = tuple(sorted(bonds))
        states: dict[tuple, list[dict]] = {}
        for order in itertools.permutations(forward):
            ir = ir0
            per_round = [_round_state(ir, bonds, scans)]
            for rnd in range(1, args.rounds + 1):
                for b in order:
                    ir, _ = _refit_bond(ir, b, scans[b], {"b3_round": rnd})
                per_round.append(_round_state(ir, bonds, scans))
            states[order] = per_round
            for rnd, st in enumerate(per_round):
                for b in bonds:
                    order_rows.append({"molecule": name, "order": " ".join(_label(x) for x in order),
                                       "round": rnd, "bond": _label(b), "bestshift_rmse": st["best"][b],
                                       "molecule_pooled_rms": st["pooled"],
                                       "scan_quality_flag": (name, _label(b)) in SCAN_QUALITY_FLAG})
        for rnd in range(args.rounds + 1):
            pooled = [states[o][rnd]["pooled"] for o in states]
            for b in bonds:
                best = [states[o][rnd]["best"][b] for o in states]
                ref_terms = states[forward][rnd]["terms"][b]
                term_dev = max(_bestshift_rmse(states[o][rnd]["terms"][b], ref_terms) for o in states)
                summary_rows.append({
                    "molecule": name, "bond": _label(b), "round": rnd, "n_orders": len(states),
                    "bestshift_min": min(best), "bestshift_max": max(best),
                    "bestshift_spread": max(best) - min(best),
                    "term_profile_max_rms_dev_vs_forward": term_dev,
                    "molecule_pooled_min": min(pooled), "molecule_pooled_max": max(pooled),
                    "molecule_pooled_spread": max(pooled) - min(pooled),
                    "scan_quality_flag": (name, _label(b)) in SCAN_QUALITY_FLAG,
                })
        print(f"[{time.time() - t0:6.1f}s] {name}: {len(states)} orders", flush=True)

    for fname, rows in (("per_order.csv", order_rows), ("summary.csv", summary_rows)):
        with open(args.outdir / fname, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    sources = [Path(__file__), ROOT / "scripts" / "refit_torsions_sequential.py",
               ROOT / "scripts" / "diagnose_torsion_profiles.py",
               *sorted((ROOT / "src" / "ffcompiler").rglob("*.py"))]
    manifest = {
        "task": "B3 order dependence check (2026-10-15 plan §4 B3)",
        "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
        "wall_s": time.time() - t0, "g2_run": str(g2), "molecules": names, "rounds": args.rounds,
        "inputs_sha256": {str(p): _sha256(p) for p in [*sorted((g2 / "ir").glob("*.json")),
                                                      *sorted((g2 / "scans").rglob("*.npz"))]},
        "source_sha256": {str(p.relative_to(ROOT)): _sha256(p) for p in sources},
        "git_commit": None, "note": "sha256 只作身份记录，不用于复用判断",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
