"""B3 机制实验：受控顺序重拟合 MLP 二面角（2026-10-15 计划 §4 B3；用户 2026-10-17 选定）。

只读一个 G2 运行目录（IR + scans/），不加载 MACE、不调 xtb、不重扫，不改任何 provider 默认值。

锁定规格（运行前定，结果出来后不改）：
- 第 0 轮 = G2 的单遍结果（用缓存按 MLPTorsionProvider 流程重算）。
- 每轮按分子内可转键的固定顺序（report.rotors，按键号）逐键更新：去掉 r 上的全部 proper，
  y = E_ref_r − E_MM[当前 IR \\ r](x_r)，fit_fourier(φ_r, y, n_max=6)，项放在 r 的扫描四元组上；
  更新后立即替换（后面的键用新值）。**每次只用 r 自己的扫描**——这一般不是联合最小二乘
  的块坐标下降，只作机制实验，不预设收敛到联合解。
- 跑满 3 轮，不设收敛判据，不因好坏提前停。
- 每轮报告 12 根键：更新时的拟合残差、轮末 IR 在各自扫描上的剖面误差（G2 原口径 min-shift、
  最佳常数平移），以及全局残差 = 全部扫描逐帧（各扣自身常数）合并 RMS（即联合目标值，等权逐帧），
  另给不含 nma 1-3 的版本（扫描滞后 57 kJ/mol，单列为扫描质量问题）。

    python scripts/refit_torsions_sequential.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_b3_sequential_20261017 [--rounds 3]
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_torsion_profiles import (  # noqa: E402
    _bestshift_rmse,
    _load_fixture,
    _minshift_rmse,
    _on_bond,
    _openmm_profile,
    _sha256,
)

SCAN_QUALITY_EXCLUDED = {("nma", "1-3")}  # 用户 2026-10-17：单列扫描质量问题，全局残差另给不含它的版本


def _refit_bond(ir, bond, scan, meta_extra):
    from ffcompiler.torsion.fourier import fit_fourier, to_ir_terms
    from ffcompiler.torsion.mm import mm_energy_without_bond

    stripped = ir.replace_terms(proper_torsions=tuple(t for t in ir.proper_torsions if not _on_bond(t, bond)))
    y = np.asarray(scan.energy, dtype=np.float64) - mm_energy_without_bond(stripped, bond, scan.positions)
    fit = fit_fourier(scan.phi_rad, y)
    meta = {"scan_quartet": "-".join(map(str, scan.quartet)), "fit_rmse": fit.rmse, **meta_extra}
    terms = to_ir_terms(fit, scan.quartet, source="mlp", meta=meta)
    return stripped.replace_terms(proper_torsions=tuple(list(stripped.proper_torsions) + terms)), fit.rmse


def _bond_term_profile(ir, bond, x) -> np.ndarray:
    from ffcompiler.ir.energy import evaluate

    sub = ir.replace_terms(proper_torsions=tuple(t for t in ir.proper_torsions if _on_bond(t, bond)))
    return np.array([evaluate(sub, xi)["proper"] for xi in x])


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

    mols = {}
    for name in names:
        base, pos, mol = _load_fixture(name)
        provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                      scan_cache=g2 / "scans" / name, model_path=default_model_path(),
                                      cache_only=True)
        ir0, report = provider.parameterize(base, pos, mol)
        bonds = [tuple(r["bond"]) for r in report.rotors if r["status"] == "mlp"]
        fit0 = {tuple(r["bond"]): r["fit_rmse"] for r in report.rotors if r["status"] == "mlp"}
        mols[name] = {"base": base, "ir": ir0, "bonds": bonds, "scans": report.scans,
                      "fit": fit0, "terms_prev": {b: _bond_term_profile(ir0, b, report.scans[b].positions)
                                                  for b in bonds}}

    rows, global_rows = [], []

    def record(round_no: int) -> None:
        pooled, pooled_ex = [], []
        for name, m in mols.items():
            for b in m["bonds"]:
                label = f"{b[0]}-{b[1]}"
                scan = m["scans"][b]
                e_ref = np.asarray(scan.energy, dtype=np.float64)
                e_mm = _openmm_profile(m["ir"], scan.positions)
                d = e_mm - e_ref
                d = d - d.mean()
                pooled.append(d)
                if (name, label) not in SCAN_QUALITY_EXCLUDED:
                    pooled_ex.append(d)
                terms_now = _bond_term_profile(m["ir"], b, scan.positions)
                dterm = terms_now - m["terms_prev"][b]
                m["terms_prev"][b] = terms_now
                rows.append({
                    "round": round_no, "molecule": name, "bond": label,
                    "fit_rmse_at_update": m["fit"][b],
                    "minshift_rmse": _minshift_rmse(e_mm, e_ref),
                    "bestshift_rmse": _bestshift_rmse(e_mm, e_ref),
                    "gaff2_minshift_rmse": _minshift_rmse(_openmm_profile(m["base"], scan.positions), e_ref),
                    "term_change_rms_vs_prev_round": _bestshift_rmse(dterm, 0.0 * dterm),
                    "scan_quality_excluded": (name, label) in SCAN_QUALITY_EXCLUDED,
                })
        allv, exv = np.concatenate(pooled), np.concatenate(pooled_ex)
        global_rows.append({"round": round_no, "global_rms_all12": float(np.sqrt(np.mean(allv ** 2))),
                            "global_rms_excl_nma_1_3": float(np.sqrt(np.mean(exv ** 2))),
                            "n_frames_all12": int(allv.size), "n_frames_excl": int(exv.size)})
        print(f"[{time.time() - t0:6.1f}s] round {round_no}: global {global_rows[-1]['global_rms_all12']:.4f}"
              f" / excl nma1-3 {global_rows[-1]['global_rms_excl_nma_1_3']:.4f}", flush=True)

    record(0)
    for rnd in range(1, args.rounds + 1):
        for name, m in mols.items():
            for b in m["bonds"]:
                m["ir"], m["fit"][b] = _refit_bond(m["ir"], b, m["scans"][b], {"b3_round": rnd})
        record(rnd)

    with open(args.outdir / "per_bond.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(args.outdir / "global.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(global_rows[0]))
        w.writeheader()
        w.writerows(global_rows)
    from ffcompiler.ir import io as ir_io

    (args.outdir / "ir").mkdir(exist_ok=True)
    for name, m in mols.items():
        ir_io.to_json(m["ir"], args.outdir / "ir" / f"{name}_round{args.rounds}.json")

    sources = [Path(__file__), ROOT / "scripts" / "diagnose_torsion_profiles.py",
               *sorted((ROOT / "src" / "ffcompiler").rglob("*.py"))]
    manifest = {
        "task": "B3 sequential refit mechanism experiment (2026-10-15 plan §4 B3)",
        "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
        "wall_s": time.time() - t0, "g2_run": str(g2), "molecules": names, "rounds": args.rounds,
        "order": {n: [f"{b[0]}-{b[1]}" for b in m["bonds"]] for n, m in mols.items()},
        "scan_quality_excluded_from_second_global": sorted(f"{a} {b}" for a, b in SCAN_QUALITY_EXCLUDED),
        "inputs_sha256": {str(p): _sha256(p)
                          for p in [*sorted((g2 / "ir").glob("*.json")), *sorted((g2 / "scans").rglob("*.npz"))]},
        "source_sha256": {str(p.relative_to(ROOT)): _sha256(p) for p in sources},
        "git_commit": None, "note": "sha256 只作身份记录，不用于复用判断",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
