#!/usr/bin/env python
"""Gate G3：电荷 provider（载入 / AM1-BCC / NAGL）与换电荷后的二面角重拟合（MVP-3 计划 T37、§6）。

    python scripts/gate_mvp3.py --scan-cache G2_RUN_DIR/scans [--skip-compare] [--skip-pytest] [--molecules ...]
                                # 后三个只用于调试

在 RunContext("gate_mvp3") 下运行，写 report.md、acceptance.json；跨引擎比较存
compare/<provider>/<分子>/report.json。扫描缓存用 --scan-cache 显式指向本次 G2 的 scans/（只读，cache_only；无默认值）。
判据 fail-closed；NAGL 与 AM1-BCC 的差、换电荷后的剖面 RMSE 只报告。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import numpy as np  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "gaff2"
G31_NAMES = (
    "ethanol", "dmso", "nma", "methyl_benzoate", "acetate", "methylammonium",
    "benzene", "chlorobenzene", "cyclopropane", "cyclobutane", "norbornane", "acetonitrile",
)
NAMES = G31_NAMES[:-1]  # 乙腈（近线性角）不进 compare（决策 0013）
PROVIDERS = ("loaded", "am1bcc", "nagl")
ROTOR_COUNTS = {"ethanol": 2, "dmso": 2, "nma": 3, "methyl_benzoate": 3, "acetate": 1,
                "methylammonium": 1}  # MVP-2 计划 T25 用例 1；共 12 根
SEED = 20261004
TOL_CHARGE = 1e-10  # e；依据见计划 §6（mol2 10 位小数舍入 5e-11；均摊后 float64 求和误差 ~N·1e-16）

REPORT_HEADER = (
    "# Gate G3 评估报告（MVP-3，电荷 provider）\n\n"
    "NAGL 与 AM1-BCC 的差、换电荷后的二面角剖面 RMSE 只报告、不设门槛。\n"
)


def _key(pair) -> str:
    return "/".join(pair) if isinstance(pair, tuple) else str(pair)


def build_report(*, full_matrix: bool, g31: dict, g32: dict, compare: dict, sander_notes: dict, g34: dict,
                 nagl_vs_am1bcc: dict, profiles: list[dict], pytest_summary: dict, failures: list[dict],
                 expected_g31: list[str], expected_compare: list[tuple[str, str]],
                 expected_g34: dict[str, int]) -> tuple[str, dict]:
    """纯函数：各项结果 → (report.md, acceptance.json 内容)。缺数据一律判未通过。"""
    criteria = {
        "G3.1": all(g31.get(n) is not None and g31[n] <= TOL_CHARGE for n in expected_g31),
        "G3.2": all(g32.get(k) is not None and g32[k] <= TOL_CHARGE for k in expected_compare),
        "G3.3": all(compare.get(k) is True for k in expected_compare),
        "G3.4": bool(expected_g34) and all(
            n in g34 and g34[n].get("identical") is True
            and g34[n].get("nagl_statuses") == ["mlp"] * count
            for n, count in expected_g34.items()),
        "G3.5": all(pytest_summary.get(m, {}).get("accepted") is True for m in ("not external", "external")),
    }
    lines = [REPORT_HEADER]
    if not full_matrix:
        lines.append("\n> 子集运行（scope=subset）：判据结果只供调试，不是正式通过（决策 0019）。\n")
    lines += ["\n## G3.1 AM1-BCC vs 固件 mol2（逐原子最大差，e）\n", "| 分子 | 最大差 |", "|---|---|"]
    lines += [f"| {n} | {g31.get(n)} |" for n in expected_g31]
    lines += ["\n## G3.2 / G3.3 总电荷偏差与跨引擎比较\n", "| provider/分子 | \\|Σq−Q\\| | compare |", "|---|---|---|"]
    lines += [f"| {_key(k)} | {g32.get(k)} | {compare.get(k)} |" for k in expected_compare]
    lines.append("\n### 已知局限：sander 腿只报告、不计失败（决策 0017）\n")
    lines += [f"- {k}: {note}" for k in sorted(sander_notes) for note in sander_notes[k]] or ["- 无"]
    lines += ["\n## G3.4 二面角重拟合\n", "| 分子 | AM1-BCC 重跑后逐项相等 | NAGL 电荷下各键状态 |", "|---|---|---|"]
    lines += [f"| {n} | {g34.get(n, {}).get('identical')} | {g34.get(n, {}).get('nagl_statuses')} |"
              for n in expected_g34]
    lines += ["\n## 只报告：NAGL vs AM1-BCC（e）\n", "| 分子 | 最大差 | RMS |", "|---|---|---|"]
    lines += [f"| {n} | {v['max_abs']:.4g} | {v['rms']:.4g} |" for n, v in sorted(nagl_vs_am1bcc.items())]
    lines += ["\n## 只报告：剖面 RMSE（kJ/mol，新 IR 的 MM vs E_MLP，各减自身最小值）\n",
              "| 分子 | 键 | 电荷 | 剖面 RMSE |", "|---|---|---|---|"]
    lines += [f"| {p['molecule']} | {p['bond']} | {p['charges']} | {p['profile_rmse']:.4g} |" for p in profiles]
    lines += ["\n## 判据\n", "| 判据 | 结果 |", "|---|---|"]
    lines += [f"| {k} | {'通过' if ok else '未通过'} |" for k, ok in criteria.items()]
    lines.append("\n## 失败项\n")
    lines += [f"- {f.get('stage')} {f.get('molecule')}: {f.get('reason')}" for f in failures] or ["- 无"]
    acceptance = {
        "gate": "G3",
        "scope": "full" if full_matrix else "subset",
        "passed": full_matrix and all(criteria.values()) and not failures,
        "criteria": criteria,
        "g31": g31,
        "g32": {_key(k): v for k, v in g32.items()},
        "compare": {_key(k): v for k, v in compare.items()},
        "sander_known_limitation": sander_notes,
        "g34": g34,
        "nagl_vs_am1bcc": nagl_vs_am1bcc,
        "profiles": profiles,
        "failures": failures,
        "pytest_summary": pytest_summary,
    }
    return "\n".join(lines) + "\n", acceptance


# ---------------------------------------------------------------------------
# 运行部分


def _q(ir) -> np.ndarray:
    return np.array([a.partial_charge for a in ir.atoms])


def main() -> int:
    parser = argparse.ArgumentParser(prog="gate_mvp3")
    parser.add_argument("--scan-cache", type=Path, required=True,
                        help="本次 G2 运行目录下的 scans/（只读；必须显式给出）")
    parser.add_argument("--skip-compare", action="store_true", help="调试用：跳过跨引擎比较（G3.3 判为未通过）")
    parser.add_argument("--skip-pytest", action="store_true", help="调试用：不跑 pytest（G3.5 判为未通过）")
    parser.add_argument("--molecules", nargs="+", default=list(G31_NAMES), choices=G31_NAMES,
                        help="调试用：只跑部分分子（正式 Gate 必须全部）（子集运行 passed 恒为 false，决策 0019）")
    args = parser.parse_args()
    if not args.scan_cache.is_dir():
        parser.error(f"--scan-cache {args.scan_cache}：不是目录")

    from gate_mvp1 import _pytest_summary
    from gate_mvp2 import _load, _mm_profile, _shifted_rmse, judge_compare

    from ffcompiler.charge.am1bcc import AM1BCCCharges
    from ffcompiler.charge.apply import assign_charges
    from ffcompiler.charge.loaded import LoadedCharges, read_mol2_atoms
    from ffcompiler.charge.nagl import NAGLCharges
    from ffcompiler.runtime import RunContext
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider, _no_reference, model_file_identity
    from ffcompiler.torsion.rotors import find_rotors

    with RunContext("gate_mvp3", root=REPO_ROOT / "jobs", seed=SEED, inputs=[
        FIXTURE_DIR / n / f"{n}.{ext}" for n in G31_NAMES for ext in ("prmtop", "inpcrd", "sdf", "mol2")
    ]) as ctx:
        run_dir, log = ctx.dir, ctx.log
        model_path = default_model_path()
        ctx.manifest["scan_cache"] = str(args.scan_cache)
        ctx.manifest["model_file"] = model_file_identity(model_path)
        names31 = [n for n in G31_NAMES if n in args.molecules]
        names = [n for n in NAMES if n in args.molecules]
        systems = {n: _load(n) for n in names31}
        failures: list[dict] = []

        # G3.1：AM1-BCC 与固件 mol2
        g31: dict[str, float | None] = {}
        am1bcc_ir = {}
        for n in names31:
            ir, pos, mol = systems[n]
            try:
                new_ir, _ = assign_charges(ir, AM1BCCCharges(run_dir / "am1bcc" / n), pos, mol)
                ref = np.array([q for _, q in read_mol2_atoms(FIXTURE_DIR / n / f"{n}.mol2")])
                g31[n] = float(np.max(np.abs(_q(new_ir) - ref)))
                am1bcc_ir[n] = new_ir
            except Exception as exc:  # noqa: BLE001
                g31[n] = None
                failures.append({"stage": "g31", "molecule": n, "reason": repr(exc)})
            log.info("G3.1 %s: max|Δq| = %s", n, g31[n])

        # G3.2 / G3.3：三种 provider × 11 分子
        g32: dict = {}
        compare_result: dict = {}
        sander_notes: dict[str, list[str]] = {}
        nagl_vs: dict[str, dict] = {}
        from ffcompiler.validation.engines import compare
        for n in names:
            ir, pos, mol = systems[n]
            for method in PROVIDERS:
                key = (method, n)
                provider = {"loaded": lambda: LoadedCharges(),
                            "am1bcc": lambda: AM1BCCCharges(run_dir / "am1bcc_g32" / n),
                            "nagl": lambda: NAGLCharges()}[method]()
                try:
                    new_ir, _ = assign_charges(ir, provider, pos, mol)
                except Exception as exc:  # noqa: BLE001
                    g32[key] = compare_result[key] = None
                    failures.append({"stage": f"g32/{method}", "molecule": n, "reason": repr(exc)})
                    continue
                g32[key] = float(abs(_q(new_ir).sum() - new_ir.total_charge))
                if method == "nagl" and n in am1bcc_ir:
                    d = _q(new_ir) - _q(am1bcc_ir[n])
                    nagl_vs[n] = {"max_abs": float(np.max(np.abs(d))), "rms": float(np.sqrt(np.mean(d * d)))}
                if args.skip_compare:
                    compare_result[key] = None
                    continue
                workdir = run_dir / "compare" / method / n
                try:
                    xreport = compare(new_ir, pos, workdir).to_json()
                    (workdir / "report.json").write_text(json.dumps(xreport, indent=2, ensure_ascii=False) + "\n")
                    compare_result[key], notes = judge_compare(xreport)
                    if notes:
                        sander_notes[_key(key)] = notes
                except Exception as exc:  # noqa: BLE001
                    compare_result[key] = False
                    failures.append({"stage": f"compare/{method}", "molecule": n, "reason": repr(exc)})
                log.info("G3.2/3.3 %s/%s: |Σq−Q| = %.3g, compare = %s", method, n, g32[key], compare_result[key])

        # G3.4：二面角重拟合（G2 扫描缓存，cache_only）
        g34: dict[str, dict] = {}
        expected_g34 = {n: c for n, c in ROTOR_COUNTS.items() if n in names}
        profiles: list[dict] = []
        for n, count in expected_g34.items():
            ir, pos, mol = systems[n]
            try:
                if len(find_rotors(ir, mol, pos)) != count:
                    raise RuntimeError(f"find_rotors 给出 {len(find_rotors(ir, mol, pos))} 根，预期 {count}")
                q_ir, _ = assign_charges(ir, AM1BCCCharges(run_dir / "am1bcc_g34a" / n), pos, mol)
                provider = MLPTorsionProvider(_no_reference, supported=frozenset(range(1, 119)),
                                              scan_cache=args.scan_cache / n, model_path=model_path,
                                              cache_only=True)
                mlp_ir, derive_report = provider.parameterize(q_ir, pos, mol)
                if [r["status"] for r in derive_report.rotors] != ["mlp"] * count:
                    raise RuntimeError(f"derive from cache: {[r['reason'] for r in derive_report.rotors]}")
                again, _ = assign_charges(mlp_ir, AM1BCCCharges(run_dir / "am1bcc_g34b" / n), pos, mol,
                                          scan_cache=args.scan_cache / n)
                nagl_ir, nagl_report = assign_charges(mlp_ir, NAGLCharges(), pos, mol,
                                                      scan_cache=args.scan_cache / n)
                statuses = [r["status"] for r in nagl_report.torsion_refit["rotors"]]
                g34[n] = {"identical": again.proper_torsions == mlp_ir.proper_torsions,
                          "nagl_statuses": statuses}
                for bond, scan in derive_report.scans.items():
                    label = f"{bond[0]}-{bond[1]}"
                    for charges, cur in (("am1bcc", mlp_ir), ("nagl", nagl_ir)):
                        profiles.append({"molecule": n, "bond": label, "charges": charges,
                                         "profile_rmse": _shifted_rmse(_mm_profile(cur, scan.positions), scan.energy)})
            except Exception as exc:  # noqa: BLE001
                failures.append({"stage": "g34", "molecule": n, "reason": repr(exc)})
            log.info("G3.4 %s: %s", n, g34.get(n))

        pytest_summary = {} if args.skip_pytest else _pytest_summary(run_dir)[0]
        report_md, acceptance = build_report(
            full_matrix=set(args.molecules) == set(G31_NAMES),
            g31=g31, g32=g32, compare=compare_result, sander_notes=sander_notes, g34=g34,
            nagl_vs_am1bcc=nagl_vs, profiles=profiles, pytest_summary=pytest_summary, failures=failures,
            expected_g31=names31, expected_compare=[(m, n) for n in names for m in PROVIDERS],
            expected_g34=expected_g34,
        )
        (run_dir / "report.md").write_text(report_md)
        (run_dir / "acceptance.json").write_text(json.dumps(acceptance, indent=2, ensure_ascii=False) + "\n")
        ctx.manifest["gate"] = "G3"
        ctx.manifest["debug_flags"] = {"skip_compare": args.skip_compare, "skip_pytest": args.skip_pytest,
                                       "molecules": list(args.molecules)}
        print(f"G3 passed={acceptance['passed']} failures={len(failures)}")
        return 0 if acceptance["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
