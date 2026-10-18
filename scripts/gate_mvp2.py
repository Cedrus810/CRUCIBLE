#!/usr/bin/env python
"""Gate G2：MLP（MACE-POLAR-1-M）二面角评估与验收（MVP-2 计划 T30）。

在 RunContext("gate_mvp2") 下运行，写 report.md、metrics.csv、acceptance.json，
扫描剖面存 scans/<分子>/<j>-<k>.npz（供以后 torsion-profile 回归比对）。

    python scripts/gate_mvp2.py [--jobs N] [--threads T] [--device cpu] [--resume RUN_DIR]
                                [--skip-compare] [--skip-pytest] [--molecules ...]   # 后三个只用于调试

并行单位 = (分子, 可转键)：先用进程池把全部扫描写进缓存，再逐分子用 MLPTorsionProvider
从缓存拟合、写回 IR。--resume 指向旧运行目录时，先把其 scans/ 复制进本次运行目录再复用（旧目录只读，扫描参数逐项相等才复用）。
拟合 RMSE、剖面复现、滞后、键长变化只报告、不设门槛；passed 由 G2.1–G2.5 逐条推出
（fail-closed）。--skip-compare 只用于调试：跳过的跨引擎比较不算通过。
G2.3 以 OpenMM 为准：sander 腿照常运行，超差只写进报告的已知局限，不计失败（决策 0017）。
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import numpy as np  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "dmso", "nma", "methyl_benzoate", "acetate", "methylammonium",
    "benzene", "chlorobenzene", "cyclopropane", "cyclobutane", "norbornane",
)  # 乙腈（近线性角）不进 Gate（决策 0013）；它没有可转键由 RF1 测试覆盖
SEED = 20261004

SELF_CONSISTENCY_TESTS = {
    "G2.1 fourier": "tests/unit/test_torsion_fourier.py::test_recovers_known_coefficients_on_nonuniform_grid",
    "G2.1 mm": "tests/unit/test_torsion_mlp.py::test_self_consistency_single_quartet_bond",
}

# G2.4：计划 §4 Review Focus 1–5 对应的测试（RF6 = G2.3 跨引擎比较本身）
REVIEW_FOCUS_TESTS: dict[str, tuple[str, ...]] = {
    "RF1 不该扫的键": tuple(
        f"tests/unit/test_torsion_rotors.py::test_rotor_counts[{n}]"
        for n in ("acetonitrile", "benzene", "chlorobenzene", "cyclobutane", "cyclopropane", "norbornane")
    ),
    "RF2 带电分子": (
        "tests/unit/test_torsion_mlp.py::test_charged_molecule_is_scanned_with_its_charge",
        "tests/unit/test_torsion_mlp.py::test_odd_electrons_falls_back_whole_molecule",
        "tests/unit/test_torsion_scan.py::test_info_is_passed_to_calculator",
        "tests/unit/test_torsion_mace.py::test_acetate_charges_sum_to_minus_one",
    ),
    "RF3 元素不支持": (
        "tests/unit/test_torsion_mlp.py::test_unsupported_element_falls_back_whole_molecule",
    ),
    "RF4 SDF/IR 不一致": (
        "tests/unit/test_torsion_rotors.py::test_atomic_number_mismatch_raises",
        "tests/unit/test_torsion_rotors.py::test_missing_bond_raises",
    ),
    "RF5 对称等价键": (
        "tests/unit/test_torsion_rotors.py::test_dmso_two_rotors_same_class_nma_methyls_distinct",
        "tests/unit/test_torsion_mlp.py::test_dmso_both_rotors_parameterized",
    ),
}

CSV_HEADER = (
    "molecule,bond,rotor_class,status,reason,fit_rmse,profile_rmse_fit,profile_rmse_gaff2,"
    "max_hysteresis,max_bond_change_pm"
)

REPORT_HEADER = (
    "# Gate G2 评估报告（MVP-2，MACE-POLAR-1-M 二面角）\n\n"
    "拟合与剖面指标只报告、不设门槛。\n\n"
    "| 指标 | 说明 |\n|---|---|\n"
    "| fit_rmse | Fourier 拟合 vs 目标 y = E_MLP − E_MM[IR \\ 键]（kJ/mol） |\n"
    "| profile_rmse_fit | 扫描几何上，新 IR 的 MM 能量 vs E_MLP，各减自身最小值后的 RMSE（kJ/mol） |\n"
    "| profile_rmse_gaff2 | 同上，原 GAFF2 IR vs E_MLP（kJ/mol） |\n"
    "| max_hysteresis | 正反两遍同一网格点的最大能量差（kJ/mol） |\n"
    "| max_bond_change_pm | 扫描中相对起点的最大键长变化（pm） |\n"
)


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _is_nan(value) -> bool:
    return isinstance(value, float) and math.isnan(value)


def judge_compare(report: dict) -> tuple[bool, list[str]]:
    """跨引擎比较结果（`CrossEngineReport.to_json()`）→ (G2.3 是否通过, sander 超差说明)。

    决策 0017：sander 腿（sander-*）只报告；其余腿（numpy、OpenMM、GROMACS）全部通过、
    没有出错，且至少有一行非 sander 的能量比较，才算通过。
    """
    def is_sander(leg: str) -> bool:
        return leg.startswith("sander")

    rows = [r for r in report["rows"] if not is_sander(r["leg"])]
    forces = [r for r in report["force_rows"] if not is_sander(r["leg"])]
    errors = [leg for leg in report["errors"] if not is_sander(leg)]
    passed = bool(rows) and not errors and all(r["passed"] for r in rows + forces)
    notes = [f"{r['leg']} {r['category']}: delta {r['delta']:.3g} > tol {r['tol']:.3g}"
             for r in report["rows"] if is_sander(r["leg"]) and not r["passed"]]
    notes += [f"{leg} error: {msg}" for leg, msg in report["errors"].items() if is_sander(leg)]
    return passed, notes


def build_report(
    rows: list[dict],
    failures: list[dict],
    *,
    full_matrix: bool,
    pytest_summary: dict,
    review_focus: dict[str, list[tuple[str, str]]],
    self_consistency: dict[str, str],
    expected: list[tuple[str, str]],
    compare: dict[str, bool | None],
    sander_notes: dict[str, list[str]] | None = None,
) -> tuple[str, str, dict]:
    """纯函数：指标行 → (report.md, metrics.csv, acceptance.json 内容)。

    failures 不丢弃；NaN 指标、缺失的 (分子, 键) 会补进 failures。兜底（status=fallback）
    不是失败，单列。passed = G2.1–G2.5 全部为真。sander_notes 只报告（决策 0017）。
    """
    failures = list(failures)
    sander_notes = {m: list(n) for m, n in (sander_notes or {}).items() if n}
    for row in rows:
        nan_keys = [k for k, v in row.items() if _is_nan(v)]
        if nan_keys:
            failures.append({"molecule": row["molecule"], "bond": row["bond"],
                             "reason": f"指标为 NaN: {', '.join(nan_keys)}"})
    covered = {(r["molecule"], r["bond"]) for r in rows} | {(f.get("molecule"), f.get("bond")) for f in failures}
    for molecule, bond in expected:
        if (molecule, bond) not in covered:
            failures.append({"molecule": molecule, "bond": bond, "reason": "可转键既无结果也无失败记录"})
    fallbacks = [{"molecule": r["molecule"], "bond": r["bond"], "reason": r["reason"]}
                 for r in rows if r.get("status") == "fallback"]

    csv_lines = [CSV_HEADER] + [
        ",".join(_fmt(row.get(key)) for key in CSV_HEADER.split(",")) for row in rows
    ]

    lines = [REPORT_HEADER]
    if not full_matrix:
        lines.append("\n> 子集运行（scope=subset）：判据结果只供调试，不是正式通过（决策 0019）。\n")
    lines += ["\n## 每根可转键\n",
             "| 分子 | 键 | 类 | 状态 | 拟合 RMSE | 剖面 RMSE（新） | 剖面 RMSE（GAFF2） | 滞后 | 键长变化 pm |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        status = r["status"] if r["status"] == "mlp" else f"{r['status']}（{r['reason']}）"
        lines.append(f"| {r['molecule']} | {r['bond']} | {r.get('rotor_class')} | {status} | "
                     f"{_fmt(r.get('fit_rmse'))} | {_fmt(r.get('profile_rmse_fit'))} | "
                     f"{_fmt(r.get('profile_rmse_gaff2'))} | {_fmt(r.get('max_hysteresis'))} | "
                     f"{_fmt(r.get('max_bond_change_pm'))} |")
    rotor_counts: dict[str, int] = {}
    for molecule, _ in expected:
        rotor_counts[molecule] = rotor_counts.get(molecule, 0) + 1
    for molecule in compare:
        rotor_counts.setdefault(molecule, 0)

    lines += ["\n## 跨引擎比较（G2.3）\n", "| 分子 | 可转键数 | compare |", "|---|---|---|"]
    for molecule in sorted(compare):
        lines.append(f"| {molecule} | {rotor_counts.get(molecule, 0)} | {compare[molecule]} |")
    lines.append("\n### 已知局限：sander 腿只报告、不计失败（决策 0017）\n")
    lines += [f"- {m}: {note}" for m in sorted(sander_notes) for note in sander_notes[m]] or ["- 无"]

    lines += ["\n## 自洽测试（G2.1）\n", "| 判据 | 结果 |", "|---|---|"]
    for label, outcome in self_consistency.items():
        lines.append(f"| {label} | {outcome} |")
    lines += ["\n## Review Focus 测试（G2.4；RF6 即 G2.3）\n", "| Review Focus | 测试 | 结果 |", "|---|---|---|"]
    for label, tests in review_focus.items():
        for test_id, outcome in tests:
            lines.append(f"| {label} | `{test_id}` | {outcome} |")

    criteria = {
        "G2.1": bool(self_consistency) and all(v == "passed" for v in self_consistency.values()),
        "G2.2": not failures,
        "G2.3": bool(compare) and all(v is True for v in compare.values()),
        "G2.4": bool(review_focus) and all(
            outcome == "passed" for tests in review_focus.values() for _, outcome in tests),
        "G2.5": all(pytest_summary.get(m, {}).get("accepted") is True for m in ("not external", "external")),
    }
    lines += ["\n## 判据\n", "| 判据 | 结果 |", "|---|---|"]
    lines += [f"| {k} | {'通过' if ok else '未通过'} |" for k, ok in criteria.items()]
    lines.append("\n## 兜底（GAFF2，不计入失败）\n")
    lines += [f"- {f['molecule']} {f['bond']}: {f['reason']}" for f in fallbacks] or ["- 无"]
    lines.append("\n## 失败项\n")
    lines += [f"- {f.get('molecule')} {f.get('bond')}: {f.get('reason')}" for f in failures] or ["- 无"]

    acceptance = {
        "gate": "G2",
        "scope": "full" if full_matrix else "subset",
        "passed": full_matrix and all(criteria.values()),
        "criteria": criteria,
        "rotor_counts": rotor_counts,
        "compare": compare,
        "sander_known_limitation": sander_notes,
        "self_consistency": self_consistency,
        "fallbacks": fallbacks,
        "failures": failures,
        "pytest_summary": pytest_summary,
    }
    return "\n".join(lines) + "\n", "\n".join(csv_lines) + "\n", acceptance


# ---------------------------------------------------------------------------
# 运行部分


def _load(name: str):
    from rdkit import Chem

    from ffcompiler.importers import from_amber

    mdir = FIXTURE_DIR / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


def _init_worker(threads: int) -> None:
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[var] = str(threads)


def _scan_job(name: str, bond: tuple[int, int], cache_root: str, scan_kwargs: dict,
              threads: int, device: str, model_path: str) -> dict:
    """子进程：扫一根可转键写进缓存。已有且键逐项相等就跳过。"""
    from ffcompiler.torsion.mace import MODEL_NAME, molecule_info, polar_calculator
    from ffcompiler.torsion.mlp import load_scan, model_file_identity, save_scan, scan_cache_key
    from ffcompiler.torsion.rotors import find_rotors
    from ffcompiler.torsion.scan import ScanError, relaxed_scan

    t0 = time.monotonic()
    ir, pos, mol = _load(name)
    rotor = {r.bond: r for r in find_rotors(ir, mol, pos)}[tuple(bond)]
    info = molecule_info(ir)
    numbers = [a.atomic_number for a in ir.atoms]
    bonds = [(b.i, b.j) for b in ir.bonds]
    key = scan_cache_key(MODEL_NAME, rotor.scan_quartet, info, scan_kwargs,
                         model_file=model_file_identity(Path(model_path)), numbers=numbers, bonds=bonds)
    path = Path(cache_root) / name / f"{bond[0]}-{bond[1]}.npz"
    if load_scan(path, key) is not None:
        return {"molecule": name, "bond": bond, "status": "cached", "seconds": 0.0}
    try:
        scan = relaxed_scan(polar_calculator(model_path=Path(model_path), device=device, threads=threads),
                            numbers, pos, rotor.scan_quartet, bonds, info=info, **scan_kwargs)
    except ScanError as exc:
        return {"molecule": name, "bond": bond, "status": "scan_error", "error": str(exc),
                "seconds": time.monotonic() - t0}
    save_scan(path, key, scan)
    return {"molecule": name, "bond": bond, "status": "scanned", "seconds": time.monotonic() - t0}


def _mm_profile(ir, positions_stack) -> np.ndarray:
    from openmm import unit

    from ffcompiler.exporters.openmm import to_openmm_system
    from ffcompiler.validation.engines import reference_context

    ctx = reference_context(to_openmm_system(ir))
    out = []
    for x in positions_stack:
        ctx.setPositions(x)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    return np.asarray(out)


def _shifted_rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(((a - a.min()) - (b - b.min())) ** 2)))


def seed_scan_cache(old_root: Path, new_root: Path, names: Iterable[str]) -> list[str]:
    """--resume：把旧运行目录 scans/<分子>/*.npz 复制进本次运行的 scans/；旧目录只读。
    键不一致的文件照常由 _scan_job 判为未命中并重扫，新结果只写本次目录。"""
    copied: list[str] = []
    for name in names:
        src_dir = Path(old_root) / name
        if not src_dir.is_dir():
            continue
        dst_dir = Path(new_root) / name
        dst_dir.mkdir(parents=True, exist_ok=True)
        for src in sorted(src_dir.glob("*.npz")):
            if src.name.endswith(".tmp.npz"):
                continue
            shutil.copy2(src, dst_dir / src.name)
            copied.append(f"{name}/{src.name}")
    return copied


def main() -> int:
    parser = argparse.ArgumentParser(prog="gate_mvp2")
    parser.add_argument("--jobs", type=int, default=1, help="并行扫描进程数")
    parser.add_argument("--threads", type=int, default=4, help="每个扫描进程的线程数")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--resume", type=Path, default=None, help="旧运行目录：复制其 scans/ 到本次目录后复用（旧目录不写）")
    parser.add_argument("--skip-compare", action="store_true",
                        help="调试用：跳过跨引擎比较（Gate 将判为未通过）")
    parser.add_argument("--skip-pytest", action="store_true",
                        help="调试用：不跑 pytest（G2.1/G2.4/G2.5 判为未通过）")
    parser.add_argument("--molecules", nargs="+", default=list(NAMES), choices=NAMES,
                        help="调试用：只跑部分分子（正式 Gate 必须全部）（子集运行 passed 恒为 false，决策 0019）")
    args = parser.parse_args()
    if args.resume is not None and not (args.resume / "scans").is_dir():
        parser.error(f"--resume {args.resume}：没有 scans/ 目录（路径写错？要全新扫描就去掉 --resume）")
    _init_worker(args.threads)

    from gate_mvp1 import _pytest_summary, junit_outcomes  # 复用 G1 的 pytest/JUnit 汇总

    from ffcompiler.runtime import RunContext
    from ffcompiler.torsion.mace import default_model_path, supported_atomic_numbers
    from ffcompiler.torsion.mlp import (SCAN_DEFAULTS, MLPTorsionProvider, model_file_identity,
                                        molecule_fallback_reason, package_versions)
    from ffcompiler.torsion.rotors import find_rotors

    with RunContext("gate_mvp2", root=REPO_ROOT / "jobs", seed=SEED, inputs=[
        FIXTURE_DIR / n / f"{n}.{ext}" for n in NAMES for ext in ("prmtop", "inpcrd", "sdf")
    ]) as ctx:
        run_dir = ctx.dir
        names = tuple(args.molecules)
        cache_root = run_dir / "scans"
        ctx.manifest["scan_cache"] = str(cache_root)
        log = ctx.log
        if args.resume:
            seeded = seed_scan_cache(args.resume / "scans", cache_root, names)
            ctx.manifest["scan_cache_seeded_from"] = {"dir": str(args.resume / "scans"), "files": seeded}
            if not seeded:
                log.warning("--resume %s/scans 里没有所选分子的扫描，全部重扫", args.resume)
        model_path = default_model_path()
        ctx.manifest["model_file"] = model_file_identity(model_path)
        ctx.manifest["mlp_package_versions"] = package_versions()
        supported = supported_atomic_numbers(model_path)
        systems = {n: _load(n) for n in names}
        rotors = {n: find_rotors(systems[n][0], systems[n][2], systems[n][1]) for n in names}
        expected = [(n, f"{r.bond[0]}-{r.bond[1]}") for n in names for r in rotors[n]]
        jobs = [(n, r.bond) for n in names for r in rotors[n]
                if molecule_fallback_reason(systems[n][0], supported) is None]
        log.info("可转键 %d 根，需扫描 %d 根；进程 %d × 线程 %d；缓存 %s",
                 len(expected), len(jobs), args.jobs, args.threads, cache_root)

        failures: list[dict] = []
        with ProcessPoolExecutor(max_workers=max(1, args.jobs), mp_context=mp.get_context("spawn"),
                                 initializer=_init_worker, initargs=(args.threads,)) as pool:
            futures = {job: pool.submit(_scan_job, job[0], job[1], str(cache_root), dict(SCAN_DEFAULTS),
                                        args.threads, args.device, str(model_path)) for job in jobs}
            scan_errors: dict[tuple[str, str], str] = {}
            for (name, bond), fut in futures.items():
                try:
                    res = fut.result()
                    log.info("扫描 %s %s：%s（%.0f s）%s", name, bond, res["status"], res["seconds"],
                             res.get("error", ""))
                    if res["status"] == "scan_error":
                        scan_errors[(name, f"{bond[0]}-{bond[1]}")] = res["error"]
                except Exception as exc:  # noqa: BLE001  子进程崩溃：记失败，不丢
                    failures.append({"molecule": name, "bond": f"{bond[0]}-{bond[1]}",
                                     "reason": f"扫描进程异常: {exc!r}"})
                    log.error("扫描 %s %s 进程异常：%r", name, bond, exc)

        rows: list[dict] = []
        compare_result: dict[str, bool | None] = {}
        sander_notes: dict[str, list[str]] = {}
        for name in names:
            ir, pos, mol = systems[name]
            try:
                # 拟合阶段只读缓存：子进程没扫成的键记为兜底，不在主进程串行重扫（评审 Important 2）
                provider = MLPTorsionProvider(lambda ir_: None, supported=supported,
                                              scan_cache=cache_root / name, model_path=model_path,
                                              cache_only=True)
                new_ir, report = provider.parameterize(ir, pos, mol)
            except Exception as exc:  # noqa: BLE001
                failures.append({"molecule": name, "bond": "*", "reason": f"参数化失败: {exc!r}"})
                compare_result[name] = None
                continue
            from ffcompiler.ir import io as ir_io

            (run_dir / "ir").mkdir(exist_ok=True)
            ir_io.to_json(new_ir, run_dir / "ir" / f"{name}.json")
            for row in report.rotors:
                bond = tuple(row["bond"])
                bond_label = f"{bond[0]}-{bond[1]}"
                reason = row["reason"] or ""
                if reason == "no_cached_scan" and (name, bond_label) in scan_errors:
                    reason = "mlp_error"  # 子进程扫描时参考势出错（错误文本见 run.log）
                out = {"molecule": name, "bond": bond_label, "rotor_class": row["rotor_class"],
                       "status": row["status"], "reason": reason,
                       "fit_rmse": row.get("fit_rmse"), "max_hysteresis": row.get("max_hysteresis"),
                       "max_bond_change_pm": row.get("max_bond_change_pm")}
                scan = report.scans.get(bond)
                if scan is not None:
                    e_gaff2 = _mm_profile(ir, scan.positions)
                    out["profile_rmse_gaff2"] = _shifted_rmse(e_gaff2, scan.energy)
                    if row["status"] == "mlp":
                        out["profile_rmse_fit"] = _shifted_rmse(_mm_profile(new_ir, scan.positions), scan.energy)
                rows.append(out)
            if args.skip_compare:
                compare_result[name] = None
            else:
                from ffcompiler.validation.engines import compare

                try:
                    xreport = compare(new_ir, pos, run_dir / "compare" / name).to_json()
                    (run_dir / "compare" / name / "report.json").write_text(
                        json.dumps(xreport, indent=2, ensure_ascii=False) + "\n")
                    compare_result[name], sander_notes[name] = judge_compare(xreport)
                except Exception as exc:  # noqa: BLE001
                    compare_result[name] = False
                    failures.append({"molecule": name, "bond": "*", "reason": f"跨引擎比较出错: {exc!r}"})

        if args.skip_pytest:
            pytest_summary, outcomes = {}, {}
        else:
            pytest_summary, xmls = _pytest_summary(run_dir)
            outcomes = junit_outcomes(xmls)
        review_focus = {label: [(t, outcomes.get(t, "missing")) for t in tests]
                        for label, tests in REVIEW_FOCUS_TESTS.items()}
        self_consistency = {label: outcomes.get(t, "missing") for label, t in SELF_CONSISTENCY_TESTS.items()}
        report_md, metrics_csv, acceptance = build_report(
            rows, failures, full_matrix=set(names) == set(NAMES),
            pytest_summary=pytest_summary, review_focus=review_focus,
            self_consistency=self_consistency, expected=expected, compare=compare_result,
            sander_notes=sander_notes,
        )
        (run_dir / "report.md").write_text(report_md)
        (run_dir / "metrics.csv").write_text(metrics_csv)
        (run_dir / "acceptance.json").write_text(json.dumps(acceptance, indent=2, ensure_ascii=False) + "\n")
        ctx.manifest["gate"] = "G2"
        ctx.manifest["debug_flags"] = {"skip_compare": args.skip_compare, "skip_pytest": args.skip_pytest,
                                       "molecules": list(names)}
        print(f"G2 passed={acceptance['passed']} failures={len(acceptance['failures'])}")
        return 0 if acceptance["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
