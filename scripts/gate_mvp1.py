#!/usr/bin/env python
"""Gate G1：Hessian 方法评估矩阵与验收（MVP-1 计划 T23）。

在 RunContext("gate_mvp1") 下运行，写 report.md、metrics.csv、acceptance.json。
评估矩阵：11 个固件（乙腈只作诊断，决策 0013） × Hessian provider（gfnff/gfn2）× 方法
（gaff2 基线、mseminario、m2seminario、drih、global_fit）。

    python scripts/gate_mvp1.py [--providers gfnff gfn2] [--jobs N] [--skip-compare]

频率/几何指标只报告、不设门槛（设计 §24 Contribution B 是研究问题）；
验收判据见计划 §6。`passed` 由 G1.1–G1.7 逐条推出（fail-closed，见
acceptance.json 的 `criteria`）；G1.8（理论文档审阅）是人工判据，不进自动
判定，报告里单列。--skip-compare 只用于调试：跳过的跨引擎比较不算通过。
报告生成函数 build_report 是纯函数，单测用假指标行驱动（test_gate_mvp1_report.py）。
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import numpy as np  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "scripts"))
from sobtop_oracle import GATE_CASES, METHODS as ORACLE_METHODS  # noqa: E402
from sobtop_oracle import REFS, load_case, oracle_mismatches  # noqa: E402

from ffcompiler.bonded.assign import NonPositiveForceConstantError  # noqa: E402
from ffcompiler.hessian.base import check_hessian, harmonic_frequencies  # noqa: E402
from ffcompiler.hessian.xtb import run_xtb_hessian  # noqa: E402
from ffcompiler.importers import from_amber  # noqa: E402
from ffcompiler.runtime import RunContext  # noqa: E402
from ffcompiler.validation.pytest_results import summarize_pytest  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
)  # 乙腈（近线性角）不进 Gate，只作诊断（决策 0013）
METHODS = ("gaff2", "mseminario", "m2seminario", "drih", "global_fit")
SEED = 20261004
SELF_CONSISTENCY_TOL = 1e-5  # G1.3

# G1.6：计划 §4 Review Focus 1–5 对应的测试（JUnit 里必须存在且通过）
REVIEW_FOCUS_TESTS: dict[str, tuple[str, ...]] = {
    "RF1 原子顺序/拓扑": (
        "tests/unit/test_assign.py::test_shuffled_atom_order_rejected",
        "tests/unit/test_assign.py::test_stretched_bond_rejected",
        "tests/unit/test_global_fit.py::test_topology_checks",
    ),
    "RF2 线性键角": (
        "tests/unit/test_internal.py::test_acetonitrile_linear_angle",
        "tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[mseminario]",
        "tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[m2seminario]",
        "tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[drih]",
        "tests/unit/test_drih.py::test_linear_angle_in_input_raises",
    ),
    "RF3 非驻点": (
        "tests/unit/test_hessian_base.py::test_nonstationary_gradient",
        "tests/unit/test_drih.py::test_nonstationary_rejected",
        "tests/unit/test_global_fit.py::test_nonstationary_rejected",
    ),
    "RF4 鞍点/虚频": (
        "tests/unit/test_hessian_base.py::test_imaginary_frequency_error_and_warning",
    ),
    "RF5 带电分子": (
        "tests/unit/test_hessian_xtb.py::test_run_charged_molecules[acetate--1]",
        "tests/unit/test_hessian_xtb.py::test_run_charged_molecules[methylammonium-1]",
    ),
}

CSV_HEADER = (
    "molecule,provider,method,freq_rmse_cm,freq_mae_cm,freq_max_cm,"
    "freq_rmse_lt1000,freq_rmse_ge1000,kabsch_rmsd_pm,bond_mae_pm,"
    "angle_mae_deg,max_class_spread,gf_alpha,gf_lam,gf_cond,"
    "gf_rel_res_fit,gf_rel_res_prior,gf_n_at_bound,gf_max_force,compare_passed"
)

REPORT_HEADER = (
    "# Gate G1 评估报告（MVP-1）\n\n"
    "频率与几何指标只报告、不设门槛。\n\n"
    "| 指标 | 说明 |\n|---|---|\n"
    "| freq_rmse_cm | FF Hessian 频率 vs 参考 Hessian 频率的 RMSE（cm⁻¹） |\n"
    "| kabsch_rmsd_pm | FF 优化几何 vs 参考几何的 Kabsch RMSD（pm） |\n"
    "| bond_mae_pm | 优化后键长 MAE（pm） |\n"
    "| angle_mae_deg | 优化后键角 MAE（deg） |\n"
)


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _is_nan(value) -> bool:
    return isinstance(value, float) and math.isnan(value)


def _full_match(match: str | None) -> bool:
    if not match or "/" not in match:
        return False
    n_ok, n_all = (int(x) for x in match.split("/"))
    return n_all > 0 and n_ok == n_all


def build_report(
    rows: list[dict],
    failures: list[dict],
    *,
    pytest_summary: dict,
    sobtop_match: dict[str, str],
    m2_naming: str,
    global_fit_self_consistency: float | None,
    review_focus: dict[str, list[tuple[str, str]]],
    expected: list[tuple[str, str, str]] | None = None,
) -> tuple[str, str, dict]:
    """纯函数：指标行 → (report.md, metrics.csv, acceptance.json 内容)。

    failures 不会被丢弃；NaN 指标、缺失的 (分子, provider, 方法) 组合会被补进
    failures。passed = G1.1–G1.7 全部为真。

    `category == "method_limitation"` 且 `evidence.sobtop_confirmed` 为真的条目
    （方法给出非正 k，Sobtop 在同一 Hessian 上给出同样结果，决策 0014）单列为
    已知方法局限，照样写进报告，但不计入 G1.5；没有证据的仍是失败。
    """
    limitations = [
        f for f in failures
        if f.get("category") == "method_limitation"
        and (f.get("evidence") or {}).get("sobtop_confirmed") is True
    ]
    failures = [f for f in failures if f not in limitations]
    for row in rows:
        nan_keys = [k for k, v in row.items() if _is_nan(v)]
        if nan_keys:
            failures.append({**{k: row[k] for k in ("molecule", "provider", "method")},
                             "reason": f"指标为 NaN: {', '.join(nan_keys)}"})
    if expected is not None:
        covered = {(r["molecule"], r["provider"], r["method"]) for r in rows}
        covered |= {(f.get("molecule"), f.get("provider"), f.get("method"))
                    for f in failures + limitations}
        for combo in expected:
            if combo not in covered:
                failures.append({"molecule": combo[0], "provider": combo[1],
                                 "method": combo[2], "reason": "组合既无结果也无失败记录"})

    csv_lines = [CSV_HEADER]
    for row in rows:
        csv_lines.append(",".join(_fmt(row.get(key)) for key in CSV_HEADER.split(",")))

    providers = sorted({row["provider"] for row in rows})
    report_lines = [REPORT_HEADER]
    metric_keys = [
        ("freq_rmse_cm", "freq RMSE (cm⁻¹)"),
        ("kabsch_rmsd_pm", "Kabsch RMSD (pm)"),
        ("bond_mae_pm", "bond MAE (pm)"),
        ("angle_mae_deg", "angle MAE (deg)"),
        ("max_class_spread", "max class spread"),
    ]
    for provider in providers:
        report_lines.append(f"\n## Hessian: {provider}\n")
        report_lines.append("| method | 指标 | 中位数 | 最大值 |")
        report_lines.append("|---|---|---|---|")
        for method in METHODS:
            subset = [
                r for r in rows
                if r["provider"] == provider and r["method"] == method
            ]
            for key, label in metric_keys:
                values = [r[key] for r in subset
                          if isinstance(r.get(key), (int, float)) and not _is_nan(r[key])]
                if not values:
                    continue
                report_lines.append(
                    f"| {method} | {label} | {np.median(values):.4g} | {np.max(values):.4g} |"
                )

    report_lines.append("\n## Sobtop 比对\n")
    report_lines.append("| 方法 | 一致 case 数 / 总数 |")
    report_lines.append("|---|---|")
    for method, match in sorted(sobtop_match.items()):
        report_lines.append(f"| {method} | {match} |")

    report_lines.append("\n## Review Focus 测试（G1.6）\n")
    report_lines.append("| Review Focus | 测试 | 结果 |")
    report_lines.append("|---|---|---|")
    for label, tests in review_focus.items():
        for test_id, outcome in tests:
            report_lines.append(f"| {label} | `{test_id}` | {outcome} |")

    cross_engine = bool(rows) and all(row.get("compare_passed") is True for row in rows)
    gf = global_fit_self_consistency
    criteria = {
        "G1.1": _full_match(sobtop_match.get("mseminario")) and _full_match(sobtop_match.get("drih")),
        "G1.2": (m2_naming == "m2seminario" and _full_match(sobtop_match.get("m2seminario")))
        or m2_naming == "symseminario",
        "G1.3": gf is not None and not _is_nan(gf) and gf <= SELF_CONSISTENCY_TOL,
        "G1.4": cross_engine,
        "G1.5": not failures,
        "G1.6": bool(review_focus) and all(
            outcome == "passed" for tests in review_focus.values() for _, outcome in tests
        ),
        "G1.7": all(
            pytest_summary.get(marker, {}).get("accepted") is True
            for marker in ("not external", "external")
        ),
    }
    report_lines.append("\n## 判据\n")
    report_lines.append("| 判据 | 结果 |")
    report_lines.append("|---|---|")
    for key, ok in criteria.items():
        report_lines.append(f"| {key} | {'通过' if ok else '未通过'} |")
    report_lines.append("| G1.8 | 人工审阅（理论文档），不进自动判定 |")

    report_lines.append("\n## 已知方法局限（决策 0014，不计入 G1.5）\n")
    if limitations:
        for item in limitations:
            terms = "; ".join(
                f"{t['atoms']}: ours={t['ours']:.4g} sobtop={t['sobtop']:.4g}"
                for t in item["evidence"].get("terms", [])
            )
            report_lines.append(
                f"- {item['molecule']}/{item['provider']}/{item['method']}: "
                f"{item['reason']}（Sobtop 同 Hessian 对照：{terms}）"
            )
    else:
        report_lines.append("- 无")

    report_lines.append("\n## 失败项\n")
    if failures:
        for failure in failures:
            report_lines.append(
                f"- {failure.get('molecule')}/{failure.get('provider')}/"
                f"{failure.get('method')}: {failure.get('reason')}"
            )
    else:
        report_lines.append("- 无")

    acceptance = {
        "gate": "G1",
        "passed": all(criteria.values()),
        "criteria": criteria,
        "sobtop_match": sobtop_match,
        "m2_naming": m2_naming,
        "global_fit_self_consistency": gf,
        "cross_engine_all_passed": cross_engine,
        "failures": failures,
        "method_limitations": limitations,
        "pytest_summary": pytest_summary,
    }
    return "\n".join(report_lines) + "\n", "\n".join(csv_lines) + "\n", acceptance


def junit_outcomes(xml_paths: list[Path]) -> dict[str, str]:
    """JUnit testcase → {"tests/unit/x.py::test_name[param]": passed/failed/error/skipped}。"""
    outcomes: dict[str, str] = {}
    for path in xml_paths:
        try:
            cases = ET.parse(path).getroot().iter("testcase")
        except (OSError, ET.ParseError):
            continue
        for case in cases:
            node = case.get("classname", "").replace(".", "/") + ".py::" + case.get("name", "")
            if case.find("error") is not None:
                outcomes[node] = "error"
            elif case.find("failure") is not None:
                outcomes[node] = "failed"
            elif case.find("skipped") is not None:
                outcomes[node] = "skipped"
            else:
                outcomes[node] = "passed"
    return outcomes


def _frequency_metrics(ir, hessian, positions, reference_frequencies) -> dict:
    masses = np.array([a.mass for a in ir.atoms], dtype=np.float64)
    computed = np.sort(harmonic_frequencies(hessian, positions, masses))
    reference = np.sort(np.asarray(reference_frequencies))
    if computed.size != reference.size:
        raise ValueError(f"频率个数不同：FF {computed.size} vs 参考 {reference.size}")
    diff = computed - reference
    below = reference < 1000.0
    return {
        "freq_rmse_cm": float(np.sqrt((diff**2).mean())),
        "freq_mae_cm": float(np.abs(diff).mean()),
        "freq_max_cm": float(np.abs(diff).max()),
        "freq_rmse_lt1000": float(np.sqrt((diff[below] ** 2).mean())) if below.any() else None,
        "freq_rmse_ge1000": float(np.sqrt((diff[~below] ** 2).mean())) if (~below).any() else None,
    }


def _geometry_metrics(ir, reference_positions) -> dict:
    from ffcompiler.bonded.internal import angle_value, bond_value
    from ffcompiler.validation.vibrations import kabsch_rmsd, minimize_ff

    optimized = minimize_ff(ir, reference_positions, tolerance_kj_nm=1e-3)
    bond_errors = [
        abs(bond_value(optimized, b.i, b.j) - bond_value(reference_positions, b.i, b.j)) * 1000.0
        for b in ir.bonds
    ]
    angle_errors = [
        abs(np.degrees(angle_value(optimized, a.i, a.j, a.k_atom)
                       - angle_value(reference_positions, a.i, a.j, a.k_atom)))
        for a in ir.angles
    ]
    return {
        "kabsch_rmsd_pm": kabsch_rmsd(optimized, reference_positions) * 1000.0,
        "bond_mae_pm": float(np.mean(bond_errors)) if bond_errors else None,
        "angle_mae_deg": float(np.mean(angle_errors)) if angle_errors else None,
    }


def _pytest_summary(run_dir: Path) -> tuple[dict, list[Path]]:
    summary: dict = {}
    xmls = []
    for marker, stem in (("not external", "pytest_unit"), ("external", "pytest_external")):
        xml = (run_dir / f"{stem}.xml").resolve()
        xmls.append(xml)
        try:
            proc = subprocess.run(
                ["bash", str(REPO_ROOT / "scripts/run_logged.sh"), sys.executable,
                 "-m", "pytest", "-m", marker, "-q", "--durations=0", f"--junitxml={xml}"],
                capture_output=True, text=True, cwd=REPO_ROOT, timeout=7200,
            )
        except subprocess.TimeoutExpired as exc:
            (run_dir / f"{stem}.log").write_text(f"TIMEOUT after {exc.timeout} s\n")
            summary[marker] = {"accepted": False, "parse_error": f"timeout {exc.timeout} s"}
            continue
        (run_dir / f"{stem}.log").write_text(proc.stdout + proc.stderr)
        summary[marker] = summarize_pytest(xml, proc.returncode)
    return summary, xmls


def self_consistency_error() -> float:
    """G1.3：T22 用例 2 的自洽恢复误差（symmetrize 关/开取最大）。

    ethanol GAFF2 IR 的 bond/angle 平衡值设成 x₀（xtb_ethanol_gfn2 固件几何），
    H_ref = 该 IR 的完整 FF Hessian；先验 = 真值 × U(0.5, 1.5)（种子 20261004），α=0。
    """
    from ffcompiler.bonded.global_fit import global_fit
    from ffcompiler.bonded.internal import angle_value, bond_value
    from ffcompiler.hessian.base import CartesianHessian
    from ffcompiler.hessian.xtb import read_xtb_hessian
    from ffcompiler.ir import Angle, Bond
    from ffcompiler.validation.vibrations import ff_hessian

    fixture = REPO_ROOT / "tests/fixtures/hessian/xtb_ethanol_gfn2"
    x0 = read_xtb_hessian(fixture / "hessian", fixture / "xtbopt.xyz", method="gfn2").positions_nm
    ir, _ = from_amber(*(FIXTURE_DIR / "ethanol" / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    ir0 = ir.replace_terms(
        bonds=tuple(Bond(i=b.i, j=b.j, r0=bond_value(x0, b.i, b.j), k=b.k, source=b.source)
                    for b in ir.bonds),
        angles=tuple(Angle(i=a.i, j=a.j, k_atom=a.k_atom,
                           theta0=angle_value(x0, a.i, a.j, a.k_atom),
                           force_constant=a.force_constant, source=a.source)
                     for a in ir.angles),
    )
    h = CartesianHessian(
        atomic_numbers=np.array([a.atomic_number for a in ir.atoms]), positions_nm=x0,
        hessian=ff_hessian(ir0, x0), gradient=None, energy_kj=None, provider="ff", method="gaff2",
    )
    rng = np.random.default_rng(SEED)
    prior = ir0.replace_terms(
        bonds=tuple(Bond(i=b.i, j=b.j, r0=b.r0, k=b.k * rng.uniform(0.5, 1.5), source=b.source)
                    for b in ir0.bonds),
        angles=tuple(Angle(i=a.i, j=a.j, k_atom=a.k_atom, theta0=a.theta0,
                           force_constant=a.force_constant * rng.uniform(0.5, 1.5),
                           source=a.source)
                     for a in ir0.angles),
    )
    truth = [b.k for b in ir0.bonds] + [a.force_constant for a in ir0.angles]
    worst = 0.0
    for sym in (False, True):
        fitted, _, _ = global_fit(ir0, h, prior_ir=prior, alpha=0.0, symmetrize=sym)
        got = [b.k for b in fitted.bonds] + [a.force_constant for a in fitted.angles]
        worst = max(worst, max(abs(g - t) / t for g, t in zip(got, truth)))
    return float(worst)


def _sobtop_match(failures: list[dict]) -> dict[str, str]:
    """§2.8 逐 case 比对（k、键长、键角；除 sf6_fchk）。"""
    from ffcompiler.bonded.sobtop_ref import read_sobtop_params

    match = {}
    for method, (fname, fn) in ORACLE_METHODS.items():
        n_ok = 0
        for case in GATE_CASES:
            try:
                h, bonds, angles = load_case(case)
                bad = oracle_mismatches(fn(h, bonds, angles),
                                        read_sobtop_params(REFS / case / fname))
            except Exception as exc:  # noqa: BLE001
                bad = [f"异常: {exc}"]
            if bad:
                failures.append({"molecule": case, "provider": "sobtop", "method": method,
                                 "reason": f"{len(bad)} 项超差，例如 {bad[0]}"})
            else:
                n_ok += 1
        match[method] = f"{n_ok}/{len(GATE_CASES)}"
    return match


_SOBTOP_METHOD_NUMBER = {"mseminario": 2, "m2seminario": 3, "drih": 4}


def sobtop_confirm(ir, h, hessian_path: Path, error, workdir: Path) -> dict:
    """决策 0014：用 Sobtop 在同一 Hessian 上复算方法给出的非正 k。

    sobtop_confirmed = Sobtop 对每一个非正项也给出非正 k，且与我们一致（§2.8 容差）。
    Sobtop 失败或任何一项不符 → False（照常算 Gate 失败）。
    """
    import make_sobtop_refs as refs
    from ffcompiler.bonded.sobtop_ref import read_sobtop_params

    workdir.mkdir(parents=True, exist_ok=True)
    try:
        refs.prepare_sobtop()
        mol2 = workdir / "input.mol2"
        refs.write_ir_mol2(ir, h.positions_nm, mol2)
        params = read_sobtop_params(refs.run_sobtop(
            workdir, mol2, _SOBTOP_METHOD_NUMBER[error.method], Path(hessian_path).resolve()
        ))
    except (Exception, SystemExit) as exc:  # noqa: BLE001 - run_sobtop 以 SystemExit 报错
        return {"sobtop_confirmed": False, "reason": f"Sobtop 运行失败: {exc}", "terms": []}
    terms, confirmed = [], True
    for atoms, k in error.terms:
        if len(atoms) == 2:
            ours, ref = k / 100.0, params.bonds[tuple(sorted(atoms))][1]
        else:
            i, j, k_at = atoms
            ours, ref = k, params.angles[(min(i, k_at), j, max(i, k_at))][1]
        agree = abs(ours - ref) <= 1e-4 * abs(ref) + 5e-3
        confirmed = confirmed and agree and ref <= 0.0
        terms.append({"atoms": list(atoms), "ours": ours, "sobtop": ref})
    return {"sobtop_confirmed": confirmed, "sobtop_dir": str(workdir), "terms": terms}


def _evaluate(name, provider, method, ir, positions, h, ref_freqs, run_dir, skip_compare) -> dict:
    from ffcompiler.bonded.assign import derive_bonded
    from ffcompiler.bonded.global_fit import global_fit
    from ffcompiler.validation.vibrations import ff_hessian

    policy = "error"  # Gate 分子都没有线性键角；遇到就报错（决策 0013）
    geom = h.positions_nm  # 参考几何 = Hessian 几何（全部方法同一参考）
    row = {"molecule": name, "provider": provider, "method": method}
    if method == "gaff2":
        new_ir = ir
        row["max_class_spread"] = 0.0
    elif method == "global_fit":
        new_ir, geom, gf = global_fit(
            ir, h, linear_angle_policy=policy, input_positions_nm=positions,
        )
        row.update({
            "gf_alpha": gf.alpha, "gf_lam": gf.lam, "gf_cond": gf.cond,
            "gf_rel_res_fit": gf.rel_residual_fit, "gf_rel_res_prior": gf.rel_residual_prior,
            "gf_n_at_bound": gf.n_at_bound, "gf_max_force": gf.max_force_at_x0,
        })
    else:
        new_ir, geom, report = derive_bonded(
            ir, h, method=method, linear_angle_policy=policy, input_positions_nm=positions,
        )
        row["max_class_spread"] = max(report.class_spread.values(), default=0.0)
    row.update(_frequency_metrics(new_ir, ff_hessian(new_ir, geom), geom, ref_freqs))
    row.update(_geometry_metrics(new_ir, geom))
    if skip_compare:
        row["compare_passed"] = None
    else:
        from ffcompiler.validation.engines import compare

        mdir = FIXTURE_DIR / name
        # sander-tleap 腿比的是原 GAFF2 Hamiltonian：只对 gaff2 基线有意义
        tleap = {"tleap_prmtop": mdir / f"{name}.prmtop",
                 "tleap_inpcrd": mdir / f"{name}.inpcrd"} if method == "gaff2" else {}
        result = compare(new_ir, geom, run_dir / "compare" / f"{name}_{provider}_{method}",
                         **tleap)
        row["compare_passed"] = bool(result.passed)
    return row


def evaluate_matrix(keys, systems, hessian_for, run_dir: Path, skip_compare: bool,
                    *, evaluate=None) -> tuple[list[dict], list[dict]]:
    """逐 (分子, provider) × 方法评估；每个异常都记进 failures，不丢弃（G1.5）。

    hessian_for(key) 返回 CartesianHessian 或抛错（xtb 失败记到该组合的每个方法）。
    方法给出非正 k 时现场用 Sobtop 复算（决策 0014）。
    """
    evaluate = evaluate or _evaluate
    rows: list[dict] = []
    failures: list[dict] = []
    for name, provider in keys:
        ir, positions = systems[name]
        try:
            h = hessian_for((name, provider))
            diag = check_hessian(h, np.array([a.mass for a in ir.atoms]))
        except Exception as exc:  # noqa: BLE001
            for method in METHODS:
                failures.append({"molecule": name, "provider": provider,
                                 "method": method, "reason": f"xtb/check: {exc}"})
            continue
        for method in METHODS:
            try:
                rows.append(evaluate(name, provider, method, ir, positions, h,
                                     diag.frequencies_cm, run_dir, skip_compare))
            except NonPositiveForceConstantError as exc:
                evidence = sobtop_confirm(
                    ir, h, run_dir / "xtb" / provider / name / "hessian", exc,
                    run_dir / "sobtop_confirm" / f"{name}_{provider}_{method}",
                )
                failures.append({"molecule": name, "provider": provider,
                                 "method": method, "reason": str(exc),
                                 "category": "method_limitation", "evidence": evidence})
            except Exception as exc:  # noqa: BLE001
                failures.append({"molecule": name, "provider": provider,
                                 "method": method, "reason": f"{type(exc).__name__}: {exc}"})
    return rows, failures


def main() -> int:
    parser = argparse.ArgumentParser(prog="gate_mvp1")
    parser.add_argument("--providers", nargs="+", default=["gfnff", "gfn2"],
                        choices=["gfnff", "gfn2"])
    parser.add_argument("--jobs", type=int, default=1, help="并行 xtb 运行数")
    parser.add_argument("--skip-compare", action="store_true",
                        help="调试用：跳过跨引擎比较（Gate 将判为未通过）")
    args = parser.parse_args()

    with RunContext("gate_mvp1", root=REPO_ROOT / "jobs", seed=SEED, inputs=[
        FIXTURE_DIR / n / f"{n}.{ext}" for n in NAMES
        for ext in ("prmtop", "inpcrd", "sdf")
    ]) as ctx:
        run_dir = ctx.dir
        systems = {
            n: from_amber(*(FIXTURE_DIR / n / f"{n}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
            for n in NAMES
        }

        def run_one(key):
            name, provider = key
            ir, positions = systems[name]
            return run_xtb_hessian(
                np.array([a.atomic_number for a in ir.atoms]), positions,
                method=provider, charge=ir.total_charge,
                workdir=run_dir / "xtb" / provider / name,
            )

        keys = [(n, p) for n in NAMES for p in args.providers]
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {key: pool.submit(run_one, key) for key in keys}

        rows, failures = evaluate_matrix(
            keys, systems, lambda key: futures[key].result(), run_dir, args.skip_compare,
        )

        sobtop_match = _sobtop_match(failures)
        try:
            self_consistency = self_consistency_error()
        except Exception as exc:  # noqa: BLE001
            self_consistency = None
            failures.append({"molecule": "ethanol", "provider": "ff", "method": "global_fit",
                             "reason": f"自洽恢复计算失败: {exc}"})
        pytest_summary, xmls = _pytest_summary(run_dir)
        outcomes = junit_outcomes(xmls)
        review_focus = {
            label: [(test_id, outcomes.get(test_id, "missing")) for test_id in tests]
            for label, tests in REVIEW_FOCUS_TESTS.items()
        }
        report_md, metrics_csv, acceptance = build_report(
            rows, failures,
            pytest_summary=pytest_summary, sobtop_match=sobtop_match,
            m2_naming="m2seminario", global_fit_self_consistency=self_consistency,
            review_focus=review_focus,
            expected=[(n, p, m) for n, p in keys for m in METHODS],
        )
        (run_dir / "report.md").write_text(report_md)
        (run_dir / "metrics.csv").write_text(metrics_csv)
        (run_dir / "acceptance.json").write_text(
            json.dumps(acceptance, indent=2, ensure_ascii=False) + "\n"
        )
        ctx.manifest["gate"] = "G1"
        print(f"G1 passed={acceptance['passed']} failures={len(acceptance['failures'])}")
        return 0 if acceptance["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
