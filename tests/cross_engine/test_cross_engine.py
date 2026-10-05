"""T9：引擎适配器与跨引擎比较的必测用例。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import openmm
import openmm.app
import pytest
from openmm import unit

from ffcompiler.importers import from_amber
from ffcompiler.ir.energy import evaluate
from ffcompiler.validation.engines import (
    CrossEngineReport,
    openmm_energy_forces,
    parse_mdout_energies,
    sander_energies,
    compare,
)
from ffcompiler.validation.engines import _SANDER_E, _SANDER_EEL

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _fixture_ir(name: str):
    mdir = FIXTURES / name
    return from_amber(mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf")


def _perturbations(ir, positions, rng):
    return [positions] + [positions + rng.normal(0.0, 0.005, size=positions.shape)
                          for _ in range(3)]


@pytest.mark.external
def test_sander_self_check_ethanol(tmp_path: Path) -> None:
    """sander 的 mdin 与解析自检：对 tleap 的 ethanol prmtop，sander 与
    OpenMM `AmberPrmtopFile` 的能量满足 sander 容差。"""
    mdir = FIXTURES / "ethanol"
    ir, positions = _fixture_ir("ethanol")

    sander = sander_energies(mdir / "ethanol.prmtop", mdir / "ethanol.inpcrd", tmp_path)

    amber = openmm.app.AmberPrmtopFile(str(mdir / "ethanol.prmtop"))
    system = amber.createSystem(
        nonbondedMethod=openmm.app.NoCutoff, constraints=None, rigidWater=False,
        removeCMMotion=False,
    )
    energies, _ = openmm_energy_forces(system, positions)
    ref = evaluate(ir, positions)

    checks = {
        "bond": (sander["bond"], ref["bond"] / 4.184),
        "angle": (sander["angle"], ref["angle"] / 4.184),
        "dihed": (sander["dihed"], (ref["proper"] + ref["improper"]) / 4.184),
        "nonbonded_all": (
            sander["vdw"] + sander["eel"] + sander["vdw14"] + sander["eel14"],
            energies["nonbonded_all"] / 4.184,
        ),
    }
    for key, (s_val, ref_val) in checks.items():
        tol = _SANDER_EEL(ref_val) if "eel" in key else _SANDER_E(ref_val)
        assert abs(s_val - ref_val) <= tol, f"{key}: sander {s_val} vs ref {ref_val}"


def test_parse_mdout_truncated_raises() -> None:
    # 与 Amber 24 真实输出格式一致："BOND    = ...  ANGLE   = ...  DIHED ..."
    complete = (
        "   NSTEP           ENERGY          RMS            GMAX         NAME    NUMBER\n"
        "       1      -2.0192E+00     4.3248E+00     1.6780E+01     C2          2\n"
        "\n"
        " BOND    =        2.5000  ANGLE   =        3.5000  DIHED      =        1.0000\n"
        " VDWAALS =        0.5000  EEL     =       -4.0000  HBOND      =        0.0000\n"
        " 1-4 VDW =        0.1000  1-4 EEL =        0.2000  RESTRAINT  =        0.0000\n"
    )
    energies = parse_mdout_energies(complete)
    assert energies["bond"] == 2.5 and energies["eel"] == -4.0 and energies["eel14"] == 0.2

    truncated = complete[: complete.index("VDWAALS")]
    with pytest.raises(ValueError, match="missing energies"):
        parse_mdout_energies(truncated)
    with pytest.raises(ValueError, match="no energy block"):
        parse_mdout_energies("nothing here")


def test_report_markdown_and_json() -> None:
    report = CrossEngineReport()
    report.rows.append({
        "leg": "openmm-ir", "reference": "numpy-ir", "category": "bond",
        "E_leg": 1.0, "E_ref": 1.0, "delta": 0.0, "tol": 1e-6, "passed": True,
    })
    report.force_rows.append({
        "leg": "openmm-xml", "reference": "openmm-ir", "max_df": 1e-9,
        "tol": 1e-6, "passed": True,
    })
    report.skipped["gmx"] = "gmx: not installed (roadmap decision D1)"

    md = report.to_markdown()
    assert md.splitlines()[0].startswith("| leg | reference | category |")
    assert md.splitlines()[1] == "|---|---|---|---|---|---|---|---|"
    assert "SKIPPED gmx" in md
    payload = report.to_json()
    json.dumps(payload)  # 可序列化
    assert payload["passed"] is True


@pytest.mark.external
@pytest.mark.parametrize("name", [
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
])
def test_compare_all_fixtures_pass(name: str, tmp_path: Path) -> None:
    ir, positions = _fixture_ir(name)
    rng = np.random.default_rng(20261004)
    for n, pos in enumerate(_perturbations(ir, positions, rng)):
        report = compare(
            ir, pos, tmp_path / f"cfg{n}",
            tleap_prmtop=FIXTURES / name / f"{name}.prmtop",
            tleap_inpcrd=FIXTURES / name / f"{name}.inpcrd",
        )
        assert report.passed, (
            f"{name} config {n} failed:\n"
            + "\n".join(str(r) for r in report.rows if not r["passed"])
            + "\n".join(str(r) for r in report.force_rows if not r["passed"])
        )


@pytest.mark.external
def test_compare_detects_tampered_prmtop(tmp_path: Path) -> None:
    """T9 用例 4：故意坏 IR（bond k×1.001 只进 prmtop）必须被定位到 bond。"""
    ir, positions = _fixture_ir("ethanol")
    bad_bond = replace(ir.bonds[0], k=ir.bonds[0].k * 1.001)
    bad_ir = ir.replace_terms(
        bonds=(bad_bond,) + tuple(b for b in ir.bonds if b is not ir.bonds[0])
    )
    assert bad_ir.bonds[0].k == pytest.approx(ir.bonds[0].k * 1.001)

    report = compare(ir, positions, tmp_path, amber_ir=bad_ir)
    assert not report.passed
    failing = {(r["leg"], r["category"]) for r in report.rows if not r["passed"]}
    # 决策 0004.3：sander 的绝对容差恰好容纳 1.14e-4 kcal 的篡改；以
    # openmm-prmtop 行为准（计划原文：失败行指向 bond 类别）。
    assert any(leg == "openmm-prmtop" and cat == "bond" for leg, cat in failing), failing


@pytest.mark.external
def test_compare_detects_tampered_gromacs_topology(tmp_path: Path) -> None:
    """决策 0013：放宽到 mixed 容差后，gmx 腿仍能定位真实拓扑错误（bond k×1.001 只进 itp）。"""
    from ffcompiler.validation.gromacs import find_gromacs

    if find_gromacs() is None:
        pytest.skip("gmx not installed")
    ir, positions = _fixture_ir("ethanol")
    bad_bond = replace(ir.bonds[0], k=ir.bonds[0].k * 1.001)
    bad_ir = ir.replace_terms(
        bonds=(bad_bond,) + tuple(b for b in ir.bonds if b is not ir.bonds[0])
    )
    rng = np.random.default_rng(20261004)
    pos = _perturbations(ir, positions, rng)[1]
    report = compare(ir, pos, tmp_path, gromacs_ir=bad_ir)
    failing = {(r["leg"], r["category"]) for r in report.rows if not r["passed"]}
    assert ("gmx", "bond") in failing, failing
    clean = compare(ir, pos, tmp_path / "clean")
    assert not any(r["leg"] == "gmx" and not r["passed"] for r in clean.rows)


@pytest.mark.external
def test_acetonitrile_diagnostic(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """乙腈不进 Gate（决策 0013，近线性角 acos 病态）：只要求各引擎跑通，超差行打印不断言。"""
    ir, positions = _fixture_ir("acetonitrile")
    rng = np.random.default_rng(20261004)
    for n, pos in enumerate(_perturbations(ir, positions, rng)):
        report = compare(ir, pos, tmp_path / f"cfg{n}")
        assert not report.errors, report.errors
        for row in report.rows:
            if not row["passed"]:
                print(f"acetonitrile cfg{n}: {row}")
