"""T8：GROMACS 导出器的必测用例。"""

from __future__ import annotations

import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import openmm
import openmm.app
import pytest
from openmm import unit

from ffcompiler.exporters.gromacs import write_gromacs
from ffcompiler.ir import ForceFieldIR, ImproperTorsion
from ffcompiler.ir.energy import evaluate
from ffcompiler.ir.topology import expected_exclusions_and_pairs
from ffcompiler.importers import from_amber
from ffcompiler.validation.engines import reference_context

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
    "acetonitrile",
)


def _fixture_ir(name: str):
    mdir = FIXTURES / name
    return from_amber(mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf")


def _read_back_openmm(top: Path, positions: np.ndarray) -> dict[str, float]:
    gmx_top = openmm.app.GromacsTopFile(str(top))
    system = gmx_top.createSystem(
        nonbondedMethod=openmm.app.NoCutoff, constraints=None, rigidWater=False,
        removeCMMotion=False,
    )
    # 力组必须在建 Context 之前设好（Context 快照 System，之后的修改无效）
    for idx in range(system.getNumForces()):
        system.getForce(idx).setForceGroup(idx)
    context = reference_context(system)
    context.setPositions([openmm.Vec3(*p) for p in positions] * unit.nanometer)
    categories = {"bond": 0.0, "angle": 0.0, "torsion_all": 0.0, "nonbonded_all": 0.0}
    for idx in range(system.getNumForces()):
        cls = system.getForce(idx).__class__.__name__
        key = {
            "HarmonicBondForce": "bond",
            "HarmonicAngleForce": "angle",
            "PeriodicTorsionForce": "torsion_all",
            "NonbondedForce": "nonbonded_all",
            "CustomTorsionForce": "torsion_all",
        }.get(cls)
        if key is None:
            continue
        state = context.getState(getEnergy=True, groups=1 << idx)
        categories[key] += state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    del context
    categories["total"] = total
    return categories


TOL_FILE_E = lambda e: 1e-4 + 1e-7 * abs(e)


@pytest.mark.parametrize("name", NAMES)
def test_roundtrip_energy_matches(name: str, tmp_path: Path) -> None:
    ir, positions = _fixture_ir(name)
    out = write_gromacs(ir, positions, tmp_path / name)
    assert set(out) == {"itp", "top", "gro"}
    categories = _read_back_openmm(out["top"], positions)

    ref = evaluate(ir, positions)
    checks = {
        "bond": ref["bond"],
        "angle": ref["angle"],
        "torsion_all": ref["proper"] + ref["improper"],
        "nonbonded_all": ref["lj"] + ref["coulomb"] + ref["lj14"] + ref["coulomb14"],
        "total": ref["total"],
    }
    for key, expected in checks.items():
        assert abs(categories[key] - expected) <= TOL_FILE_E(expected), (
            f"{name}: {key} readback {categories[key]} vs IR {expected}"
        )


def test_multi_term_proper_three_lines(tmp_path: Path) -> None:
    """Review Focus 4：同一四元组 n=1,2,3 → itp 里 3 行 funct 9，能量一致。"""
    from ffcompiler.ir import ProperTorsion
    from tests.unit.test_energy import _chain_ir

    ir = _chain_ir().replace_terms(
        proper_torsions=tuple(
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=n, phase=0.1 * n,
                          force_constant=1.0 - 0.1 * n, source="t")
            for n in (1, 2, 3)
        )
    )
    rng = np.random.default_rng(20261004)
    positions = rng.normal(0.0, 0.15, size=(4, 3))
    out = write_gromacs(ir, positions, tmp_path)

    dihedral_lines = [
        line.split()
        for line in out["itp"].read_text().splitlines()
        if line.strip() and not line.startswith((";", "["))
    ]
    proper_lines = [f for f in dihedral_lines if len(f) >= 8 and f[4] == "9"]
    assert len(proper_lines) == 3
    assert [int(f[7]) for f in proper_lines] == [1, 2, 3]

    categories = _read_back_openmm(out["top"], positions)
    ref = evaluate(ir, positions)
    assert abs(categories["torsion_all"] - ref["proper"]) <= TOL_FILE_E(ref["proper"])


def test_exported_coordinates_keep_reference_precision(tmp_path):
    ir, pos = _fixture_ir("ethanol")
    pos = pos + np.random.default_rng(20261004).normal(0, 0.005, pos.shape)
    out = write_gromacs(ir, pos, tmp_path)
    gro = openmm.app.GromacsGroFile(str(out["gro"]))
    readback = np.array(gro.positions.value_in_unit(unit.nanometer))
    assert np.max(np.abs(readback - pos)) <= 5.1e-11


def test_two_impropers_same_quartet_raises(tmp_path: Path) -> None:
    ir, positions = _fixture_ir("nma")
    # 复制第一个 improper 换一个 periodicity → 同四元组两项（校验允许，GROMACS 不允许）
    second = replace(
        ir.impropers[0],
        periodicity=3 if ir.impropers[0].periodicity != 3 else 5,
    )
    ir2 = ir.replace_terms(impropers=ir.impropers + (second,))
    with pytest.raises(ValueError, match=r"quartet .* has 2 terms"):
        write_gromacs(ir2, positions, tmp_path)


def test_charged_atom_block_charge_sum(tmp_path: Path) -> None:
    for name, charge in (("acetate", -1), ("methylammonium", 1)):
        ir, positions = _fixture_ir(name)
        out = write_gromacs(ir, positions, tmp_path / name)
        lines = out["itp"].read_text().splitlines()
        start = lines.index("[ atoms ]")
        total = 0.0
        for line in lines[start + 1:]:
            if line.startswith("["):
                break
            if not line.strip() or line.startswith(";"):
                continue
            total += float(line.split()[6])
        assert abs(total - charge) <= 1e-8, name


def test_pairs_row_count_and_norbornane(tmp_path: Path) -> None:
    for name in ("ethanol", "norbornane"):
        ir, positions = _fixture_ir(name)
        out = write_gromacs(ir, positions, tmp_path / name)
        lines = out["itp"].read_text().splitlines()
        start = lines.index("[ pairs ]")
        n_rows = 0
        for line in lines[start + 1:]:
            if line.startswith("[") or not line.strip() or line.startswith(";"):
                break
            n_rows += 1
        assert n_rows == len(ir.pairs14), name
    # norbornane 能量一致（Review Focus 1，多路径 1-4）
    ir, positions = _fixture_ir("norbornane")
    out = write_gromacs(ir, positions, tmp_path / "nb")
    categories = _read_back_openmm(out["top"], positions)
    ref = evaluate(ir, positions)
    nonb = ref["lj"] + ref["coulomb"] + ref["lj14"] + ref["coulomb14"]
    assert abs(categories["nonbonded_all"] - nonb) <= TOL_FILE_E(nonb)


@pytest.mark.external
def test_real_gromacs_rerun(tmp_path: Path) -> None:
    """真 gmx 腿：没有 gmx 时 skip（决策 D1）。"""
    from ffcompiler.validation.gromacs import find_gromacs, energies
    if find_gromacs() is None:
        pytest.skip("gmx not installed")
    # 有 gmx 的机器上：grompp + mdrun -rerun 单点能量对照
    ir, positions = _fixture_ir("ethanol")
    out = write_gromacs(ir, positions, tmp_path)
    observed, metadata = energies(ir, positions, out, tmp_path / "real_gmx")
    ref = evaluate(ir, positions)
    assert metadata["precision"] in ("mixed", "double", "single")
    for category in ref:
        assert abs(observed[category] - ref[category]) <= TOL_FILE_E(ref[category])
