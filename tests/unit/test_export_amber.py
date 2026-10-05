"""T7：AMBER 导出器的必测用例。"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import openmm
import openmm.app
import pytest
from openmm import unit

from ffcompiler.exporters.amber import write_amber
from ffcompiler.formats.amber import read_prmtop
from ffcompiler.ir import ForceFieldIR
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


def _openmm_by_category(prmtop: Path, positions: np.ndarray) -> dict[str, float]:
    amber = openmm.app.AmberPrmtopFile(str(prmtop))
    system = amber.createSystem(
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
    out = write_amber(ir, positions, tmp_path / name)
    categories = _openmm_by_category(out["prmtop"], positions)

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


def test_no_torsion_14_pair_still_scaled(tmp_path: Path) -> None:
    """Review Focus 2：没有 torsion 项的 1-4 对不能丢（载体条目兜底）。"""
    from ffcompiler.ir import Angle, Atom, Bond

    def atom(i, q, elem="C", z=6):
        return Atom(index=i, name=f"A{i}", element=elem, atomic_number=z, mass=12.0,
                    atom_type=f"T{i}", formal_charge=0, partial_charge=q,
                    sigma=0.3, epsilon=0.4)

    atoms = (atom(0, -0.1), atom(1, 0.1), atom(2, 0.1), atom(3, -0.1))
    ir = ForceFieldIR(
        name="no_torsion_14", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.15, k=3.0e5, source="t"),
            Bond(i=1, j=2, r0=0.15, k=3.0e5, source="t"),
            Bond(i=2, j=3, r0=0.15, k=3.0e5, source="t"),
        ),
        angles=(
            Angle(i=0, j=1, k_atom=2, theta0=1.9, force_constant=500.0, source="t"),
            Angle(i=1, j=2, k_atom=3, theta0=1.9, force_constant=500.0, source="t"),
        ),
        proper_torsions=(),  # 全部删掉
        impropers=(), exclusions=(), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    exclusions, pairs14 = expected_exclusions_and_pairs(ir)
    assert pairs14 == ((0, 3),)
    ir = ForceFieldIR(
        name=ir.name, atoms=ir.atoms, bonds=ir.bonds, angles=ir.angles,
        proper_torsions=ir.proper_torsions, impropers=ir.impropers,
        exclusions=exclusions, pairs14=pairs14, total_charge=ir.total_charge,
        forcefield_family=ir.forcefield_family, combination_rule=ir.combination_rule,
        coulomb_14_scale=ir.coulomb_14_scale, vdw_14_scale=ir.vdw_14_scale,
    )

    rng = np.random.default_rng(20261004)
    positions = rng.normal(0.0, 0.15, size=(4, 3))
    out = write_amber(ir, positions, tmp_path)
    categories = _openmm_by_category(out["prmtop"], positions)
    ref = evaluate(ir, positions)
    nonb = ref["lj"] + ref["coulomb"] + ref["lj14"] + ref["coulomb14"]
    assert ref["lj14"] != 0.0  # 1-4 确实有贡献
    assert abs(categories["nonbonded_all"] - nonb) <= TOL_FILE_E(nonb)
    assert categories["torsion_all"] == 0.0


@pytest.mark.parametrize("name", ["cyclobutane", "norbornane"])
def test_multipath_14_counted_once(name: str, tmp_path: Path) -> None:
    ir, positions = _fixture_ir(name)
    out = write_amber(ir, positions, tmp_path / name)
    categories = _openmm_by_category(out["prmtop"], positions)
    ref = evaluate(ir, positions)
    nonb = ref["lj"] + ref["coulomb"] + ref["lj14"] + ref["coulomb14"]
    assert abs(categories["nonbonded_all"] - nonb) <= TOL_FILE_E(nonb)


@pytest.mark.parametrize("name", ["acetate", "methylammonium"])
def test_charges_roundtrip(name: str, tmp_path: Path) -> None:
    ir, positions = _fixture_ir(name)
    out = write_amber(ir, positions, tmp_path / name)
    charges = [q / 18.2223 for q in read_prmtop(out["prmtop"])["CHARGE"]]
    assert abs(sum(charges) - ir.total_charge) <= 1e-6


def test_multi_term_propers_counted(tmp_path: Path) -> None:
    ir, positions = _fixture_ir("nma")
    out = write_amber(ir, positions, tmp_path / "nma")
    prm = read_prmtop(out["prmtop"])
    entries = prm["DIHEDRALS_WITHOUT_HYDROGEN"] + prm["DIHEDRALS_INC_HYDROGEN"]
    n_proper = sum(1 for k in range(0, len(entries), 5) if entries[k + 3] > 0)
    assert n_proper == len(ir.proper_torsions)


def _assert_same_ir(a, b, path="ir"):
    import dataclasses
    import math

    if dataclasses.is_dataclass(a):
        assert type(a) is type(b), path
        for f in dataclasses.fields(a):
            if f.name in ("provenance", "name"):
                continue
            _assert_same_ir(getattr(a, f.name), getattr(b, f.name), f"{path}.{f.name}")
    elif isinstance(a, (tuple, list)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_same_ir(x, y, f"{path}[{i}]")
    elif isinstance(a, dict):
        assert set(a) == set(b), path
        for key in a:
            _assert_same_ir(a[key], b[key], f"{path}[{key!r}]")
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-8, abs_tol=1e-12), f"{path}: {a!r} vs {b!r}"
    else:
        assert a == b, f"{path}: {a!r} vs {b!r}"


@pytest.mark.parametrize("name", NAMES)
def test_ir_roundtrip_through_native_prmtop(name: str, tmp_path: Path) -> None:
    """决策 0015：IR → 原生 prmtop/rst7 → 原生导入，逐项还原（E16.8 精度内），含顺序、
    方向、1-4 载体与排除表。"""
    ir, positions = _fixture_ir(name)
    out = write_amber(ir, positions, tmp_path / name)
    back, pos_back = from_amber(out["prmtop"], out["inpcrd"], FIXTURES / name / f"{name}.sdf")
    _assert_same_ir(ir, back)
    np.testing.assert_allclose(pos_back, positions, rtol=0, atol=1e-8)


def test_atom_zero_never_carries_a_sign(tmp_path: Path) -> None:
    """原子 0 落在第 3/4 位时整体反序写出（负号无法作用于 0）；能量不变。

    把全部二面角整体反序（φ 不变），强制制造原子 0 在第 4 位的条目。
    """
    ir, positions = _fixture_ir("ethanol")
    flipped = ir.replace_terms(
        proper_torsions=tuple(replace(t, i=t.l, j=t.k, k=t.j, l=t.i) for t in ir.proper_torsions)
    )
    assert any(t.l == 0 for t in flipped.proper_torsions)
    out = write_amber(flipped, positions, tmp_path / "e")
    prm = read_prmtop(out["prmtop"])
    entries = prm["DIHEDRALS_WITHOUT_HYDROGEN"] + prm["DIHEDRALS_INC_HYDROGEN"]
    for k in range(0, len(entries), 5):
        assert entries[k + 2] != 0 and entries[k + 3] != 0
    categories = _openmm_by_category(out["prmtop"], positions)
    ref = evaluate(ir, positions)
    assert abs(categories["total"] - ref["total"]) <= TOL_FILE_E(ref["total"])
