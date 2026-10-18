"""T26：ASE 约束松弛扫描与 MM 工具（MVP-2 计划 §2.2、§2.3）。参考势用 GAFF2 的 MM，不用 MACE。"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes
from rdkit import Chem

from ffcompiler.exporters.openmm import to_openmm_system
from ffcompiler.importers import from_amber
from ffcompiler.ir.energy import forces_fd
from ffcompiler.torsion.mm import OpenMMCalculator, mm_energy_without_bond
from ffcompiler.torsion.rotors import find_rotors
from ffcompiler.torsion.scan import ScanError, relaxed_scan
from ffcompiler.units import EV_KJ_MOL
from ffcompiler.validation.engines import reference_context

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _ethanol():
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    rotor = {r.bond: r for r in find_rotors(ir, mol, pos)}[(1, 2)]  # C–O
    return ir, pos, rotor


def _openmm_energy(ir, pos_nm) -> float:
    import openmm
    from openmm import unit

    ctx = reference_context(to_openmm_system(ir))
    ctx.setPositions(pos_nm)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


@pytest.fixture(scope="module")
def ethanol_scan():
    ir, pos, rotor = _ethanol()
    numbers = [a.atomic_number for a in ir.atoms]
    bonds = [(b.i, b.j) for b in ir.bonds]
    result = relaxed_scan(OpenMMCalculator(ir), numbers, pos, rotor.scan_quartet, bonds)
    return ir, result


def test_constraint_holds_on_all_grid_points(ethanol_scan) -> None:
    _, res = ethanol_scan
    assert res.phi_rad.shape == (36,) and res.energy.shape == (36,)
    dev = np.degrees(np.angle(np.exp(1j * (res.phi_rad - res.phi_target_rad))))
    assert np.max(np.abs(dev)) <= 1e-6
    steps = np.degrees(np.angle(np.exp(1j * np.diff(res.phi_target_rad))))
    assert np.allclose(np.abs(steps), 10.0, atol=1e-6)
    assert res.converged.all()


def test_energy_is_reference_potential_itself(ethanol_scan) -> None:
    ir, res = ethanol_scan
    for s in (0, 7, 19):
        assert res.energy[s] == pytest.approx(_openmm_energy(ir, res.positions[s]), abs=1e-8)


class _NaNCalc(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {"energy": float("nan"), "forces": np.full((len(atoms), 3), np.nan)}


def test_nan_reference_raises_scan_error() -> None:
    ir, pos, rotor = _ethanol()
    with pytest.raises(ScanError):
        relaxed_scan(_NaNCalc(), [a.atomic_number for a in ir.atoms], pos,
                     rotor.scan_quartet, [(b.i, b.j) for b in ir.bonds])


def test_openmm_calculator_forces_match_numpy_fd() -> None:
    from ase import Atoms

    ir, pos, _ = _ethanol()
    atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos * 10.0)
    atoms.calc = OpenMMCalculator(ir)
    f_kj_nm = atoms.get_forces() * EV_KJ_MOL * 10.0
    assert np.max(np.abs(f_kj_nm - forces_fd(ir, pos))) <= 1e-4
    assert atoms.get_potential_energy() * EV_KJ_MOL == pytest.approx(_openmm_energy(ir, pos), abs=1e-8)


class _InfoRecorder(OpenMMCalculator):
    seen: list = []

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        type(self).seen.append(dict(atoms.info))
        super().calculate(atoms, properties, system_changes)


def test_info_is_passed_to_calculator() -> None:
    ir, pos, rotor = _ethanol()
    _InfoRecorder.seen = []
    relaxed_scan(_InfoRecorder(ir), [a.atomic_number for a in ir.atoms], pos, rotor.scan_quartet,
                 [(b.i, b.j) for b in ir.bonds], info={"charge": -1, "spin": 1}, step_deg=120.0)
    assert _InfoRecorder.seen and all(i.get("charge") == -1 and i.get("spin") == 1 for i in _InfoRecorder.seen)


def test_mm_energy_without_bond_drops_only_that_bond() -> None:
    ir, pos, rotor = _ethanol()
    j, k = rotor.bond
    full = _openmm_energy(ir, pos)
    on_bond = [t for t in ir.proper_torsions if {t.j, t.k} == {j, k}]
    assert on_bond
    from ffcompiler.ir.energy import dihedral_angle

    e_bond = sum(t.force_constant * (1 + math.cos(t.periodicity * dihedral_angle(*(pos[a] for a in (t.i, t.j, t.k, t.l))) - t.phase))
                 for t in on_bond)
    assert mm_energy_without_bond(ir, rotor.bond, pos) == pytest.approx(full - e_bond, abs=1e-8)
    stack = np.stack([pos, pos])
    assert mm_energy_without_bond(ir, rotor.bond, stack).shape == (2,)
