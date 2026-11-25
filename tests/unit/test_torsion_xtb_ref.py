"""torsion/xtb_ref.py：xtb GFN2 作二面角参考势（决策 0028）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.torsion.xtb_ref import MODEL_NAME, SUPPORTED, XtbGFN2Calculator, is_ion_or_zwitterion

GAFF2 = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _load(name):
    from ffcompiler.importers import from_amber

    d = GAFF2 / name
    return from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")


def test_ion_or_zwitterion_judgement():
    assert not is_ion_or_zwitterion(_load("ethanol")[0])
    assert is_ion_or_zwitterion(_load("acetate")[0])
    assert is_ion_or_zwitterion(_load("methylammonium")[0])


def test_identity_constants():
    assert MODEL_NAME == "xtb-gfn2"
    assert min(SUPPORTED) == 1 and max(SUPPORTED) == 86


@pytest.mark.external
def test_forces_are_negative_energy_gradient(tmp_path):
    """单位与符号：解析力 = −dE/dx（中心差分 1e-3 Å），甲铵正离子（charge +1）。"""
    from ase import Atoms

    ir, pos = _load("methylammonium")
    atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos * 10.0)
    atoms.calc = XtbGFN2Calculator(int(ir.total_charge), threads=2)
    f = atoms.get_forces()
    h = 1e-3
    for a, k in [(0, 0), (1, 2), (4, 1)]:
        p = atoms.positions.copy()
        p[a, k] += h
        atoms.positions = p
        e_plus = atoms.get_potential_energy()
        p[a, k] -= 2 * h
        atoms.positions = p
        e_minus = atoms.get_potential_energy()
        p[a, k] += h
        atoms.positions = p
        fd = -(e_plus - e_minus) / (2 * h)
        assert f[a, k] == pytest.approx(fd, abs=2e-4), (a, k)


@pytest.mark.external
def test_energy_matches_direct_singlepoint(tmp_path):
    from ase import Atoms

    from ffcompiler.hessian.xtb import run_xtb_singlepoint

    ir, pos = _load("acetate")
    z = np.array([a.atomic_number for a in ir.atoms])
    direct = run_xtb_singlepoint(z, pos, method="gfn2", charge=-1, workdir=tmp_path, threads=2)
    atoms = Atoms(numbers=z, positions=pos * 10.0)
    atoms.calc = XtbGFN2Calculator(-1, threads=2)
    assert atoms.get_potential_energy() == pytest.approx(direct["energy_hartree"] * 27.211386245988,
                                                         rel=0, abs=1e-9)
