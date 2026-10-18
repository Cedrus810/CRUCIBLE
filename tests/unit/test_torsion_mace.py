"""T27：MACE-POLAR 计算器（MVP-2 计划 §0、§1.1）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.importers import from_amber
from ffcompiler.torsion.mace import default_model_path, molecule_info

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _load(name: str):
    mdir = FIXTURES / name
    return from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))


def test_default_model_path_env_override(monkeypatch, tmp_path) -> None:
    fake = tmp_path / "m.model"
    fake.write_bytes(b"x")
    monkeypatch.setenv("FFC_MACE_MODEL", str(fake))
    assert default_model_path() == fake
    monkeypatch.setenv("FFC_MACE_MODEL", str(tmp_path / "missing.model"))
    with pytest.raises(FileNotFoundError):
        default_model_path()


def test_default_model_path_is_polar_m(monkeypatch) -> None:
    monkeypatch.delenv("FFC_MACE_MODEL", raising=False)
    assert default_model_path() == Path("/home/user/MLP/mace/MACE-POLAR-1-M.model")


def test_molecule_info_carries_total_charge() -> None:
    ir, _ = _load("acetate")
    assert molecule_info(ir) == {"charge": -1, "spin": 1}


def _energy_and_charges(name: str):
    from ase import Atoms

    from ffcompiler.torsion.mace import polar_calculator

    ir, pos = _load(name)
    atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos * 10.0)
    atoms.info.update(molecule_info(ir))
    atoms.calc = polar_calculator()
    from ffcompiler.units import EV_KJ_MOL

    return atoms.get_potential_energy() * EV_KJ_MOL, atoms.calc.results["charges"]


@pytest.mark.external
@pytest.mark.mlp
def test_ethanol_energy_matches_smoke() -> None:
    e, _ = _energy_and_charges("ethanol")
    assert e == pytest.approx(-407018.399, abs=1e-2)


@pytest.mark.external
@pytest.mark.mlp
def test_acetate_charges_sum_to_minus_one() -> None:
    _, q = _energy_and_charges("acetate")
    assert float(np.sum(q)) == pytest.approx(-1.0, abs=1e-6)


@pytest.mark.external
@pytest.mark.mlp
def test_supported_atomic_numbers() -> None:
    from ffcompiler.torsion.mace import supported_atomic_numbers

    z = supported_atomic_numbers()
    assert len(z) == 83
    assert {1, 6, 7, 8, 9, 15, 16, 17, 35, 53} <= z
