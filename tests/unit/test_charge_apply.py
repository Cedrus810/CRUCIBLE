"""T31：电荷写回框架（MVP-3 计划 §2.1、§2.2）。"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeError, ChargeProvider, equalize_residual
from ffcompiler.charge.loaded import LoadedCharges
from ffcompiler.importers import from_amber
from ffcompiler.ir import validate

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _load(name: str):
    mdir = FIXTURES / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


class _Fixed(ChargeProvider):
    method = "fixed-test"

    def __init__(self, q):
        self.q = np.asarray(q, dtype=float)

    def compute(self, ir, positions_nm, mol):
        return self.q.copy()

    def source_info(self):
        return {"kind": "test"}


def test_equalize_residual_spreads_and_last_atom_absorbs() -> None:
    q = np.array([0.1, -0.2, 0.3, -0.1, 0.05, -0.05, 0.2, -0.299, 0.0])  # Σ = 1e-3
    out, r = equalize_residual(q, 0)
    assert r == pytest.approx(-1e-3, abs=1e-15)
    assert abs(out.sum()) <= 1e-15
    np.testing.assert_array_equal(out[:-1], q[:-1] + r / len(q))


def test_equalize_residual_no_op_below_1e12() -> None:
    q = np.array([0.5, -0.5 + 5e-13])
    out, _ = equalize_residual(q, 0)
    np.testing.assert_array_equal(out, q)


def test_assign_writes_only_partial_charges_and_provenance() -> None:
    ir, pos, mol = _load("ethanol")
    q = np.linspace(-0.4, 0.4, len(ir.atoms))
    q[-1] += 1e-3
    new_ir, report = assign_charges(ir, _Fixed(q), pos, mol)
    expected, r = equalize_residual(q, ir.total_charge)
    assert [a.partial_charge for a in new_ir.atoms] == list(expected)
    assert [replace(a, partial_charge=0.0) for a in new_ir.atoms] == \
        [replace(a, partial_charge=0.0) for a in ir.atoms]
    for field in ("bonds", "angles", "proper_torsions", "impropers", "total_charge"):
        assert getattr(new_ir, field) == getattr(ir, field)
    validate(new_ir)
    prov = new_ir.provenance["charge"]
    assert prov["method"] == "fixed-test" and prov["residual"] == r
    assert prov["raw_sum"] == pytest.approx(q.sum()) and prov["source"] == {"kind": "test"}
    assert report.method == "fixed-test" and report.torsion_refit is None


@pytest.mark.parametrize("name", ["acetate", "methylammonium"])
def test_charged_molecules_conserve_total_charge(name: str) -> None:
    ir, pos, mol = _load(name)
    n = len(ir.atoms)
    q = np.full(n, ir.total_charge / n) + np.linspace(-0.1, 0.1, n)
    q[0] += 2e-7
    new_ir, _ = assign_charges(ir, _Fixed(q), pos, mol)
    assert abs(sum(a.partial_charge for a in new_ir.atoms) - ir.total_charge) <= 1e-10


@pytest.mark.parametrize("delta", [1, -1])
def test_wrong_length_from_provider_raises(delta: int) -> None:
    ir, pos, mol = _load("ethanol")
    with pytest.raises(ChargeError, match="形状"):
        assign_charges(ir, _Fixed(np.zeros(len(ir.atoms) + delta)), pos, mol)


def test_non_finite_charge_raises() -> None:
    ir, pos, mol = _load("ethanol")
    q = np.zeros(len(ir.atoms))
    q[3] = np.nan
    with pytest.raises(ChargeError, match="非有限"):
        assign_charges(ir, _Fixed(q), pos, mol)


def test_reloading_keeps_previous_charge_source() -> None:
    ir, pos, mol = _load("ethanol")
    assert "charge" not in ir.provenance
    first, _ = assign_charges(ir, _Fixed(np.zeros(len(ir.atoms))), pos, mol)
    assert first.provenance["charge"]["previous"] is None
    again, _ = assign_charges(first, LoadedCharges(), pos, mol)
    prov = again.provenance["charge"]
    assert prov["method"] == "loaded"
    assert prov["previous"]["method"] == "fixed-test"
    assert prov["previous"]["source"] == {"kind": "test"}
