"""T32：LoadedCharges（MVP-3 计划 §2.3）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.charge.loaded import LoadedCharges
from ffcompiler.importers import from_amber

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _load(name: str):
    mdir = FIXTURES / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


def _ir_q(ir):
    return np.array([a.partial_charge for a in ir.atoms])


def test_none_returns_current_ir_charges() -> None:
    ir, pos, mol = _load("ethanol")
    np.testing.assert_array_equal(LoadedCharges().compute(ir, pos, mol), _ir_q(ir))
    assert LoadedCharges().method == "loaded" and LoadedCharges().source_info() == {"source": "ir"}


def test_fixture_mol2_matches_prmtop_charges() -> None:
    ir, pos, mol = _load("ethanol")
    q = LoadedCharges(FIXTURES / "ethanol" / "ethanol.mol2").compute(ir, pos, mol)
    assert np.max(np.abs(q - _ir_q(ir))) <= 1e-8


def test_text_file_with_blank_and_comment_lines(tmp_path: Path) -> None:
    ir, pos, mol = _load("ethanol")
    q = np.linspace(-0.3, 0.3, len(ir.atoms))
    path = tmp_path / "q.txt"
    path.write_text("# charges\n\n" + "\n".join(f"{x:.12f}" for x in q) + "\n\n")
    np.testing.assert_allclose(LoadedCharges(path).compute(ir, pos, mol), q, atol=1e-12)
    assert LoadedCharges(path).source_info() == {"source": str(path)}


def test_short_text_file_raises(tmp_path: Path) -> None:
    ir, pos, mol = _load("ethanol")
    path = tmp_path / "q.txt"
    path.write_text("\n".join("0.0" for _ in range(len(ir.atoms) - 1)))
    with pytest.raises(TopologyMismatchError):
        LoadedCharges(path).compute(ir, pos, mol)


def test_mol2_with_renamed_atom_raises(tmp_path: Path) -> None:
    ir, pos, mol = _load("ethanol")
    text = (FIXTURES / "ethanol" / "ethanol.mol2").read_text()
    bad = tmp_path / "bad.mol2"
    bad.write_text(text.replace(f" {ir.atoms[0].name} ", " XX9 ", 1))
    with pytest.raises(TopologyMismatchError):
        LoadedCharges(bad).compute(ir, pos, mol)
