"""T33：AM1BCCCharges（MVP-3 计划 §2.3；调用 antechamber，标 external）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.charge.am1bcc import AM1BCCCharges
from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeError
from ffcompiler.charge.loaded import read_mol2_atoms
from ffcompiler.importers import from_amber

pytestmark = pytest.mark.external

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ALL12 = ("ethanol", "dmso", "nma", "methyl_benzoate", "acetate", "methylammonium", "benzene",
         "chlorobenzene", "cyclopropane", "cyclobutane", "norbornane", "acetonitrile")


def _load(name: str):
    mdir = FIXTURES / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


@pytest.mark.parametrize("name", ALL12)
def test_reproduces_fixture_charges(name: str, tmp_path: Path) -> None:
    ir, pos, mol = _load(name)
    new_ir, report = assign_charges(ir, AM1BCCCharges(tmp_path), pos, mol)
    ref = np.array([q for _, q in read_mol2_atoms(FIXTURES / name / f"{name}.mol2")])
    got = np.array([a.partial_charge for a in new_ir.atoms])
    assert np.max(np.abs(got - ref)) <= 1e-10
    assert report.method == "am1bcc" and report.source["antechamber"]


@pytest.mark.parametrize("name", ["acetate", "methylammonium"])
def test_net_charge_flag_is_ir_total_charge(name: str, tmp_path: Path) -> None:
    ir, pos, mol = _load(name)
    AM1BCCCharges(tmp_path).compute(ir, pos, mol)
    argv = (tmp_path / "command.txt").read_text().split()
    assert argv[argv.index("-nc") + 1] == str(ir.total_charge)


def test_failing_antechamber_raises_with_log(tmp_path: Path, monkeypatch) -> None:
    fake = tmp_path / "fake_antechamber.sh"
    fake.write_text("#!/bin/bash\necho 'boom from fake antechamber'\nexit 3\n")
    fake.chmod(0o755)
    monkeypatch.setenv("FFC_ANTECHAMBER", str(fake))
    ir, pos, mol = _load("ethanol")
    with pytest.raises(ChargeError, match="boom from fake antechamber"):
        AM1BCCCharges(tmp_path / "work").compute(ir, pos, mol)


def test_swapped_hydrogens_raise_before_running_antechamber(tmp_path: Path) -> None:
    # 审阅 Important 1：元素序列相同、键集合不同的 SDF 必须在调用 antechamber 之前被拒
    from ffcompiler.bonded.assign import TopologyMismatchError

    ir, pos, mol = _load("ethanol")
    order = list(range(mol.GetNumAtoms()))
    order[3], order[8] = order[8], order[3]
    with pytest.raises(TopologyMismatchError):
        AM1BCCCharges(tmp_path).compute(ir, pos, Chem.RenumberAtoms(mol, order))
    assert not (tmp_path / "command.txt").exists()
