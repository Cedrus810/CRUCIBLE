"""T34：NAGLCharges（MVP-3 计划 §2.3；import torch，标 external）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeError
from ffcompiler.charge.loaded import read_mol2_atoms
from ffcompiler.charge.nagl import NAGLCharges, nagl_raw_charges
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
def test_all_fixtures_conserve_charge_and_stay_near_am1bcc(name: str) -> None:
    ir, pos, mol = _load(name)
    new_ir, report = assign_charges(ir, NAGLCharges(), pos, mol)
    got = np.array([a.partial_charge for a in new_ir.atoms])
    assert abs(got.sum() - ir.total_charge) <= 1e-10
    ref = np.array([q for _, q in read_mol2_atoms(FIXTURES / name / f"{name}.mol2")])
    assert np.max(np.abs(got - ref)) <= 0.1  # 只抓原子错位（冒烟实测最大 0.040 e）
    src = report.source
    assert src["model"] == "openff-gnn-am1bcc-1.0.0.pt" and src["size"] > 0 and src["openff_nagl"]


def test_out_of_domain_raises() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("C[Si](C)(C)C"))
    AllChem.EmbedMolecule(mol, randomSeed=20261004)
    with pytest.raises(ChargeError):
        nagl_raw_charges(mol, "openff-gnn-am1bcc-1.0.0.pt")


def test_reordered_mol_raises() -> None:
    ir, pos, mol = _load("ethanol")
    order = list(range(mol.GetNumAtoms()))[::-1]
    with pytest.raises(TopologyMismatchError):
        NAGLCharges().compute(ir, pos, Chem.RenumberAtoms(mol, order))


def test_swapped_hydrogens_raise() -> None:
    # 审阅 Important 1：甲基 H（3）与羟基 H（8）对调，元素序列不变但键集合变了
    ir, pos, mol = _load("ethanol")
    order = list(range(mol.GetNumAtoms()))
    order[3], order[8] = order[8], order[3]
    with pytest.raises(TopologyMismatchError):
        NAGLCharges().compute(ir, pos, Chem.RenumberAtoms(mol, order))
