"""molecule/io.py 的 §3.1 入口几何合同（MVP-6，2026-11-11 计划 §2）。"""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def test_read_sdf_valid_ethanol():
    from ffcompiler.molecule.io import read_sdf

    mi = read_sdf(FIXTURES / "ethanol" / "ethanol.sdf")
    assert mi.n_atoms == 9
    assert mi.net_charge == 0
    assert mi.name == "ethanol"
    assert len(mi.sha256) == 64


def test_read_sdf_net_charge_from_sdf():
    from ffcompiler.molecule.io import read_sdf

    assert read_sdf(FIXTURES / "acetate" / "acetate.sdf").net_charge == -1
    assert read_sdf(FIXTURES / "methylammonium" / "methylammonium.sdf").net_charge == 1


def test_read_sdf_missing_file():
    from ffcompiler.molecule.io import MoleculeInputError, read_sdf

    with pytest.raises(MoleculeInputError, match="不存在"):
        read_sdf(FIXTURES / "ethanol" / "nope.sdf")


def test_read_sdf_unparsable(tmp_path):
    from ffcompiler.molecule.io import MoleculeInputError, read_sdf

    bad = tmp_path / "bad.sdf"
    bad.write_text("not an sdf at all\n", encoding="utf-8")
    with pytest.raises(MoleculeInputError, match="解析失败"):
        read_sdf(bad)


def test_read_sdf_rejects_no_3d(tmp_path):
    """无构象/全零坐标的 SDF 必须被拒（§3.1：不静默生成构象）。"""
    from rdkit import Chem

    from ffcompiler.molecule.io import MoleculeInputError, read_sdf

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))  # 不嵌入构象 → 无 3D
    sdf = tmp_path / "flat.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    with pytest.raises(MoleculeInputError, match="3D"):
        read_sdf(sdf)


def test_read_sdf_charge_mismatch():
    from ffcompiler.molecule.io import MoleculeInputError, read_sdf

    with pytest.raises(MoleculeInputError, match="不一致"):
        read_sdf(FIXTURES / "acetate" / "acetate.sdf", net_charge=0)


def test_read_sdf_rejects_clash(tmp_path):
    """任意两原子 < 0.5 Å 视为严重穿距（§3.1）。"""
    from rdkit import Chem

    from ffcompiler.molecule.io import CLASH_MIN_A, MoleculeInputError, read_sdf

    mol = Chem.MolFromMolFile(str(FIXTURES / "ethanol" / "ethanol.sdf"), removeHs=False)
    conf = mol.GetConformer()
    conf.SetAtomPosition(5, conf.GetAtomPosition(4))  # 两个原子叠在一起
    sdf = tmp_path / "clash.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    with pytest.raises(MoleculeInputError, match=f"{CLASH_MIN_A} Å"):
        read_sdf(sdf)
