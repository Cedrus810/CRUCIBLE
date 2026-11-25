"""覆盖集清单 benchmarks/coverage_set/molecules.json（计划 §2.1，用户 2026-11-24 认可 34 个）。"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from rdkit import Chem

from ffcompiler.validation.coverage import load_table

TABLE = Path(__file__).resolve().parents[2] / "benchmarks" / "coverage_set" / "molecules.json"

# 计划 §2.1 表：每类（已有 + 新增）的分子数
PLAN_COUNTS = {
    "alcohol_ether": 4, "carbonyl_ester_amide": 4, "amine": 2, "heteroaromatic": 3,
    "aromatic_conjugated_fused": 6, "small_ring": 3, "heteroatom_halogen": 9,
    "ion_zwitterion": 3,
}


def test_table_has_34_molecules_matching_plan_categories():
    rows = load_table(TABLE)
    assert len(rows) == 34
    assert Counter(r["category"] for r in rows) == PLAN_COUNTS


def test_every_sdf_exists_and_formal_charge_matches():
    for row in load_table(TABLE):
        mol = Chem.MolFromMolFile(str(row["sdf"]), removeHs=False)
        assert mol is not None, row["name"]
        assert mol.GetNumConformers() == 1 and mol.GetConformer().Is3D(), row["name"]
        assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == row["net_charge"], row["name"]


def test_fixtures_are_referenced_not_copied():
    rows = {r["name"]: r for r in load_table(TABLE)}
    fixtures = [n for n, r in rows.items() if "tests/fixtures/gaff2" in str(r["sdf"])]
    assert len(fixtures) == 12
