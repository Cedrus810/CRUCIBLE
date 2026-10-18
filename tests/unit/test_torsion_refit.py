"""T35：二面角 replaced_terms 与换电荷后的重拟合（MVP-3 计划 §2.4）。参考势用 MM，不加载 MACE。"""

from __future__ import annotations

import copy
import shutil
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeProvider
from ffcompiler.charge.loaded import LoadedCharges
from ffcompiler.importers import from_amber
from ffcompiler.torsion.mlp import MLPTorsionProvider, TorsionRefitError
from ffcompiler.torsion.mm import OpenMMCalculator
from ffcompiler.torsion.scan import ScanError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ALL_Z = frozenset(range(1, 119))
FAST = {"step_deg": 20.0}  # 18 点 ≥ 13 个拟合参数；只为省单元测试时间


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    """整个模块只扫描一次（乙醇两根键）；各测试用自己的缓存副本。"""
    cache = tmp_path_factory.mktemp("scans")
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    base, _ = assign_charges(ir, LoadedCharges(), pos, mol)  # 先均摊，之后再均摊是 no-op
    provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=ALL_Z,
                                  scan_kwargs=FAST, scan_cache=cache)
    mlp_ir, report = provider.parameterize(base, pos, mol)
    return base, mlp_ir, report, pos, mol, cache


def _mlp_ethanol(scanned, tmp_path: Path):
    base, mlp_ir, report, pos, mol, cache = scanned
    own = tmp_path / "scans"
    shutil.copytree(cache, own)
    return base, mlp_ir, report, pos, mol, own


def test_replaced_terms_recorded(scanned, tmp_path: Path) -> None:
    base, mlp_ir, report, _, _, _ = _mlp_ethanol(scanned, tmp_path)
    rows = [r for r in report.rotors if r["status"] == "mlp"]
    assert len(rows) == 2
    for row in rows:
        bond = set(row["bond"])
        expected = [(i, t) for i, t in enumerate(base.proper_torsions) if {t.j, t.k} == bond]
        got = row["replaced_terms"]
        assert [r["index"] for r in got] == [i for i, _ in expected]
        assert [(r["i"], r["j"], r["k"], r["l"], r["periodicity"], r["phase"], r["force_constant"])
                for r in got] == [(t.i, t.j, t.k, t.l, t.periodicity, t.phase, t.force_constant)
                                  for _, t in expected]
    assert mlp_ir.provenance["torsion"]["charge_method"] == "loaded"


def test_refit_with_unchanged_charges_is_identical_and_idempotent(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    again, report = assign_charges(mlp_ir, LoadedCharges(), pos, mol, scan_cache=cache)
    assert again.proper_torsions == mlp_ir.proper_torsions
    assert report.torsion_refit is not None
    third, _ = assign_charges(again, LoadedCharges(), pos, mol, scan_cache=cache)
    assert third.proper_torsions == mlp_ir.proper_torsions
    assert [a.partial_charge for a in third.atoms] == [a.partial_charge for a in mlp_ir.atoms]


class _Shifted(ChargeProvider):
    method = "shifted-test"

    def __init__(self, q):
        self.q = q

    def compute(self, ir, positions_nm, mol):
        return self.q.copy()

    def source_info(self):
        return {"kind": "test"}


def test_refit_with_new_charges_changes_terms(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    assert mlp_ir.provenance["torsion"]["charge_method"] == "loaded"
    q = np.array([a.partial_charge for a in mlp_ir.atoms])
    q[0] += 0.05
    q[1] -= 0.05
    new_ir, report = assign_charges(mlp_ir, _Shifted(q), pos, mol, scan_cache=cache)
    rows = new_ir.provenance["torsion"]["rotors"]
    assert [r["status"] for r in rows] == ["mlp", "mlp"]
    assert new_ir.proper_torsions != mlp_ir.proper_torsions
    assert new_ir.provenance["torsion"]["charge_method"] == "shifted-test"
    assert report.torsion_refit["rotors"] == rows


def test_refit_without_cache_raises(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    with pytest.raises(TorsionRefitError):
        assign_charges(mlp_ir, LoadedCharges(), pos, mol)


def test_refit_with_missing_scan_file_names_the_bond(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    victim = sorted(cache.glob("*.npz"))[0]
    victim.unlink()
    with pytest.raises(TorsionRefitError, match=victim.stem):
        assign_charges(mlp_ir, LoadedCharges(), pos, mol, scan_cache=cache)


def test_refit_without_replaced_terms_raises(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    torsion_prov = copy.deepcopy(dict(mlp_ir.provenance["torsion"]))
    for row in torsion_prov["rotors"]:
        row.pop("replaced_terms", None)
    stripped = mlp_ir.replace_terms(provenance_update={"torsion": torsion_prov})
    with pytest.raises(TorsionRefitError, match="derive-torsions"):
        assign_charges(stripped, LoadedCharges(), pos, mol, scan_cache=cache)


def test_refit_keeps_original_fallback_reason_and_cache_bytes(tmp_path: Path) -> None:
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    base, _ = assign_charges(ir, LoadedCharges(), pos, mol)
    calls = []

    def factory(ir_):
        calls.append(1)
        if len(calls) == 2:
            raise ScanError("boom")  # 第二根键模拟参考势出错 → mlp_error
        return OpenMMCalculator(ir_)

    cache = tmp_path / "scans"
    mlp_ir, report = MLPTorsionProvider(factory, supported=ALL_Z, scan_kwargs=FAST,
                                        scan_cache=cache).parameterize(base, pos, mol)
    reasons = {tuple(r["bond"]): r["reason"] for r in report.rotors}
    assert sorted(r["status"] for r in report.rotors) == ["fallback", "mlp"]
    failed = next(b for b, why in reasons.items() if why == "mlp_error")
    before = {p.name: p.read_bytes() for p in cache.iterdir()}

    new_ir, chg = assign_charges(mlp_ir, LoadedCharges(), pos, mol, scan_cache=cache)
    rows = {tuple(r["bond"]): r for r in chg.torsion_refit["rotors"]}
    assert rows[failed]["reason"] == "mlp_error"
    assert {tuple(r["bond"]): r["reason"] for r in new_ir.provenance["torsion"]["rotors"]}[failed] == "mlp_error"
    assert {p.name: p.read_bytes() for p in cache.iterdir()} == before


def test_refit_keeps_original_fallback_even_if_cache_gains_scan(tmp_path: Path) -> None:
    """原兜底键之后缓存里补齐了匹配扫描，重拟合仍原样沿用原行与原项，不升级为 mlp（2026-10-13 用户定）。"""
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    base, _ = assign_charges(ir, LoadedCharges(), pos, mol)
    calls = []

    def factory(ir_):
        calls.append(1)
        if len(calls) == 2:
            raise ScanError("boom")
        return OpenMMCalculator(ir_)

    cache = tmp_path / "scans"
    mlp_ir, report = MLPTorsionProvider(factory, supported=ALL_Z, scan_kwargs=FAST,
                                        scan_cache=cache).parameterize(base, pos, mol)
    failed = next(tuple(r["bond"]) for r in report.rotors if r["reason"] == "mlp_error")
    # 补扫原兜底键：缓存里出现与重拟合缓存键一致的扫描
    MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=ALL_Z, scan_kwargs=FAST,
                       scan_cache=cache).parameterize(base, pos, mol)
    assert (cache / f"{failed[0]}-{failed[1]}.npz").exists()

    def on_bond(ir_, bond):
        return sorted((t.i, t.j, t.k, t.l, t.periodicity, t.force_constant, t.phase, t.source)
                      for t in ir_.proper_torsions if {t.j, t.k} == set(bond))

    new_ir, chg = assign_charges(mlp_ir, LoadedCharges(), pos, mol, scan_cache=cache)
    rows = {tuple(r["bond"]): r for r in new_ir.provenance["torsion"]["rotors"]}
    assert rows[failed]["status"] == "fallback" and rows[failed]["reason"] == "mlp_error"
    assert on_bond(new_ir, failed) == on_bond(mlp_ir, failed)
    assert {tuple(r["bond"]): r["status"] for r in chg.torsion_refit["rotors"]}[failed] == "fallback"
