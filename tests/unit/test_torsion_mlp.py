"""T28：MLPTorsionProvider（MVP-2 计划 §2.3–§2.7）。参考势用 MM（OpenMMCalculator），不加载 MACE。"""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.calculator import all_changes
from rdkit import Chem

from ffcompiler.importers import from_amber
from ffcompiler.ir import validate
from ffcompiler.ir.energy import dihedral_angle
from ffcompiler.torsion.mlp import MLPTorsionProvider
from ffcompiler.torsion.mm import OpenMMCalculator

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ALL_Z = frozenset(range(1, 119))


def _load(name: str):
    mdir = FIXTURES / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


def _torsion_energy(terms, pos) -> float:
    return sum(
        t.force_constant * (1 + math.cos(t.periodicity * dihedral_angle(*(pos[a] for a in (t.i, t.j, t.k, t.l))) - t.phase))
        for t in terms
    )


def _same_quartet(t, q) -> bool:
    return (t.i, t.j, t.k, t.l) in (tuple(q), tuple(q)[::-1])


def test_self_consistency_single_quartet_bond() -> None:
    ir, pos, mol = _load("ethanol")
    provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=ALL_Z)
    # 先拿到 C–O 键的扫描四元组，再把这根键上的项删到只剩这一个四元组
    from ffcompiler.torsion.rotors import find_rotors

    co = {r.bond: r for r in find_rotors(ir, mol, pos)}[(1, 2)]
    q = co.scan_quartet
    kept = [t for t in ir.proper_torsions if {t.j, t.k} != {1, 2} or _same_quartet(t, q)]
    orig = [t for t in kept if _same_quartet(t, q)]
    assert orig, "GAFF2 在扫描四元组上应当有项"
    toy = ir.replace_terms(proper_torsions=tuple(kept))

    new_ir, report = provider.parameterize(toy, pos, mol)
    row = {tuple(r["bond"]): r for r in report.rotors}[(1, 2)]
    assert row["status"] == "mlp"
    new_terms = [t for t in new_ir.proper_torsions if {t.j, t.k} == {1, 2}]
    assert len(new_terms) == 6 and all(t.source == "mlp" and _same_quartet(t, q) for t in new_terms)
    scan = report.scans[(1, 2)]
    diff = np.array([_torsion_energy(new_terms, x) - _torsion_energy(orig, x) for x in scan.positions])
    assert float(np.std(diff)) <= 1e-8
    validate(new_ir)


class _ChargeRecorder(OpenMMCalculator):
    charges: list = []

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        type(self).charges.append(atoms.info.get("charge"))
        super().calculate(atoms, properties, system_changes)


def test_charged_molecule_is_scanned_with_its_charge() -> None:
    ir, pos, mol = _load("acetate")
    _ChargeRecorder.charges = []
    new_ir, report = MLPTorsionProvider(lambda ir_: _ChargeRecorder(ir_), supported=ALL_Z).parameterize(ir, pos, mol)
    assert report.molecule_fallback is None
    assert [r["status"] for r in report.rotors] == ["mlp"]
    assert _ChargeRecorder.charges and set(_ChargeRecorder.charges) == {-1}


def test_odd_electrons_falls_back_whole_molecule() -> None:
    ir, pos, mol = _load("ethanol")
    atoms = list(ir.atoms)
    atoms[0] = replace(atoms[0], partial_charge=atoms[0].partial_charge + 1.0)  # 电荷守恒，IR 仍合法
    odd = replace(ir, atoms=tuple(atoms), total_charge=1)  # 26 − 1 = 25 个电子
    new_ir, report = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir), supported=ALL_Z).parameterize(odd, pos, mol)
    assert report.molecule_fallback == "odd_electrons"
    assert new_ir.proper_torsions == odd.proper_torsions
    assert all(r["status"] == "fallback" and r["reason"] == "odd_electrons" for r in report.rotors)


def test_unsupported_element_falls_back_whole_molecule() -> None:
    ir, pos, mol = _load("chlorobenzene")
    supported = frozenset({1, 6, 7, 8})
    new_ir, report = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=supported).parameterize(ir, pos, mol)
    assert report.molecule_fallback == "unsupported_element"
    assert new_ir.proper_torsions == ir.proper_torsions


class _Boom(OpenMMCalculator):
    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        raise RuntimeError("模型推理失败（测试注入）")


def test_error_on_second_rotor_only_falls_back_that_rotor() -> None:
    ir, pos, mol = _load("ethanol")
    calls = []

    def factory(ir_):
        calls.append(1)
        return OpenMMCalculator(ir_) if len(calls) == 1 else _Boom(ir_)

    new_ir, report = MLPTorsionProvider(factory, supported=ALL_Z).parameterize(ir, pos, mol)
    status = [(r["status"], r["reason"]) for r in report.rotors]
    assert status == [("mlp", None), ("fallback", "mlp_error")]
    second = tuple(report.rotors[1]["bond"])
    assert [t for t in new_ir.proper_torsions if {t.j, t.k} == set(second)] == \
           [t for t in ir.proper_torsions if {t.j, t.k} == set(second)]


def test_not_converged_falls_back() -> None:
    ir, pos, mol = _load("ethanol")
    provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=ALL_Z,
                                  scan_kwargs={"max_steps": 1})
    _, report = provider.parameterize(ir, pos, mol)
    assert {r["reason"] for r in report.rotors} == {"scan_not_converged"}


def test_dmso_both_rotors_parameterized() -> None:
    ir, pos, mol = _load("dmso")
    # 20° 网格（18 点 ≥ 13 个拟合参数），只为省单元测试时间（MVP-0 §1.5：not external ≤ 1 min）
    new_ir, report = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_), supported=ALL_Z,
                                        scan_kwargs={"step_deg": 20.0}).parameterize(ir, pos, mol)
    assert [r["status"] for r in report.rotors] == ["mlp", "mlp"]
    assert len({r["rotor_class"] for r in report.rotors}) == 1
    for r in report.rotors:
        terms = [t for t in new_ir.proper_torsions if {t.j, t.k} == set(r["bond"])]
        assert len(terms) == 6
    validate(new_ir)
    assert new_ir.provenance["torsion"]["rotors"][0]["status"] == "mlp"


def test_scan_cache_reused(tmp_path) -> None:
    ir, pos, mol = _load("ethanol")
    calls = []

    def factory(ir_):
        calls.append(1)
        return OpenMMCalculator(ir_)

    fast = {"step_deg": 20.0}
    p = MLPTorsionProvider(factory, supported=ALL_Z, scan_cache=tmp_path, scan_kwargs=fast)
    first, _ = p.parameterize(ir, pos, mol)
    n_first = len(calls)
    second, _ = p.parameterize(ir, pos, mol)
    assert len(calls) == n_first  # 第二次全部从缓存读
    assert first.proper_torsions == second.proper_torsions
    p2 = MLPTorsionProvider(factory, supported=ALL_Z, scan_cache=tmp_path, scan_kwargs={**fast, "fmax": 0.004})
    p2.parameterize(ir, pos, mol)
    assert len(calls) > n_first  # 扫描参数变了就重扫


def test_cache_key_includes_model_file_and_topology() -> None:
    from ffcompiler.torsion.mlp import scan_cache_key

    base = dict(model="M", quartet=(0, 1, 2, 3), info={"charge": 0, "spin": 1},
                scan_kwargs={"step_deg": 10.0}, model_file={"path": "/a", "size": 1, "mtime": 1.0},
                numbers=[6, 6, 8], bonds=[(0, 1), (1, 2)])
    k0 = scan_cache_key(**base)
    assert scan_cache_key(**{**base, "model_file": {"path": "/b", "size": 1, "mtime": 1.0}}) != k0
    assert scan_cache_key(**{**base, "numbers": [6, 8, 6]}) != k0
    assert scan_cache_key(**{**base, "bonds": [(0, 1), (0, 2)]}) != k0


def test_model_file_recorded_and_change_forces_rescan(tmp_path) -> None:
    ir, pos, mol = _load("ethanol")
    model = tmp_path / "fake.model"
    model.write_bytes(b"v1")
    calls = []

    def factory(ir_):
        calls.append(1)
        return OpenMMCalculator(ir_)

    fast = {"step_deg": 20.0}
    p = MLPTorsionProvider(factory, supported=ALL_Z, scan_cache=tmp_path / "c", scan_kwargs=fast,
                           model_path=model)
    new_ir, _ = p.parameterize(ir, pos, mol)
    prov = new_ir.provenance["torsion"]
    assert prov["model_file"]["path"] == str(model) and prov["model_file"]["size"] == 2
    assert {"mace-torch", "ase"} <= set(prov["package_versions"])
    n = len(calls)
    model.write_bytes(b"v2-longer")  # 换了模型文件：缓存不能复用
    p.parameterize(ir, pos, mol)
    assert len(calls) > n


def test_cache_only_missing_scan_is_recorded_fallback(tmp_path) -> None:
    ir, pos, mol = _load("ethanol")
    calls = []
    p = MLPTorsionProvider(lambda ir_: calls.append(1) or OpenMMCalculator(ir_), supported=ALL_Z,
                           scan_cache=tmp_path, cache_only=True)
    new_ir, report = p.parameterize(ir, pos, mol)
    assert calls == []  # 不扫描
    assert {(r["status"], r["reason"]) for r in report.rotors} == {("fallback", "no_cached_scan")}
    assert new_ir.proper_torsions == ir.proper_torsions


def test_scan_cache_key_without_reference_is_unchanged() -> None:
    import json

    from ffcompiler.torsion.mlp import scan_cache_key

    # 保护已有 MACE 缓存（jobs/gate_mvp2_20261008_3/scans）：reference=None 时键与旧公式逐字相同
    kw = dict(model_file={"path": "m", "size": 1, "mtime": 2.0}, numbers=[6, 1], bonds=[(1, 0)])
    old = json.dumps({"model": "M", "model_file": {"path": "m", "size": 1, "mtime": 2.0},
                      "quartet": [0, 1, 2, 3], "info": {"charge": 0, "spin": 1},
                      "scan": {"step_deg": 10.0}, "numbers": [6, 1], "bonds": [[0, 1]]}, sort_keys=True)
    assert scan_cache_key("M", (0, 1, 2, 3), {"charge": 0, "spin": 1}, {"step_deg": 10.0}, **kw) == old
    with_ref = scan_cache_key("M", (0, 1, 2, 3), {"charge": 0, "spin": 1}, {"step_deg": 10.0}, **kw,
                              reference={"kind": "mm-self", "charges": [0.1, -0.1]})
    assert json.loads(with_ref)["reference"] == {"kind": "mm-self", "charges": [0.1, -0.1]}
