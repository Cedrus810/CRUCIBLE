"""T25：可转键识别（MVP-2 计划 §2.1）。"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from rdkit import Chem

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.importers import from_amber
from ffcompiler.torsion.rotors import find_rotors

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"

EXPECTED = {
    "ethanol": 2, "dmso": 2, "nma": 3, "methyl_benzoate": 3, "acetate": 1, "methylammonium": 1,
    "benzene": 0, "chlorobenzene": 0, "cyclopropane": 0, "cyclobutane": 0, "norbornane": 0,
    "acetonitrile": 0,
}


def _load(name: str):
    mdir = FIXTURES / name
    ir, pos = from_amber(mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf")
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_rotor_counts(name: str) -> None:
    ir, pos, mol = _load(name)
    rotors = find_rotors(ir, mol, pos)
    assert len(rotors) == EXPECTED[name]
    graph = {a.index: set() for a in ir.atoms}
    for b in ir.bonds:
        graph[b.i].add(b.j)
        graph[b.j].add(b.i)
    for r in rotors:
        j, k = r.bond
        assert j < k and k in graph[j]
        i, j2, k2, l = r.scan_quartet
        assert (j2, k2) == (j, k) and i in graph[j] - {k} and l in graph[k] - {j}
        central = {(t.j, t.k) for t in ir.proper_torsions} | {(t.k, t.j) for t in ir.proper_torsions}
        assert all(tuple(q[1:3]) in central for q in r.ir_quartets)
        assert all({q[1], q[2]} == {j, k} for q in r.ir_quartets)


def test_scan_quartet_prefers_heavy_ends() -> None:
    ir, pos, mol = _load("ethanol")
    heavy = {a.index for a in ir.atoms if a.atomic_number > 1}
    by_bond = {r.bond: r for r in find_rotors(ir, mol, pos)}
    cc = by_bond[(0, 1)]
    co = by_bond[(1, 2)]
    # C-O：i 可以选 C(0)；l 只能是羟基 H
    assert co.scan_quartet[0] in heavy
    # C-C：甲基端只有 H，O 端选 O
    assert cc.scan_quartet[3] in heavy


def test_dmso_two_rotors_same_class_nma_methyls_distinct() -> None:
    ir, pos, mol = _load("dmso")
    rotors = find_rotors(ir, mol, pos)
    assert len({r.rotor_class for r in rotors}) == 1
    ir, pos, mol = _load("nma")
    rotors = find_rotors(ir, mol, pos)
    assert len({r.rotor_class for r in rotors}) == 3


def test_atomic_number_mismatch_raises() -> None:
    ir, pos, mol = _load("ethanol")
    rw = Chem.RWMol(mol)
    rw.GetAtomWithIdx(2).SetAtomicNum(16)
    with pytest.raises(TopologyMismatchError):
        find_rotors(ir, rw.GetMol(), pos)


def test_missing_bond_raises() -> None:
    ir, pos, mol = _load("ethanol")
    short = replace(ir, bonds=ir.bonds[:-1])
    with pytest.raises(TopologyMismatchError):
        find_rotors(short, mol, pos)
