"""可转键识别、扫描四元组与等价类（MVP-2 计划 §2.1，T25）。

键 (j,k) 可转 ⟺ RDKit 键类型 SINGLE、不在环内、两端度数（含氢）都 ≥ 2、
两端都不是线性中心（以它为顶点的任一键角 > 175°，与 MVP-1 §2.4 相同）。
等价类只记录、不共享参数（T25 裁定：图自同构在几何上可能是镜面）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from ..bonded.assign import TopologyMismatchError
from ..bonded.symmetry import term_classes
from ..ir.schema import ForceFieldIR
from ..ir.topology import bond_graph

LINEAR_ANGLE_DEG = 175.0


@dataclass(frozen=True)
class Rotor:
    bond: tuple[int, int]                              # j < k
    scan_quartet: tuple[int, int, int, int]            # (i, j, k, l)
    ir_quartets: tuple[tuple[int, int, int, int], ...]  # IR 里中心键为 {j,k} 的不同四元组
    rotor_class: int


def check_mol_topology(ir: ForceFieldIR, mol: Chem.Mol) -> None:
    """SDF 分子与 IR：原子数、逐原子原子序数、键集合都一致，否则 TopologyMismatchError。"""
    if mol.GetNumAtoms() != len(ir.atoms):
        raise TopologyMismatchError(
            f"SDF 有 {mol.GetNumAtoms()} 个原子，IR 有 {len(ir.atoms)} 个")
    for atom in ir.atoms:
        z = mol.GetAtomWithIdx(atom.index).GetAtomicNum()
        if z != atom.atomic_number:
            raise TopologyMismatchError(
                f"原子 {atom.index}：SDF 原子序数 {z}，IR {atom.atomic_number}")
    ir_bonds = {frozenset((b.i, b.j)) for b in ir.bonds}
    mol_bonds = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()}
    if ir_bonds != mol_bonds:
        only_ir = sorted(tuple(sorted(b)) for b in ir_bonds - mol_bonds)
        only_mol = sorted(tuple(sorted(b)) for b in mol_bonds - ir_bonds)
        raise TopologyMismatchError(f"键集合不一致：只在 IR {only_ir}，只在 SDF {only_mol}")


def _is_linear_center(center: int, graph, pos: np.ndarray) -> bool:
    nbrs = sorted(graph[center])
    for a in range(len(nbrs)):
        for b in range(a + 1, len(nbrs)):
            u = pos[nbrs[a]] - pos[center]
            v = pos[nbrs[b]] - pos[center]
            theta = math.degrees(math.atan2(np.linalg.norm(np.cross(u, v)), float(np.dot(u, v))))
            if theta > LINEAR_ANGLE_DEG:
                return True
    return False


def _scan_quartet(j: int, k: int, graph, heavy: set[int]) -> tuple[int, int, int, int]:
    candidates = [
        (i, j, k, l)
        for i in graph[j] - {k}
        for l in graph[k] - {j}
        if i != l
    ]
    return min(candidates, key=lambda q: (-((q[0] in heavy) + (q[3] in heavy)), q[0], q[3]))


def find_rotors(ir: ForceFieldIR, mol: Chem.Mol, positions_nm: np.ndarray) -> list[Rotor]:
    check_mol_topology(ir, mol)
    pos = np.asarray(positions_nm, dtype=np.float64)
    graph = bond_graph(ir)
    heavy = {a.index for a in ir.atoms if a.atomic_number > 1}
    found: list[tuple[tuple[int, int], tuple[int, int, int, int], tuple]] = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        j, k = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        if len(graph[j]) < 2 or len(graph[k]) < 2:
            continue
        if _is_linear_center(j, graph, pos) or _is_linear_center(k, graph, pos):
            continue
        quartets = sorted({
            (t.i, t.j, t.k, t.l) for t in ir.proper_torsions if {t.j, t.k} == {j, k}
        })
        found.append(((j, k), _scan_quartet(j, k, graph, heavy), tuple(quartets)))
    found.sort()
    classes = term_classes(ir, [q for _, q, _ in found]) if found else ()
    return [
        Rotor(bond=b, scan_quartet=q, ir_quartets=irq, rotor_class=int(c))
        for (b, q, irq), c in zip(found, classes)
    ]
