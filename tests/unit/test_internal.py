"""T14：内坐标、Wilson B、线性键角与 proper_quartets 的必测用例。"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.internal import (
    InternalCoord,
    LinearAngleError,
    coord_value,
    linear_angles,
    proper_quartets,
    wilson_b,
)
from ffcompiler.importers import from_amber
from ffcompiler.ir import ForceFieldIR
from ffcompiler.ir.topology import bond_graph

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
FD_STEP_NM = 1e-6


def _load(name: str) -> tuple[ForceFieldIR, np.ndarray]:
    mdir = FIXTURES / name
    return from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )


def _ir_coords(ir: ForceFieldIR) -> list[InternalCoord]:
    """IR 全部 bond/angle + 键图全部 proper 四元组，组成内坐标集合。"""
    coords: list[InternalCoord] = [
        ("bond", (b.i, b.j)) for b in ir.bonds
    ]
    coords += [("angle", (a.i, a.j, a.k_atom)) for a in ir.angles]
    coords += [("dihedral", q) for q in proper_quartets(bond_graph(ir))]
    return coords


def _numeric_row(
    pos: np.ndarray, coord: InternalCoord, n3: int
) -> np.ndarray:
    """coord_value 的中心差分；二面角用 mod 2π 处理 ±π 处的跳变。"""
    row = np.zeros(n3)
    for a in range(n3):
        p_plus = pos.copy()
        p_plus[a // 3, a % 3] += FD_STEP_NM
        p_minus = pos.copy()
        p_minus[a // 3, a % 3] -= FD_STEP_NM
        v_plus = coord_value(p_plus, coord)
        v_minus = coord_value(p_minus, coord)
        diff = v_plus - v_minus
        if coord[0] == "dihedral":
            diff = (diff + math.pi) % (2.0 * math.pi) - math.pi
        row[a] = diff / (2.0 * FD_STEP_NM)
    return row


def test_wilson_b_matches_central_difference_ethanol() -> None:
    ir, pos = _load("ethanol")
    coords = _ir_coords(ir)
    n3 = 3 * pos.shape[0]
    b_matrix = wilson_b(pos, coords)
    assert b_matrix.shape == (len(coords), n3)
    for row_idx, coord in enumerate(coords):
        analytic = b_matrix[row_idx]
        numeric = _numeric_row(pos, coord, n3)
        scale = max(float(np.abs(numeric).max()), 1e-12)
        assert float(np.abs(analytic - numeric).max()) / scale <= 1e-7, coord


def test_wilson_b_translation_invariance_ethanol() -> None:
    ir, pos = _load("ethanol")
    b_matrix = wilson_b(pos, _ir_coords(ir))
    n = pos.shape[0]
    for row in b_matrix:
        axis_sums = row.reshape(n, 3).sum(axis=0)
        assert float(np.abs(axis_sums).max()) <= 1e-10


def test_acetonitrile_linear_angle() -> None:
    ir, pos = _load("acetonitrile")
    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]
    linear = linear_angles(pos, angles)
    assert len(linear) == 1
    i, j, k = linear[0]
    assert {ir.atoms[i].element, ir.atoms[j].element, ir.atoms[k].element} == {"C", "N"}
    assert ir.atoms[j].element == "C"  # 顶点是 sp 碳

    with pytest.raises(LinearAngleError) as exc_info:
        wilson_b(pos, [("angle", (i, j, k))])
    assert exc_info.value.angles == [(i, j, k)]

    # 经过该线性键角的二面角同样报错（H–C2–C1≡N）
    graph = bond_graph(ir)
    c_sp = j if ir.atoms[i].element == "C" else i  # 甲基碳（连着 H 的那个）
    methyl = i if c_sp == j else k
    substituent = next(a for a in graph[methyl] if a != c_sp)
    quartet = (substituent, methyl, c_sp, k) if methyl == i else (substituent, methyl, c_sp, i)
    # 统一为 (X, methyl, sp-C, N) 形式：顶点方向按实际 angle 顺序
    if (i, j, k) == (methyl, c_sp, k):
        quartet = (substituent, methyl, c_sp, k)
    else:
        quartet = (k, c_sp, methyl, substituent)
    with pytest.raises(LinearAngleError) as exc_info:
        wilson_b(pos, [("dihedral", quartet)])
    assert (i, j, k) in exc_info.value.angles


def test_proper_quartets_counts() -> None:
    # 正丁烷骨架：C1–C2–C3–C4 只有一个 proper 四元组
    butane = {
        0: frozenset({1}),
        1: frozenset({0, 2}),
        2: frozenset({1, 3}),
        3: frozenset({2}),
    }
    quartets = proper_quartets(butane)
    assert quartets == [(0, 1, 2, 3)]

    # 苯（含 H）：6 个 C 环 + 每个 C 一个 H，共 24 个；无重复、全部 i<l
    benzene = {a: set() for a in range(12)}
    for c in range(6):
        benzene[c].add((c + 1) % 6)
        benzene[(c + 1) % 6].add(c)
        benzene[c].add(6 + c)
        benzene[6 + c].add(c)
    benzene = {a: frozenset(nbrs) for a, nbrs in benzene.items()}
    quartets = proper_quartets(benzene)
    assert len(quartets) == 24
    assert len(set(quartets)) == 24
    assert all(i < l for (i, j, k, l) in quartets)


def test_dihedral_b_near_trans_branch() -> None:
    # 反式丁烷样几何：φ = π，δ 绕中心键小幅旋转跨过 ±π 分支
    r_j = np.array([0.0, 0.0, 0.0])
    r_k = np.array([0.154, 0.0, 0.0])
    r_i = np.array([-0.077, 0.12, 0.02])
    b1 = r_j - r_i
    r_l_trans = r_k + b1  # b3 = +b1 → φ = π（反式）

    for delta in (-0.01, 0.0, 0.01):
        c, s = math.cos(delta), math.sin(delta)
        rot = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
        r_l = r_k + rot @ (r_l_trans - r_k)
        pos = np.array([r_i, r_j, r_k, r_l])
        phi = coord_value(pos, ("dihedral", (0, 1, 2, 3)))
        if delta == 0.0:
            assert abs(abs(phi) - math.pi) <= 1e-12
        else:
            assert abs(abs(phi) - (math.pi - abs(delta))) <= 1e-12

        n3 = 12
        analytic = wilson_b(pos, [("dihedral", (0, 1, 2, 3))])[0]
        numeric = _numeric_row(pos, ("dihedral", (0, 1, 2, 3)), n3)
        scale = max(float(np.abs(numeric).max()), 1e-12)
        assert float(np.abs(analytic - numeric).max()) / scale <= 1e-7, delta
