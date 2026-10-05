"""DRIH：冗余内坐标 Hessian 投影（MVP-1 计划 T19，§2.5 锁定）。

冗余坐标集合 = bond ∪ 非线性 angle ∪ proper dihedral（不含经过线性键角的
四元组）。G = BBᵀ 伪逆（特征值 < 1e-8·λmax 的方向置零），A = G⁻B，
H_q = A·H·Aᵀ；bond/angle 的 k 取 H_q 对角元，平衡值取参考几何上的坐标值。
忽略梯度项 → 必须先通过驻点检查（check_hessian）。

**内部单位：原子单位**（决策 0012，用户已确认）。几何换成 Bohr、Cartesian
Hessian 换成 E_h/a₀² 后再做 Wilson B 变换；H_q 的单位随坐标类型而定：bond 对角元
E_h/a₀²（×HARTREE/BOHR² → kJ/mol/nm²），angle 对角元 E_h/rad²（×HARTREE →
kJ/mol/rad²，不再除 a₀²）。冗余关系同时牵涉 bond 与 angle 时（环系），G 的伪逆
依赖长度单位：Bohr 与 Sobtop DRIH 一致，nm/Å 不一致。能量单位不影响结果。
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ffcompiler.bonded.internal import (
    InternalCoord,
    LinearAngleError,
    angle_value,
    bond_value,
    linear_angles,
    wilson_b,
)
from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.hessian.base import CartesianHessian, check_hessian
from ffcompiler.units import BOHR_NM, HARTREE_KJ_MOL

_EIGENVALUE_CUTOFF = 1e-8
# 原子单位（决策 0012）：E_h/a₀² ↔ kJ mol⁻¹ nm⁻²，E_h/rad² ↔ kJ mol⁻¹ rad⁻²
_STRETCH_AU_TO_IR = HARTREE_KJ_MOL / BOHR_NM**2
_ANGLE_AU_TO_IR = HARTREE_KJ_MOL


def _pseudo_inverse(g: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(g)
    cutoff = _EIGENVALUE_CUTOFF * float(eigenvalues.max())
    inverse = np.where(eigenvalues >= cutoff, 1.0 / np.where(eigenvalues >= cutoff, eigenvalues, 1.0), 0.0)
    return (eigenvectors * inverse) @ eigenvectors.T


def drih(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
    *,
    dihedrals: Sequence[tuple[int, int, int, int]],
    allow_nonstationary: bool = False,
) -> list[DerivedTerm]:
    """DRIH 推导 bond/angle 参数（½k 约定；先 bond 后 angle，按输入顺序）。

    输入含线性键角 → LinearAngleError（§2.4，由调用方按 linear_angle_policy
    先剔除）。经过线性键角的 dihedral 四元组按几何判定后剔除（§2.5 第 1 条），
    不依赖 angle 列表里是否还有该线性角。`dihedrals` 只用来组成冗余坐标
    集合，不返回。
    """
    check_hessian(h, allow_nonstationary=allow_nonstationary)
    pos = h.positions_nm
    linear = linear_angles(pos, list(angles))
    if linear:
        raise LinearAngleError(linear)
    kept_dihedrals = [
        q for q in dihedrals
        if not linear_angles(pos, [(q[0], q[1], q[2]), (q[1], q[2], q[3])])
    ]

    coords: list[InternalCoord] = [("bond", (i, j)) for i, j in bonds]
    coords += [("angle", a) for a in angles]
    coords += [("dihedral", q) for q in kept_dihedrals]

    pos_bohr = pos / BOHR_NM
    h_au = 0.5 * (h.hessian + h.hessian.T) / _STRETCH_AU_TO_IR  # E_h/a₀²
    b_matrix = wilson_b(pos_bohr, coords)
    g = b_matrix @ b_matrix.T
    a_matrix = _pseudo_inverse(g) @ b_matrix
    h_q = a_matrix @ h_au @ a_matrix.T  # 对角元：bond E_h/a₀²，angle E_h/rad²

    terms: list[DerivedTerm] = []
    for row, (i, j) in enumerate(bonds):
        terms.append(
            DerivedTerm(
                atoms=(i, j), value=bond_value(pos, i, j),
                k=float(h_q[row, row]) * _STRETCH_AU_TO_IR,
            )
        )
    offset = len(bonds)
    for row, (i, j, k_atom) in enumerate(angles):
        terms.append(
            DerivedTerm(
                atoms=(i, j, k_atom),
                value=angle_value(pos, i, j, k_atom),
                k=float(h_q[offset + row, offset + row]) * _ANGLE_AU_TO_IR,
            )
        )
    return terms
