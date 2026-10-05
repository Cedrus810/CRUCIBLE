"""mSeminario（Allen et al., J. Chem. Theory Comput. 2018, 14, 274；MVP-1 计划 T17）。

记号同计划 §2.6：K_AB = −H_AB（3×3 原子间块，不对称化），λ_i、v_i 是 K_AB
的本征值与本征向量（取实部）。投影力常数 k(A,B;u) = Σ_i λ_i |u·v_i|。

- bond (i, j)：k = k(i,j; u_ij)，u_ij 从 i 指向 j。**不对 k_ij、k_ji 取平均**，
  键方向即输入顺序（Sobtop mSeminario 的行为，决策 0011）。
- angle A-B-C：u_N = û(u_CB×u_AB)、u_PA = u_N×u_AB、u_PC = u_CB×u_N，
  1/kθ = f_A/(R_AB²·k(A,B;u_PA)) + f_C/(R_CB²·k(C,B;u_PC))，再对
  (A,B,C) 与 (C,B,A) 两种写法取平均。f 是相邻键角缩放因子（见
  `adjacent_angle_factor`）。

公式按原文结构实现，并以 Sobtop mSeminario 输出逐项核对（12 个 oracle
case 全部在 §2.8 容差内，决策 0011）；原文公式编号待原文到位后补进理论文档。
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from ffcompiler.bonded.internal import angle_value, linear_angles, LinearAngleError
from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.hessian.base import CartesianHessian

# 投影力常数：(h, a, b, u) → 沿 u 的 A-B 间力常数
Projector = Callable[[CartesianHessian, int, int, np.ndarray], float]


def _k_block(h: CartesianHessian, a: int, b: int) -> np.ndarray:
    return -h.hessian[3 * a:3 * a + 3, 3 * b:3 * b + 3]


def eigen_projection(h: CartesianHessian, a: int, b: int, u: np.ndarray) -> float:
    """k(A,B;u) = Σ_i λ_i |u·v_i|，λ、v 为 K_AB 的本征对（实部）。"""
    eigenvalues, eigenvectors = np.linalg.eig(_k_block(h, a, b))
    return float(np.sum(eigenvalues.real * np.abs(eigenvectors.real.T @ u)))


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.linalg.norm(x)


def _u_pa(pos: np.ndarray, a: int, b: int, c: int) -> np.ndarray:
    """A-B-C 平面内、垂直于 B→A 的单位向量 u_PA = u_N × u_AB。"""
    u_ab = _unit(pos[a] - pos[b])
    u_cb = _unit(pos[c] - pos[b])
    u_n = _unit(np.cross(u_cb, u_ab))
    return np.cross(u_n, u_ab)


def adjacent_angle_factor(
    pos: np.ndarray, angles: Sequence[tuple[int, int, int]], a: int, b: int, c: int
) -> float:
    """相邻键角缩放因子 f_A（Allen 2018）。

    以 B 为顶点、含键 A-B 的其余键角 A-B-D（D ≠ C）共 n−1 个：
    f_A = 1 + Σ_D |u_PA(A,B,C)·u_PA(A,B,D)|² / (n−1)；没有其余键角时 f_A = 1。
    """
    others = []
    for x, y, z in angles:
        if y != b:
            continue
        if x == a and z != c:
            others.append(z)
        elif z == a and x != c:
            others.append(x)
    if not others:
        return 1.0
    u_ref = _u_pa(pos, a, b, c)
    total = sum(abs(float(u_ref @ _u_pa(pos, a, b, d))) ** 2 for d in others)
    return 1.0 + total / len(others)


def scaled_angle_k(
    h: CartesianHessian,
    angles: Sequence[tuple[int, int, int]],
    a: int,
    b: int,
    c: int,
    projector: Projector,
) -> float:
    """1/kθ = f_A/(R_AB²·k_PA) + f_C/(R_CB²·k_PC)，对两种端点顺序取平均。"""
    pos = h.positions_nm

    def one_side(first: int, last: int) -> float:
        r1 = float(np.linalg.norm(pos[first] - pos[b]))
        r2 = float(np.linalg.norm(pos[last] - pos[b]))
        k1 = projector(h, first, b, _u_pa(pos, first, b, last))
        k2 = projector(h, last, b, _u_pa(pos, last, b, first))
        f1 = adjacent_angle_factor(pos, angles, first, b, last)
        f2 = adjacent_angle_factor(pos, angles, last, b, first)
        return 1.0 / (f1 / (r1**2 * k1) + f2 / (r2**2 * k2))

    return 0.5 * (one_side(a, c) + one_side(c, a))


def mseminario(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
) -> list[DerivedTerm]:
    """mSeminario 推导 bond/angle 参数（返回顺序：先全部 bond 再全部 angle）。

    k 为二阶导（½k 约定）；平衡值取 Hessian 参考几何上的值。
    输入含线性键角 → LinearAngleError。
    """
    linear = linear_angles(h.positions_nm, list(angles))
    if linear:
        raise LinearAngleError(linear)
    pos = h.positions_nm
    terms: list[DerivedTerm] = []
    for i, j in bonds:
        terms.append(
            DerivedTerm(
                atoms=(i, j),
                value=float(np.linalg.norm(pos[j] - pos[i])),
                k=eigen_projection(h, i, j, _unit(pos[j] - pos[i])),
            )
        )
    for i, j, k_atom in angles:
        terms.append(
            DerivedTerm(
                atoms=(i, j, k_atom),
                value=angle_value(pos, i, j, k_atom),
                k=scaled_angle_k(h, angles, i, j, k_atom, eigen_projection),
            )
        )
    return terms
