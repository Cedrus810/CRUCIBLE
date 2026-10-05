"""m2Seminario（MVP-1 计划 T18、§2.6；命名见决策 0011）。

两个候选（§2.6 锁定）都在本模块实现：
- C1（直接投影，`_c1`）：bond k_AB = u_ABᵀ K_AB u_AB；angle
  1/kθ = 1/(R_AB²·u_PAᵀK_AB·u_PA) + 1/(R_CB²·u_PCᵀK_CB·u_PC)。
- C2（`m2seminario`）：bond 同 C1；angle 在 C1 上乘 mSeminario 的相邻键角
  缩放因子 f（与 `mseminario.scaled_angle_k` 同一实现，只把投影换成 uᵀKu）。

oracle 比对（决策 0011）：C2 在 12 个 oracle case 上与 Sobtop m2Seminario
全部一致；C1 的 bond 全部一致、angle 恰差因子 f。按 §2.6 命名规则，
公开方法名为 `m2seminario`。uᵀKu 是标量且 K_BA = K_ABᵀ ⇒ k_AB = k_BA。
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ffcompiler.bonded.internal import angle_value, linear_angles, LinearAngleError
from ffcompiler.bonded.mseminario import scaled_angle_k
from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.hessian.base import CartesianHessian


def quadratic_projection(h: CartesianHessian, a: int, b: int, u: np.ndarray) -> float:
    """uᵀ K_AB u，K_AB = −H_AB。"""
    k_block = -h.hessian[3 * a:3 * a + 3, 3 * b:3 * b + 3]
    return float(u @ k_block @ u)


def _c1_angle(h: CartesianHessian, i: int, j: int, k_atom: int) -> float:
    """C1：不带缩放因子的调和式（仅供候选比对，见决策 0011）。"""
    return scaled_angle_k(h, [(i, j, k_atom)], i, j, k_atom, quadratic_projection)


def _derive(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
    *,
    candidate: str,
) -> list[DerivedTerm]:
    linear = linear_angles(h.positions_nm, list(angles))
    if linear:
        raise LinearAngleError(linear)
    pos = h.positions_nm
    terms: list[DerivedTerm] = []
    for i, j in bonds:
        u = (pos[j] - pos[i]) / np.linalg.norm(pos[j] - pos[i])
        terms.append(
            DerivedTerm(
                atoms=(i, j),
                value=float(np.linalg.norm(pos[j] - pos[i])),
                k=quadratic_projection(h, i, j, u),
            )
        )
    for i, j, k_atom in angles:
        if candidate == "C1":
            k = _c1_angle(h, i, j, k_atom)
        else:
            k = scaled_angle_k(h, angles, i, j, k_atom, quadratic_projection)
        terms.append(
            DerivedTerm(atoms=(i, j, k_atom), value=angle_value(pos, i, j, k_atom), k=k)
        )
    return terms


def _c1(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
) -> list[DerivedTerm]:
    return _derive(h, bonds, angles, candidate="C1")


def m2seminario(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
) -> list[DerivedTerm]:
    """C2 推导 bond/angle 参数（先全部 bond 再全部 angle，½k 约定）。"""
    return _derive(h, bonds, angles, candidate="C2")
