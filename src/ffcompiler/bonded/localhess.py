"""局部 Hessian 柔量法（Sobtop 隐藏 k_method=5 的复刻，决策 0023，experimental）。

来源：Sobtop 2026.1.16 逆向（docs/theory/sobtop_reverse_engineering.md §5.1、§8.2）。
``calcbondparm/calcangleparm/caldihparm`` 的 ``k_method==5`` 分支调用
``genlochess``（抽取项所涉原子的 3n×3n Cartesian Hessian 局部块），取该项的
Wilson B 行 b，对局部块做 SVD 伪逆（Sobtop 阈值 |σ|>1e-10 绝对值），最后按
``k = bᵀ·M·b`` 形式的二次型收尾（汇编双循环 0x43601e–0x436086 确证）。

**该菜单无入口（无法生成 oracle），汇编有两处无法裁决的歧义**：

1. 二次型里的 M 是 pinv(H_loc) 还是 pinv(pinv(H_loc))（谱截断后的 H_loc）；
2. 收尾是否漏了取倒数——bᵀ·pinv(H_loc)·b 的量纲是柔量 1/k，直接当 k 用
   量纲不自洽（这也可能是作者未放出该方法的原因之一）。

本项目采用**量纲自洽**的解释（局部柔量/compliance 常数，Decius 1963 系）：

    k_t = 1 / (bᵀ · pinv(H_loc, σ>1e-8·σ_max) · b)

即局部块伪逆的柔量沿 B 行投影后取倒数。刚体平移零空间由 SVD 截断处理。
键 (i,j) 用 6×6 块、角 (i,j,k) 用 9×9 块。扭转按 sobtop 语义同样适用，但
本项目柔性扭转走 MLP 管线（设计 §10.1），这里只做 bond/angle。

单位：输入 Hessian 为 kJ/mol/nm²（CartesianHessian 约定），B 行按 nm/rad
量纲构造，k 输出 kJ/mol/nm²（键）与 kJ/mol/rad²（角），½k 约定。
对 i,j 交换严格不变（b 与 H_loc 在交换下都不变），无 mSeminario 的方向
不对称问题；线性角（θ→180°，sinθ→0）在 B 行奇异，与 mseminario 一律
LinearAngleError。
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ffcompiler.bonded.internal import angle_value, linear_angles, LinearAngleError
from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.hessian.base import CartesianHessian

# 伪逆截断：σ < 1e-8·σ_max 置零（与 drih 的 G 阈值决策 0012 同精神；
# Sobtop 原始实现是 Eh/Bohr² 下的绝对 1e-10，单位换算后不等价，此处取相对阈值）
_SV_CUTOFF = 1e-8


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.linalg.norm(x)


def _local_block(h: CartesianHessian, atoms: Sequence[int]) -> np.ndarray:
    """抽 |atoms| 个原子的 3n×3n Cartesian 块（kJ/mol/nm²）。"""
    idx = [a for a in atoms for a in (3 * a, 3 * a + 1, 3 * a + 2)]
    return h.hessian[np.ix_(idx, idx)]


def _pinv_cutoff(matrix: np.ndarray) -> np.ndarray:
    u, s, vt = np.linalg.svd(matrix)
    cutoff = _SV_CUTOFF * (s[0] if s.size else 0.0)
    s_inv = np.where(s > cutoff, 1.0 / np.where(s > cutoff, s, 1.0), 0.0)
    return (u * s_inv) @ vt


def bond_b_row(pos: np.ndarray, i: int, j: int) -> np.ndarray:
    """键长的 Wilson B 行（对 i、j 的 6 个分量）。

    r = |x_j − x_i|：∂r/∂x_i = (x_i−x_j)/r = −û_ij，∂r/∂x_j = û_ij。
    """
    u = _unit(pos[j] - pos[i])
    return np.concatenate([-u, u])


def angle_b_row(pos: np.ndarray, i: int, j: int, k: int) -> np.ndarray:
    """键角的 Wilson B 行（对 i、j、k 的 9 个分量，nm⁻¹ 量纲）。

    cosθ = u_ji·u_jk；∂θ/∂x_i = (cosθ·u_ji − u_jk)/(sinθ·|r_ji|)，
    ∂θ/∂x_k 对称，∂θ/∂x_j = −(∂θ/∂x_i + ∂θ/∂x_k)（平移不变）。
    """
    u_ji = _unit(pos[i] - pos[j])
    u_jk = _unit(pos[k] - pos[j])
    r1 = float(np.linalg.norm(pos[i] - pos[j]))
    r2 = float(np.linalg.norm(pos[k] - pos[j]))
    cos_t = float(u_ji @ u_jk)
    sin_t = float(np.sqrt(max(1.0 - cos_t * cos_t, 0.0)))
    if sin_t < 1e-8:
        raise LinearAngleError([(i, j, k)])
    bi = (cos_t * u_ji - u_jk) / (sin_t * r1)
    bk = (cos_t * u_jk - u_ji) / (sin_t * r2)
    return np.concatenate([bi, -(bi + bk), bk])


def local_compliance_k(
    h: CartesianHessian, atoms: Sequence[int], b_row: np.ndarray
) -> float:
    """k = 1/(bᵀ·pinv(H_loc)·b)；柔量非正时返回负值（由上层报方法局限）。"""
    block = _local_block(h, atoms)
    compliance = float(b_row @ _pinv_cutoff(block) @ b_row)
    if compliance == 0.0 or not np.isfinite(compliance):
        raise NonFiniteComplianceError(atoms, compliance)
    return 1.0 / compliance


class NonFiniteComplianceError(ValueError):
    """局部柔量为零/非有限（局部块全零或 B 行落入零空间）。"""

    def __init__(self, atoms, value: float):
        self.atoms = tuple(atoms)
        self.value = value
        super().__init__(
            f"localhess：项 {self.atoms} 的局部柔量为 {value!r}（局部块退化或 B 行"
            "落入截断零空间），无法取倒数。"
        )


def localhess(
    h: CartesianHessian,
    bonds: Sequence[tuple[int, int]],
    angles: Sequence[tuple[int, int, int]],
) -> list[DerivedTerm]:
    """局部柔量法推导 bond/angle（返回顺序：先全部 bond 再全部 angle）。

    k 为二阶导（½k 约定）；平衡值取 Hessian 参考几何上的值。
    输入含线性键角 → LinearAngleError（与 mseminario 相同行为）。
    """
    linear = linear_angles(h.positions_nm, list(angles))
    if linear:
        raise LinearAngleError(linear)
    pos = h.positions_nm
    terms: list[DerivedTerm] = []
    nonpositive: list[tuple[tuple[int, ...], float]] = []
    for i, j in bonds:
        k = local_compliance_k(h, (i, j), bond_b_row(pos, i, j))
        if k <= 0.0:
            nonpositive.append(((i, j), k))
        terms.append(
            DerivedTerm(
                atoms=(i, j),
                value=float(np.linalg.norm(pos[j] - pos[i])),
                k=k,
            )
        )
    for i, j, k_atom in angles:
        k = local_compliance_k(h, (i, j, k_atom), angle_b_row(pos, i, j, k_atom))
        if k <= 0.0:
            nonpositive.append(((i, j, k_atom), k))
        terms.append(
            DerivedTerm(
                atoms=(i, j, k_atom),
                value=angle_value(pos, i, j, k_atom),
                k=k,
            )
        )
    if nonpositive:
        from ffcompiler.bonded.assign import NonPositiveForceConstantError

        raise NonPositiveForceConstantError("localhess", nonpositive)
    return terms
