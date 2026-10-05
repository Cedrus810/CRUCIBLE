"""内坐标、Wilson B 矩阵与线性键角判定（MVP-1 计划 T14）。

坐标类型：bond (i,j)、angle (i,j,k)（j 为顶点）、dihedral (i,j,k,l)。
二面角定义与 MVP-0 §2.1 一致（顺式为 0），直接复用 ``ir.energy.dihedral_angle``。
B 的行是坐标对 3N 个 Cartesian 分量的一阶导数：键无量纲，角为 rad/nm。
"""

from __future__ import annotations

import math
from typing import Literal, Mapping, Sequence

import numpy as np

from ffcompiler.ir.energy import dihedral_angle

InternalCoord = tuple[Literal["bond", "angle", "dihedral"], tuple[int, ...]]
LINEAR_ANGLE_DEG = 175.0


class LinearAngleError(ValueError):
    """内坐标集合含线性键角（计划 §2.4）。

    angles 列出全部线性键角 (i, j, k)，j 为顶点。
    """

    def __init__(self, angles: Sequence[tuple[int, int, int]], message: str | None = None):
        self.angles = list(angles)
        super().__init__(
            message
            or f"内坐标集合含 {len(self.angles)} 个线性键角"
            f"（> {LINEAR_ANGLE_DEG}°）：{self.angles}"
        )


def bond_value(positions_nm: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(positions_nm[j] - positions_nm[i]))


def angle_value(positions_nm: np.ndarray, i: int, j: int, k: int) -> float:
    """顶点在 j 的键角，atan2 形式，范围 [0, π]，线性键角处 sin→0 不奇异。"""
    u1 = positions_nm[i] - positions_nm[j]
    u2 = positions_nm[k] - positions_nm[j]
    return math.atan2(
        float(np.linalg.norm(np.cross(u1, u2))), float(np.dot(u1, u2))
    )


def coord_value(positions_nm: np.ndarray, coord: InternalCoord) -> float:
    """内坐标当前值：键为 nm，角和二面角为 rad。"""
    kind, idx = coord
    pos = np.asarray(positions_nm, dtype=np.float64)
    if kind == "bond":
        i, j = idx
        return bond_value(pos, i, j)
    if kind == "angle":
        i, j, k = idx
        return angle_value(pos, i, j, k)
    if kind == "dihedral":
        i, j, k, l = idx
        return dihedral_angle(pos[i], pos[j], pos[k], pos[l])
    raise ValueError(f"未知内坐标类型 {kind!r}")


def linear_angles(
    positions_nm: np.ndarray,
    angles: Sequence[tuple[int, int, int]],
) -> list[tuple[int, int, int]]:
    """返回角度 > LINEAR_ANGLE_DEG 的键角（计划 §2.4）。"""
    pos = np.asarray(positions_nm, dtype=np.float64)
    return [
        (i, j, k)
        for (i, j, k) in angles
        if math.degrees(angle_value(pos, i, j, k)) > LINEAR_ANGLE_DEG
    ]


def _bond_row(pos: np.ndarray, i: int, j: int, n3: int) -> np.ndarray:
    """∂|r_j − r_i|/∂r：Wilson 键坐标。"""
    vec = pos[j] - pos[i]
    r = float(np.linalg.norm(vec))
    if r == 0.0:
        raise ValueError(f"原子 {i} 和 {j} 重合，键坐标奇异")
    unit = vec / r
    row = np.zeros(n3)
    row[3 * i:3 * i + 3] = -unit
    row[3 * j:3 * j + 3] = unit
    return row


def _angle_row(pos: np.ndarray, i: int, j: int, k: int, n3: int) -> np.ndarray:
    """∂θ/∂r，顶点在 j。解析公式（Wilson）：

    ∂θ/∂r_i = (cosθ·â − b̂)/(R_a sinθ)，∂θ/∂r_k = (cosθ·b̂ − â)/(R_b sinθ)，
    ∂θ/∂r_j = −(∂θ/∂r_i + ∂θ/∂r_k)（平移不变）。
    """
    u1 = pos[i] - pos[j]
    u2 = pos[k] - pos[j]
    ra = float(np.linalg.norm(u1))
    rb = float(np.linalg.norm(u2))
    cross_uv = np.cross(u1, u2)
    sin_theta = float(np.linalg.norm(cross_uv)) / (ra * rb)
    cos_theta = float(np.dot(u1, u2)) / (ra * rb)
    if sin_theta < 1e-12:
        # 线性键角（θ→0 或 π）下公式奇异；正常由 wilson_b 的检查拦截
        raise LinearAngleError([(i, j, k)])
    a_hat = u1 / ra
    b_hat = u2 / rb
    s_i = (cos_theta * a_hat - b_hat) / (ra * sin_theta)
    s_k = (cos_theta * b_hat - a_hat) / (rb * sin_theta)
    s_j = -(s_i + s_k)
    row = np.zeros(n3)
    row[3 * i:3 * i + 3] = s_i
    row[3 * j:3 * j + 3] = s_j
    row[3 * k:3 * k + 3] = s_k
    return row


def _dihedral_row(
    pos: np.ndarray, i: int, j: int, k: int, l: int, n3: int
) -> np.ndarray:
    """∂φ/∂r，Blondel–Karplus 型解析公式（链式法则展开 atan2）。

    b1 = r_j − r_i, b2 = r_k − r_j, b3 = r_l − r_k；n1 = b1×b2, n2 = b2×b3；
    x = n1·n2, y = |b2|·det[b1,b2,b3]，φ = atan2(y, x)，x²+y² = |n1|²|n2|²。
    对 b1/b2/b3 的梯度为
      g_b1 = (x·|b2|·n2 + y·(b2×n2))/(x²+y²)
      g_b3 = (x·|b2|·n1 + y·(b2×n1))/(x²+y²)
      ∇_{b2}y = D·b2/|b2| + |b2|·(b3×b1)，∇_{b2}x = (b2·b3)b1 + (b1·b2)b3 − 2(b1·b3)b2
      g_b2 = (x·∇_{b2}y − y·∇_{b2}x)/(x²+y²)
    原子行：g_i = −g_b1，g_j = g_b1 − g_b2，g_k = g_b2 − g_b3，g_l = g_b3
    （自动满足平移不变）。
    """
    b1 = pos[j] - pos[i]
    b2 = pos[k] - pos[j]
    b3 = pos[l] - pos[k]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    norm_b2 = float(np.linalg.norm(b2))
    x = float(np.dot(n1, n2))
    d = float(np.dot(b1, n2))  # det[b1,b2,b3]
    y = norm_b2 * d
    denom = x * x + y * y
    if denom == 0.0:
        raise LinearAngleError([(i, j, k), (j, k, l)])

    # ∂x/∂b1 = (b2·b3)b2 − |b2|²b3 = +b2×n2；∂y/∂b1 = |b2|·n2
    g_b1 = (x * norm_b2 * n2 - y * np.cross(b2, n2)) / denom
    # ∂x/∂b3 = (b1·b2)b2 − |b2|²b1 = −b2×n1；∂y/∂b3 = |b2|·n1
    g_b3 = (x * norm_b2 * n1 + y * np.cross(b2, n1)) / denom
    grad_y_b2 = d * b2 / norm_b2 + norm_b2 * np.cross(b3, b1)
    grad_x_b2 = (
        float(np.dot(b2, b3)) * b1
        + float(np.dot(b1, b2)) * b3
        - 2.0 * float(np.dot(b1, b3)) * b2
    )
    g_b2 = (x * grad_y_b2 - y * grad_x_b2) / denom

    row = np.zeros(n3)
    row[3 * i:3 * i + 3] = -g_b1
    row[3 * j:3 * j + 3] = g_b1 - g_b2
    row[3 * k:3 * k + 3] = g_b2 - g_b3
    row[3 * l:3 * l + 3] = g_b3
    return row


def wilson_b(
    positions_nm: np.ndarray, coords: Sequence[InternalCoord]
) -> np.ndarray:
    """内坐标集合的 Wilson B 矩阵，形状 (n_q, 3N)。

    遇到线性键角（含二面角经过的键角）抛 LinearAngleError，列出全部。
    """
    pos = np.asarray(positions_nm, dtype=np.float64)
    n3 = 3 * pos.shape[0]

    all_linear: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()
    for kind, idx in coords:
        contained: list[tuple[int, int, int]] = []
        if kind == "angle":
            contained.append((idx[0], idx[1], idx[2]))
        elif kind == "dihedral":
            contained.append((idx[0], idx[1], idx[2]))
            contained.append((idx[1], idx[2], idx[3]))
        for triple in contained:
            if triple not in seen:
                seen.add(triple)
                if angle_value(pos, *triple) > math.radians(LINEAR_ANGLE_DEG):
                    all_linear.append(triple)
    if all_linear:
        raise LinearAngleError(all_linear)

    rows = np.zeros((len(coords), n3))
    for row_idx, (kind, idx) in enumerate(coords):
        if kind == "bond":
            rows[row_idx] = _bond_row(pos, idx[0], idx[1], n3)
        elif kind == "angle":
            rows[row_idx] = _angle_row(pos, idx[0], idx[1], idx[2], n3)
        elif kind == "dihedral":
            rows[row_idx] = _dihedral_row(pos, idx[0], idx[1], idx[2], idx[3], n3)
        else:
            raise ValueError(f"未知内坐标类型 {kind!r}")
    return rows


def proper_quartets(
    bond_graph: Mapping[int, frozenset[int]],
) -> list[tuple[int, int, int, int]]:
    """键图生成的全部 proper 二面角四元组（计划 §2.5 第 1 条）。

    (i,j,k,l) 与 (l,k,j,i) 是同一个二面角，只留 i<l 的那个；
    中心键每条只处理一次，输出按字典序排序。
    """
    quartets: set[tuple[int, int, int, int]] = set()
    for j, nbrs_j in bond_graph.items():
        for k in nbrs_j:
            if j >= k:
                continue
            for i in bond_graph[j]:
                if i == k:
                    continue
                for l in bond_graph[k]:
                    if l == j or l == i:
                        continue
                    quartets.add((i, j, k, l) if i < l else (l, k, j, i))
    return sorted(quartets)
