"""几何成键感知：键 / 角 / 二面角 / improper 中心（决策 0022，吸收自 Sobtop 逆向）。

规则与 Sobtop 逐条对齐（动态实验 ground truth 见
docs/theory/sobtop_reverse_engineering.md §4、§7.3 与附录 A）：

- 键：``d < (r_i + r_j) × factor``（Cordero 2008 半径，factor 默认 1.15），
  bondcrit.dat 的元素对**距离窗口**覆盖优先于半径判据；周期体系用最近镜像距离；
- 角：每个中心原子的邻居两两组合（含氢、不去重、不剔除线性角）；
- 二面角：所有 i-j-k-l 路径；任一侧角满足 ``180° − θ < linear_deg``（默认 1°）
  即剔除；输出规范化为 i<l 方向（Sobtop 的 (i,j,k,l)/(l,k,j,i) 等价形式）；
- improper 中心：配位数恰为 3 且原子类型 ∈ 库的 $IMPROPER 第三列集合，
  邻居升序 (j₁, j₂, center, j₃)。
"""

from __future__ import annotations

import numpy as np

from .pbc import distance_vectors, nearest_image
from .radii import covalent_radius

__all__ = [
    "perceive_bonds",
    "build_angles",
    "build_dihedrals",
    "find_improper_centers",
    "linear_angle_deg",
    "DEFAULT_BOND_FACTOR",
    "LINEAR_CUTOFF_DEG",
    "ELEMENT_Z",
]

DEFAULT_BOND_FACTOR = 1.15
LINEAR_CUTOFF_DEG = 1.0

# 元素符号（大写归一）→ 原子序数；bondcrit.dat 的覆盖规则用元素符号书写
ELEMENT_Z = {
    "H": 1, "D": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15, "S": 16,
    "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22, "V": 23,
    "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29, "ZN": 30,
    "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36, "RB": 37,
    "SR": 38, "Y": 39, "ZR": 40, "NB": 41, "MO": 42, "TC": 43, "RU": 44,
    "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50, "SB": 51,
    "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57, "CE": 58,
    "PR": 59, "ND": 60, "PM": 61, "SM": 62, "EU": 63, "GD": 64, "TB": 65,
    "DY": 66, "HO": 67, "ER": 68, "TM": 69, "YB": 70, "LU": 71, "HF": 72,
    "TA": 73, "W": 74, "RE": 75, "OS": 76, "IR": 77, "PT": 78, "AU": 79,
    "HG": 80, "TL": 81, "PB": 82, "BI": 83, "PO": 84, "AT": 85, "RN": 86,
    "FR": 87, "RA": 88, "AC": 89, "TH": 90, "PA": 91, "U": 92,
}


def _neighbor_lists(n_atoms: int, bonds: list[tuple[int, int]]) -> list[list[int]]:
    nbr: list[list[int]] = [[] for _ in range(n_atoms)]
    for i, j in bonds:
        nbr[i].append(j)
        nbr[j].append(i)
    return [sorted(x) for x in nbr]


def linear_angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    """两键向量的夹角 [°]。"""
    cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _angle_deg(
    pos: np.ndarray, i: int, j: int, k: int, box_nm: np.ndarray | None
) -> float:
    u = nearest_image(pos[i] - pos[j], box_nm)
    v = nearest_image(pos[k] - pos[j], box_nm)
    return linear_angle_deg(u, v)


def perceive_bonds(
    atomic_numbers: np.ndarray | list[int],
    positions_nm: np.ndarray,
    box_nm: np.ndarray | None = None,
    factor: float = DEFAULT_BOND_FACTOR,
    pair_windows: dict[tuple[int, int], tuple[float, float]] | None = None,
) -> list[tuple[int, int]]:
    """按半径判据感知键（0 基索引，i<j，(i,j) 字典序输出）。

    ``pair_windows``：元素对 (Z_i, Z_j) 无序键 → (下限, 上限) [Å] 距离窗口
    （bondcrit.dat）；落在窗口内（含边界）即成键，优先于半径判据。
    """
    z = np.asarray(atomic_numbers, dtype=int)
    n = len(z)
    if n < 2:
        return []
    d_nm = np.linalg.norm(distance_vectors(positions_nm, box_nm), axis=-1)
    radii_a = np.array([covalent_radius(int(zi)) for zi in z])
    crit_nm = (radii_a[:, None] + radii_a[None, :]) * factor * 0.1
    bonded = d_nm < crit_nm
    bonded[np.diag_indices(n)] = False
    for (za, zb), (lo_a, hi_a) in (pair_windows or {}).items():
        mask = ((z[:, None] == za) & (z[None, :] == zb)) | (
            (z[:, None] == zb) & (z[None, :] == za)
        )
        bonded |= mask & (d_nm >= lo_a * 0.1) & (d_nm <= hi_a * 0.1)
    idx = np.argwhere(np.triu(bonded, 1))
    return [(int(i), int(j)) for i, j in idx]


def build_angles(
    bonds: list[tuple[int, int]], n_atoms: int
) -> list[tuple[int, int, int]]:
    """每个中心原子邻居的两两组合 (i, j, k)，按中心原子升序输出（Sobtop 顺序）。"""
    out: list[tuple[int, int, int]] = []
    for j, nbrs in enumerate(_neighbor_lists(n_atoms, bonds)):
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                out.append((nbrs[a], j, nbrs[b]))
    return out


def build_dihedrals(
    bonds: list[tuple[int, int]],
    n_atoms: int,
    positions_nm: np.ndarray,
    box_nm: np.ndarray | None = None,
    linear_deg: float = LINEAR_CUTOFF_DEG,
) -> list[tuple[int, int, int, int]]:
    """全部 i-j-k-l 路径；任一侧角 180°−θ < linear_deg 即剔除；规范 i<l 输出。"""
    pos = np.asarray(positions_nm, dtype=float)
    out: list[tuple[int, int, int, int]] = []
    nbrs = _neighbor_lists(n_atoms, bonds)
    for j in range(n_atoms):
        for k in nbrs[j]:
            if k <= j:
                continue  # 每条中心键只处理一次（j<k）
            for i in nbrs[j]:
                if i == k:
                    continue
                if _angle_deg(pos, i, j, k, box_nm) <= 180.0 - linear_deg:
                    for l in nbrs[k]:
                        if l == j or l == i:
                            continue
                        if _angle_deg(pos, j, k, l, box_nm) <= 180.0 - linear_deg:
                            out.append((i, j, k, l) if i < l else (l, k, j, i))
    out.sort()
    return out


def find_improper_centers(
    atom_types: list[str],
    bonds: list[tuple[int, int]],
    center_types: set[str],
) -> list[tuple[int, int, int, int]]:
    """数据驱动 improper 感知（Sobtop outitprtp 规则）。

    配位数恰为 3 且类型 ∈ center_types（fflib.improper_center_types 从
    $IMPROPER 第三列导出）→ 一条 improper，邻居升序 (j₁, j₂, center, j₃)。
    """
    n = len(atom_types)
    nbrs = _neighbor_lists(n, bonds)
    out: list[tuple[int, int, int, int]] = []
    for center in range(n):
        if len(nbrs[center]) != 3:
            continue
        if atom_types[center].strip() not in center_types:
            continue
        j1, j2, j3 = nbrs[center]
        out.append((j1, j2, center, j3))
    return out
