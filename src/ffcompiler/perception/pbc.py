"""周期性最近镜像工具（Sobtop 的 getbondvec 语义，决策 0022 吸收）。

盒约定与 GROMACS 一致：``box`` 为 3×3 矩阵，**行**是晶胞向量 a/b/c [nm]；
正交盒可传对角向量的 3×3 对角阵。``box=None`` 表示孤立体系（不施加镜像）。
"""

from __future__ import annotations

import numpy as np

__all__ = ["nearest_image", "distance_vectors"]


def nearest_image(delta_nm: np.ndarray, box_nm: np.ndarray | None) -> np.ndarray:
    """把笛卡尔位置差折回最近镜像（最小镜像约定，支持三斜）。

    对三斜盒采用分数坐标 round 法，与 Sobtop 的最近镜像搜索等价
    （常规盒向量下两者给出同一镜像）。
    """
    if box_nm is None:
        return np.asarray(delta_nm, dtype=float)
    delta = np.asarray(delta_nm, dtype=float)
    box = np.asarray(box_nm, dtype=float)
    frac = np.linalg.solve(box.T, delta)
    frac -= np.round(frac)
    return box.T @ frac


def distance_vectors(
    positions_nm: np.ndarray, box_nm: np.ndarray | None
) -> np.ndarray:
    """全部原子对的最近镜像位置差 r_j − r_i，shape (N, N, 3)。

    Sobtop 的 genneighlist 用邻居表避免 O(N²)；这里面向中小体系直接展开
    （324 原子的 SiO₂ 超胞 52k 对，NumPy 下瞬时完成）。更大的体系交给调用方分块。
    """
    pos = np.asarray(positions_nm, dtype=float)
    delta = pos[None, :, :] - pos[:, None, :]
    if box_nm is None:
        return delta
    box = np.asarray(box_nm, dtype=float)
    frac = np.linalg.solve(box.T, delta.reshape(-1, 3).T).T.reshape(delta.shape)
    frac -= np.round(frac)
    return frac @ box
