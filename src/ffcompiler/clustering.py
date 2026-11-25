"""几何聚类的基本件（2026-11-20 计划 §B/§C1 共用，§10.4 分支描述符"同源一套实现"）。

- :func:`kabsch_align`：Kabsch 最优刚体对齐（只旋转+平移，不含反射）；
- :func:`heavy_rmsd_nm`：对齐后重原子 RMSD（nm）；
- :func:`greedy_cluster`：能量升序贪心聚类（§10.4 的逐格点分支标注、§C1 的留出盆地
  划分都用它）。平局按输入序号升序，结果确定；簇代表 = 簇内最低能量成员（首个加入）。

聚类的"能量"由调用方解释：§10.4 是候选扫描能量，§C1 是逐帧参考能（平移后）。
阈值单位 nm（IR 长度单位）；不设自动判废阈值，聚类只产生诊断/划分标签。
"""

from __future__ import annotations

import numpy as np


def kabsch_align(src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """把 ``src`` (N,3) 刚体对齐到 ``dst`` (N,3)：先各自去质心，再取最优旋转（SVD，
    行列式为负时翻转最后一个奇异方向以排除反射）。返回对齐后的 ``src`` 坐标 (N,3)。"""
    src = np.asarray(src, dtype=np.float64)
    dst = np.asarray(dst, dtype=np.float64)
    src_c = src - src.mean(axis=0)
    dst_c = dst - dst.mean(axis=0)
    u, _, vt = np.linalg.svd(src_c.T @ dst_c)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return (rot @ src_c.T).T + dst.mean(axis=0)


def heavy_rmsd_nm(a: np.ndarray, b: np.ndarray, heavy_mask: np.ndarray) -> float:
    """对齐后重原子 RMSD（nm）；``heavy_mask`` 是长度 N 的布尔数组。"""
    a = np.asarray(a, dtype=np.float64)[heavy_mask]
    b = np.asarray(b, dtype=np.float64)[heavy_mask]
    aligned = kabsch_align(a, b)
    return float(np.sqrt(np.mean(np.sum((aligned - b) ** 2, axis=1))))


def greedy_cluster(
    positions: np.ndarray,
    energies: np.ndarray,
    heavy_mask: np.ndarray,
    threshold_nm: float,
) -> np.ndarray:
    """能量升序贪心聚类。

    ``positions`` (M,N,3) nm、``energies`` (M,)；每个候选归入**第一个**代表与其
    重原子对齐 RMSD ≤ threshold 的簇，没有就自成一簇。返回长度 M 的簇号数组
    （0 起编，簇号按代表加入顺序 = 能量升序，所以簇 0 是最低能量的簇）。
    """
    positions = np.asarray(positions, dtype=np.float64)
    energies = np.asarray(energies, dtype=np.float64)
    m = len(energies)
    order = np.lexsort((np.arange(m), energies))
    reps: list[np.ndarray] = []
    ids = np.empty(m, dtype=int)
    for idx in order:
        assigned = None
        for c, rep in enumerate(reps):
            if heavy_rmsd_nm(positions[idx], rep, heavy_mask) <= threshold_nm:
                assigned = c
                break
        if assigned is None:
            reps.append(positions[idx])
            assigned = len(reps) - 1
        ids[idx] = assigned
    return ids
