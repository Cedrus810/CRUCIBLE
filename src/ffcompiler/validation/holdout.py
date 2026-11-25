"""留出集的盆地划分（2026-11-20 计划 §C1 = HANDOVER 待办 B4 的落地形式）。

目的：把扫描帧按构象盆地整进整出地分成 fit / holdout，避免相邻扫描帧互相泄漏
（同一盆地内的帧高度相关，逐帧随机划分会让留出集变成插值练习）。

- 描述符与 §10.4 分支标注同源：对齐后重原子 RMSD 聚类（共用 ``ffcompiler.clustering``
  的 :func:`greedy_cluster`），阈值沿用分支阈值 0.05 nm；
- 盆地 = 全部扫描帧的聚类簇；划分规则（施工前锁定）：
  1) 按簇大小升序（平局按簇内最小能量帧的下标升序）依次尝试把整个簇划入 holdout；
  2) 约束 holdout 帧数 ≤ ``max_holdout_frac`` × 总帧数；放不下的簇留在 fit；
  3) 只有 1 个盆地（或所有簇都放不下）时 holdout 为空并写 ``note``；
- 不做拟合、不算指标：指标聚合（相对能量/力/结构 RMSD/频率）属验证报告层。

确定性：聚类与划分都是确定性的（无随机成分）。
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..clustering import greedy_cluster
from ..torsion.scan import ScanResult

DEFAULT_THRESHOLD_NM = 0.05
DEFAULT_MAX_HOLDOUT_FRAC = 0.4


def basin_split(
    scan: ScanResult,
    heavy_mask: np.ndarray,
    *,
    threshold_nm: float = DEFAULT_THRESHOLD_NM,
    max_holdout_frac: float = DEFAULT_MAX_HOLDOUT_FRAC,
) -> dict:
    """把一条扫描的帧分成 fit / holdout（整盆地进整盆地出）。

    返回 ``{"fit_indices", "holdout_indices", "n_basins", "basin_of_frame",
    "threshold_nm", "max_holdout_frac", "note"}``；索引是扫描帧下标。
    """
    n = len(scan.energy)
    if scan.positions.shape[0] != n:
        raise ValueError(f"energy ({n}) 与 positions ({scan.positions.shape[0]}) 帧数不一致")
    basin = greedy_cluster(scan.positions, scan.energy, heavy_mask, threshold_nm)
    n_basins = int(basin.max()) + 1
    members: list[list[int]] = [[] for _ in range(n_basins)]
    for s, b in enumerate(basin):
        members[int(b)].append(s)

    cap = int(np.floor(max_holdout_frac * n))
    order = sorted(range(n_basins), key=lambda b: (len(members[b]), min(members[b])))
    holdout: list[int] = []
    for b in order:
        if len(holdout) + len(members[b]) <= cap:
            holdout.extend(members[b])
    note = None
    if not holdout:
        note = (f"holdout 为空：{n_basins} 个盆地、cap={cap} 帧；"
                "单盆地扫描或最小盆地也超过上限时无留出（放宽 max_holdout_frac 或增密扫描）")
    fit = sorted(set(range(n)) - set(holdout))
    return {
        "fit_indices": fit,
        "holdout_indices": sorted(holdout),
        "n_basins": n_basins,
        "basin_of_frame": basin.tolist(),
        "threshold_nm": float(threshold_nm),
        "max_holdout_frac": float(max_holdout_frac),
        "note": note,
    }


def split_summary(split: Mapping) -> str:
    """单行摘要（报告用）：盆地数、fit/holdout 帧数。"""
    return (f"basins={split['n_basins']} fit={len(split['fit_indices'])} "
            f"holdout={len(split['holdout_indices'])}"
            + (f" [{split['note']}]" if split.get("note") else ""))
