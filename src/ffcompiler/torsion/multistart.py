"""多起点扫描的起点生成与逐格点合并（设计 §10.4，MVP-7；2026-11-20 计划 §B）。

职责边界：
- 起点集合 = {x0（输入几何，现行协议）} ∪ {ETKDG 构象经参考势优化后的极小}。
  本模块只生成 ETKDG 起点（嵌入 + 参考势优化 + 去重）；x0 腿就是现行 relaxed_scan。
- 每个起点一次 relaxed_scan（正反两遍，机制不动）；非 x0 起点用
  ``grid_origin_rad`` 对齐 x0 的 φ0，使逐格点合并良定义。
- :func:`merge_start_scans` 把全部起点 × 方向的候选逐格点聚类（对齐重原子 RMSD，
  共用 ``ffcompiler.clustering``），产出包络 :class:`ScanResult`（Fourier 拟合消费，
  与现行"取两遍较低者"语义同构）与诊断字典（逐点分支数、分支间最大能量差、
  包络相对现行单起点协议的能量变化）。

裁定（沿设计稿 §10.4 原文）：分支只进诊断、永不进拟合目标；不设自动兜底阈值；
默认扫描协议不变（是否转默认由重扫证据决定）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .scan import ScanResult
from ..clustering import greedy_cluster, heavy_rmsd_nm
from ..units import ANGSTROM_NM

# 施工裁定（2026-11-20 计划 §B；待拍板项 4 未覆盖此项，可由 spec.threshold_nm 覆盖）：
# 分支聚类阈值 = 对齐后重原子 RMSD，单位 nm。0.05 nm = 0.5 Å，与 ETKDG 去重同量级。
BRANCH_THRESHOLD_NM = 0.05
DEFAULT_MULTISTART = {
    "n_confs": 3,        # ETKDG 嵌入的构象数（含被去重跳过的）
    "seed": 202611,      # ETKDG 随机种子（固定，保证起点可复现）
    "dedup_nm": 0.05,    # 起点间/对 x0 的去重阈值（优化后重原子 RMSD）
    "threshold_nm": BRANCH_THRESHOLD_NM,
    "opt_fmax": 0.005,   # ETKDG 起点的无约束优化判据（eV/Å，与扫描同）
    "opt_max_steps": 500,
}


@dataclass(frozen=True)
class StartSpec:
    """一个非 x0 起点：优化后几何 + 缓存身份（逐项字面值，进 scan_cache_key 的 start 块）。"""

    label: str
    positions_nm: np.ndarray
    identity: dict


def etkdg_starts(
    ir,
    mol,
    positions_nm: np.ndarray,
    calc,
    *,
    n_confs: int = 3,
    seed: int = 202611,
    dedup_nm: float = 0.05,
    opt_fmax: float = 0.005,
    opt_max_steps: int = 500,
) -> list[StartSpec]:
    """ETKDG 嵌入 → 参考势无约束优化 → 对齐重原子 RMSD 去重（对 x0 与已接受起点）。

    嵌入失败（部分分子 ETKDG 出空集）返回少于 n_confs 个甚至空列表——起点数量本身
    是诊断信息的一部分。身份含种子/序号/参数/RDKit 版本：同版本同种子逐位可复现，
    缓存键因此完整（"逐项字面值"原则，不用 hash）。优化判据与扫描一致（opt_fmax）。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from ase.optimize import LBFGS

    work = Chem.Mol(mol)  # 不改调用方的分子（调用方还要用它做拓扑检查）
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    conf_ids = list(AllChem.EmbedMultipleConfs(work, numConfs=int(n_confs), params=params))
    heavy = np.array([a.atomic_number > 1 for a in ir.atoms])
    refs: list[np.ndarray] = [np.asarray(positions_nm, dtype=np.float64)]
    out: list[StartSpec] = []
    for cid in conf_ids:
        pos_nm = np.asarray(work.GetConformer(cid).GetPositions(), dtype=np.float64) / ANGSTROM_NM
        atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos_nm * ANGSTROM_NM)
        atoms.calc = calc
        if not LBFGS(atoms, logfile=None).run(fmax=float(opt_fmax), steps=int(opt_max_steps)):
            continue  # 优化不收敛的嵌入不作为起点（扫描层对收敛有同样的门）
        opt_nm = atoms.get_positions() / ANGSTROM_NM
        if any(heavy_rmsd_nm(opt_nm, ref, heavy) < float(dedup_nm) for ref in refs):
            continue
        refs.append(opt_nm)
        out.append(StartSpec(
            label=f"etkdg-{int(cid)}",
            positions_nm=opt_nm,
            identity={"kind": "etkdg", "label": f"etkdg-{int(cid)}", "index": int(cid),
                      "seed": int(seed), "n_confs": int(n_confs), "dedup_nm": float(dedup_nm),
                      "etkdg": "ETKDGv3+smallRing", "opt_fmax": float(opt_fmax),
                      "opt_max_steps": int(opt_max_steps),
                      "rdkit_version": str(Chem.rdBase.rdkitVersion)},
        ))
    return out


def merge_start_scans(
    scans: Sequence[ScanResult],
    bonds: Sequence[tuple[int, int]],
    heavy_mask: np.ndarray,
    *,
    threshold_nm: float = BRANCH_THRESHOLD_NM,
) -> tuple[ScanResult, dict]:
    """把多个起点的扫描逐格点合并成包络 + 分支诊断。

    ``scans[0]`` 必须是 x0 腿（现行单起点协议的对照）。每个格点上候选 = 各起点的
    正反两方向（旧缓存没有候选块时退回该起点的逐点最优，即现行两遍取较低的结果）。
    收敛候选优先（与 relaxed_scan 的取法一致）；包络 = 全候选最低能量，平局按
    (起点序, 方向) 序取——确定性。包络的 ``hysteresis`` = 胜出候选同起点另一方向的
    能量差（胜出路径自己的正反差）。

    返回 ``(envelope, diagnostics)``；diagnostics 含逐点 branch_counts / branch_gap /
    envelope_delta（对 x0 的能量变化，≤0）与其汇总，供诊断报告与 provenance。
    """
    if not scans:
        raise ValueError("merge_start_scans 需要至少一个起点扫描")
    x0 = scans[0]
    n = len(x0.energy)
    n_atoms = x0.positions.shape[1]
    n_scans = len(scans)
    for sc in scans[1:]:
        if len(sc.energy) != n or sc.positions.shape != x0.positions.shape:
            raise ValueError("各起点的扫描网格与原子数必须一致（grid_origin_rad 对齐）")

    # 统一候选视图 (scan_idx, direction, s)
    c_e = np.full((n_scans, 2, n), np.nan)
    c_c = np.zeros((n_scans, 2, n), dtype=bool)
    c_p = np.full((n_scans, 2, n, n_atoms, 3), np.nan)
    c_t = np.full((n_scans, 2, n), np.nan)
    c_phi = np.full((n_scans, 2, n), np.nan)
    for si, sc in enumerate(scans):
        if sc.candidates is not None:
            c_e[si] = sc.candidates["energy"]
            c_c[si] = sc.candidates["converged"]
            c_p[si] = sc.candidates["positions"]
            c_t[si] = sc.candidates["phi_target_rad"]
            c_phi[si] = sc.candidates["phi_rad"]
        else:  # 旧缓存：只有逐点最优（方向 0 占位）
            c_e[si, 0] = sc.energy
            c_c[si, 0] = sc.converged
            c_p[si, 0] = sc.positions
            c_t[si, 0] = sc.phi_target_rad
            c_phi[si, 0] = sc.phi_rad

    bond_idx = np.asarray(bonds, dtype=int)
    energy = np.empty(n)
    phi_t = np.empty(n)
    phi = np.empty(n)
    pos = np.empty((n, n_atoms, 3))
    hyst = np.zeros(n)
    converged = np.empty(n, dtype=bool)
    branch_counts = np.ones(n, dtype=int)
    branch_gap = np.zeros(n)
    winner_scan = np.zeros(n, dtype=int)
    r_all = np.empty((n, len(bond_idx)))
    r0_by_scan = np.linalg.norm(
        np.stack([sc.positions[0] for sc in scans])[:, bond_idx[:, 0]]
        - np.stack([sc.positions[0] for sc in scans])[:, bond_idx[:, 1]], axis=2)

    for s in range(n):
        cands = []
        for si in range(n_scans):
            for d in range(2):
                if np.isfinite(c_e[si, d, s]):
                    cands.append((si, d, float(c_e[si, d, s]), bool(c_c[si, d, s]),
                                  c_p[si, d, s], float(c_t[si, d, s]), float(c_phi[si, d, s])))
        pool = [c for c in cands if c[3]] or cands
        if len(pool) > 1:
            ids = greedy_cluster(np.stack([c[4] for c in pool]),
                                 np.array([c[2] for c in pool]), heavy_mask, threshold_nm)
            ncl = int(ids.max()) + 1
            mins = [min(c[2] for c, i in zip(pool, ids) if i == cl) for cl in range(ncl)]
            branch_counts[s] = ncl
            branch_gap[s] = float(max(mins) - min(mins))
        win = min(pool, key=lambda c: (c[2], c[0], c[1]))
        si = win[0]
        energy[s], pos[s] = win[2], win[4]
        phi_t[s], phi[s], converged[s] = win[5], win[6], win[3]
        winner_scan[s] = si
        other = c_e[si, 1 - win[1], s]
        hyst[s] = abs(float(other) - win[2]) if np.isfinite(other) else float(scans[si].hysteresis[s])
        r_all[s] = np.linalg.norm(pos[s][bond_idx[:, 0]] - pos[s][bond_idx[:, 1]], axis=1)

    max_bond_change_pm = float(np.max(np.abs(r_all - r0_by_scan[winner_scan])) * 100.0)
    envelope = ScanResult(
        quartet=x0.quartet,
        phi_target_rad=phi_t,
        phi_rad=phi,
        energy=energy,
        positions=pos,
        hysteresis=hyst,
        converged=converged,
        max_bond_change_pm=max_bond_change_pm,
        params={**x0.params, "n_starts": n_scans, "threshold_nm": float(threshold_nm)},
    )
    delta = energy - np.asarray(x0.energy, dtype=np.float64)
    diagnostics = {
        "threshold_nm": float(threshold_nm),
        "n_points": int(n),
        "branch_count_max": int(branch_counts.max()),
        "branch_count_mean": float(branch_counts.mean()),
        "n_points_multibranch": int((branch_counts > 1).sum()),
        "branch_gap_max": float(branch_gap.max()),
        "envelope_gain_max": float((np.asarray(x0.energy) - energy).max()),
        "n_points_improved": int((delta < -1e-12).sum()),
        "unconverged_points": [int(s) for s in range(n) if not converged[s]],
        "branch_counts": branch_counts.tolist(),
        "branch_gap_kJ_mol": branch_gap.tolist(),
        "envelope_delta_kJ_mol": delta.tolist(),
    }
    return envelope, diagnostics
