"""numpy 参考能量：严格按 MVP-0 计划 §2.1 的锁定函数形式从 IR 求值。

这是所有引擎的共同参照，用来抓住各导出器共有的约定错误。
单位：nm、kJ/mol、rad、e。
"""

from __future__ import annotations

import math

import numpy as np

from .schema import ForceFieldIR

COULOMB_KJ_NM = 138.935456

ENERGY_KEYS = ("bond", "angle", "proper", "improper", "lj", "coulomb", "lj14", "coulomb14")


def dihedral_angle(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """计划 §2.1 锁定的二面角定义：顺式为 0，范围 (-π, π]。

    b1 = r_j - r_i, b2 = r_k - r_j, b3 = r_l - r_k
    φ = atan2(|b2| b1·(b2×b3), (b1×b2)·(b2×b3))
    """
    b1 = np.asarray(p1, dtype=np.float64) - np.asarray(p0, dtype=np.float64)
    b2 = np.asarray(p2, dtype=np.float64) - np.asarray(p1, dtype=np.float64)
    b3 = np.asarray(p3, dtype=np.float64) - np.asarray(p2, dtype=np.float64)
    b2xb3 = np.cross(b2, b3)
    b1xb2 = np.cross(b1, b2)
    y = np.linalg.norm(b2) * np.dot(b1, b2xb3)
    x = np.dot(b1xb2, b2xb3)
    return math.atan2(y, x)


def _bond_angle(p_i: np.ndarray, p_j: np.ndarray, p_k: np.ndarray) -> float:
    """顶点在 j 的键角，atan2 形式，范围 [0, π]。"""
    u1 = np.asarray(p_i, dtype=np.float64) - np.asarray(p_j, dtype=np.float64)
    u2 = np.asarray(p_k, dtype=np.float64) - np.asarray(p_j, dtype=np.float64)
    return math.atan2(np.linalg.norm(np.cross(u1, u2)), np.dot(u1, u2))


def _pair_lj(sigma_ij: float, epsilon_ij: float, r: float) -> float:
    sr6 = (sigma_ij / r) ** 6
    return 4.0 * epsilon_ij * (sr6 * sr6 - sr6)


def evaluate(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    box_nm: np.ndarray | tuple | None = None,
) -> dict[str, float]:
    """返回 ENERGY_KEYS 各项和 "total"，单位 kJ/mol。

    ``box_nm``（3×3，行 = 晶胞向量）非 None 时非键与 1-4 距离用最近镜像
    （决策 0023）；缺省自动取 ``ir.box_nm``。成键项不用镜像——IR 坐标约定
    为分子内连续坐标。
    """
    pos = np.asarray(positions_nm, dtype=np.float64)
    if pos.shape != (len(ir.atoms), 3):
        raise ValueError(f"positions shape {pos.shape} != ({len(ir.atoms)}, 3)")
    if not np.all(np.isfinite(pos)):
        raise ValueError("positions contain non-finite coordinates")
    box = np.asarray(box_nm, dtype=np.float64) if box_nm is not None else (
        np.asarray(ir.box_nm, dtype=np.float64) if ir.box_nm is not None else None
    )

    def mirror(delta: np.ndarray) -> np.ndarray:
        frac = np.linalg.solve(box.T, delta)
        frac -= np.round(frac)
        return box.T @ frac

    if box is not None:
        # 周期体系：先沿键图把（可能被 wrap 打散的）坐标重整为连续分子坐标，
        # 使键/角/二面恢复真实值；对已是连续的坐标此步为恒等操作
        from .topology import bond_graph

        graph = bond_graph(ir)
        unwrapped = pos.copy()
        visited = {0}
        stack = [0]
        while stack:
            a = stack.pop()
            for b in graph.get(a, ()):  # type: ignore[union-attr]
                if b not in visited:
                    unwrapped[b] = unwrapped[a] + mirror(unwrapped[b] - unwrapped[a])
                    visited.add(b)
                    stack.append(b)
        pos = unwrapped

    def pair_distance(i: int, j: int) -> float:
        delta = pos[j] - pos[i]
        if box is not None:
            delta = mirror(delta)
        return float(np.linalg.norm(delta))

    energy = {key: 0.0 for key in ENERGY_KEYS}

    for bond in ir.bonds:
        r = float(np.linalg.norm(pos[bond.j] - pos[bond.i]))
        energy["bond"] += 0.5 * bond.k * (r - bond.r0) ** 2

    for angle in ir.angles:
        theta = _bond_angle(pos[angle.i], pos[angle.j], pos[angle.k_atom])
        energy["angle"] += 0.5 * angle.force_constant * (theta - angle.theta0) ** 2

    for term in ir.proper_torsions:
        phi = dihedral_angle(pos[term.i], pos[term.j], pos[term.k], pos[term.l])
        energy["proper"] += term.force_constant * (1.0 + math.cos(term.periodicity * phi - term.phase))

    for term in ir.impropers:
        phi = dihedral_angle(pos[term.i], pos[term.j], pos[term.k], pos[term.l])
        energy["improper"] += term.force_constant * (1.0 + math.cos(term.periodicity * phi - term.phase))

    excluded = {tuple(p) for p in ir.exclusions}
    scaled = {tuple(p) for p in ir.pairs14}
    n = len(ir.atoms)
    for i in range(n):
        atom_i = ir.atoms[i]
        for j in range(i + 1, n):
            pair = (i, j)
            if pair in excluded:
                continue
            atom_j = ir.atoms[j]
            r = pair_distance(i, j)
            if r == 0.0:
                raise ValueError(f"atoms {i} and {j} coincide; nonbonded energy diverges")
            sigma_ij = 0.5 * (atom_i.sigma + atom_j.sigma)
            epsilon_ij = math.sqrt(atom_i.epsilon * atom_j.epsilon)
            lj = _pair_lj(sigma_ij, epsilon_ij, r)
            coulomb = COULOMB_KJ_NM * atom_i.partial_charge * atom_j.partial_charge / r
            if pair in scaled:
                energy["lj14"] += ir.vdw_14_scale * lj
                energy["coulomb14"] += ir.coulomb_14_scale * coulomb
            else:
                energy["lj"] += lj
                energy["coulomb"] += coulomb

    energy["total"] = sum(energy[key] for key in ENERGY_KEYS)
    return energy


def forces_fd(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    step_nm: float = 1e-6,
) -> np.ndarray:
    """evaluate 总能量的中心差分力，(N,3) kJ/mol/nm。"""
    pos = np.asarray(positions_nm, dtype=np.float64)
    forces = np.zeros_like(pos)
    for i in range(pos.shape[0]):
        for axis in range(3):
            plus = pos.copy()
            minus = pos.copy()
            plus[i, axis] += step_nm
            minus[i, axis] -= step_nm
            e_plus = evaluate(ir, plus)["total"]
            e_minus = evaluate(ir, minus)["total"]
            forces[i, axis] = -(e_plus - e_minus) / (2.0 * step_nm)
    return forces
