"""液体基准的系统装配（2026-11-20 计划 §C3；设计 §20 Level 3 首次立项的工程前置）。

现有 ``exporters.openmm.to_openmm_system`` 用 ``NoCutoff``，周期盒只作元数据——那是
IR 等价性验证的语义（决策 0023），**保持不动**。本模块另建 MD 装配：

- :func:`pack_box` —— 单分子 IR + 均匀密度 → 周期盒子 IR（立方格点摆放 + 抖动，
  无旋转 v1；确定性：固定 seed）；
- :func:`to_liquid_system` —— PME + LJ 截断 + 色散修正的 OpenMM System
  （exceptions 语义与导出器一致：1-2/1-3 零强度、1-4 按 IR 缩放）；
- :func:`liquid_dynamics` —— System + LangevinMiddle 积分器 + MonteCarlo 恒压器
  （NPT），返回可直接 minimize/run 的 ``openmm.app.Simulation``。

ΔHvap 的气相腿用现有 ``to_openmm_system``（NoCutoff）即可，不在本模块。
协议参数（平衡/采样时长、盒数）属运行计划，按项目惯例另行报资源后执行。
"""

from __future__ import annotations

import math
from dataclasses import replace as dc_replace

import numpy as np

from ..exporters.openmm import to_openmm_topology
from ..ir.schema import ForceFieldIR
from ..ir.validation import validate

AMU_PER_NM3_TO_G_CM3 = 1.66053906660e-3  # 1 amu/nm³ = 1.66054e-3 g/cm³（1 g/cm³ = 602.21 amu/nm³）


def _extent_nm(positions_nm: np.ndarray) -> np.ndarray:
    """最大原子间距（各维）：粗略分子尺寸上界（不加缓冲）。"""
    return np.max(positions_nm, axis=0) - np.min(positions_nm, axis=0)


def pack_box(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    *,
    n_molecules: int,
    density_g_cm3: float,
    seed: int = 0,
    min_gap_nm: float = 0.05,
) -> tuple[ForceFieldIR, np.ndarray]:
    """单分子 → 周期盒子（立方、格点 + 抖动、无旋转 v1）。

    盒长 L 由目标密度决定；格点间距放不下分子（extent + 2×min_gap > 间距）时报错，
    不静默降密度。返回 ``(boxed_ir, positions)``：boxed_ir 原子/项数 = n 倍单分子、
    ``box_nm`` = diag(L,L,L)，provenance 记 pack 块。确定性：同 seed 逐位相同。
    """
    validate(ir)
    pos = np.asarray(positions_nm, dtype=np.float64)
    if pos.shape != (len(ir.atoms), 3):
        raise ValueError(f"positions {pos.shape} 与原子数 {len(ir.atoms)} 不符")
    if n_molecules < 2:
        raise ValueError("液体盒子至少 2 个分子（1 个没有界面统计意义）")
    mass_amu = sum(a.mass for a in ir.atoms)
    volume_nm3 = n_molecules * mass_amu / (density_g_cm3 / AMU_PER_NM3_TO_G_CM3)
    length = volume_nm3 ** (1.0 / 3.0)

    n_side = math.ceil(round(n_molecules ** (1.0 / 3.0), 9))
    while n_side ** 3 < n_molecules:
        n_side += 1
    spacing = length / n_side
    extent = _extent_nm(pos)
    margin = extent + 2.0 * min_gap_nm
    if np.any(margin > spacing):
        raise ValueError(
            f"格点间距 {spacing:.3f} nm 放不下分子（extent {np.round(extent, 3).tolist()} nm"
            f" + 2×{min_gap_nm} nm 间隙）；减少分子数或提高密度/换更大的盒子协议")

    # 中心化后放格点，格心剩余空间内均匀抖动（不超过 (spacing−margin)/2，保证 ≥ min_gap）
    centered = pos - 0.5 * (np.max(pos, axis=0) + np.min(pos, axis=0))
    slack = np.maximum(spacing - margin, 0.0) * 0.5
    rng = np.random.default_rng(seed)
    offset = 0.5 * (n_side - 1) * spacing
    atoms = []
    bonds, angles, propers, impropers = [], [], [], []
    exclusions, pairs14 = [], []
    positions = np.empty((n_molecules, len(ir.atoms), 3))
    placed = 0
    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                if placed >= n_molecules:
                    break
                shift = (np.array([ix, iy, iz]) * spacing - offset
                         + rng.uniform(-slack, slack))
                positions[placed] = centered + shift
                atoms.extend(
                    dc_replace(a, index=a.index + placed * len(ir.atoms))
                    for a in ir.atoms)
                off = placed * len(ir.atoms)
                bonds.extend(dc_replace(b, i=b.i + off, j=b.j + off) for b in ir.bonds)
                angles.extend(dc_replace(a, i=a.i + off, j=a.j + off, k_atom=a.k_atom + off)
                              for a in ir.angles)
                propers.extend(dc_replace(t, i=t.i + off, j=t.j + off, k=t.k + off, l=t.l + off)
                               for t in ir.proper_torsions)
                impropers.extend(dc_replace(t, i=t.i + off, j=t.j + off, k=t.k + off, l=t.l + off)
                                 for t in ir.impropers)
                exclusions.extend((i + off, j + off) for i, j in ir.exclusions)
                pairs14.extend((i + off, j + off) for i, j in ir.pairs14)
                placed += 1

    boxed = ir.replace_terms(
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        angles=tuple(angles),
        proper_torsions=tuple(propers),
        impropers=tuple(impropers),
        provenance_update={"pack": {
            "n_molecules": int(n_molecules),
            "density_g_cm3": float(density_g_cm3),
            "seed": int(seed),
            "min_gap_nm": float(min_gap_nm),
            "box_length_nm": float(length),
            "rotation": "none",
            "note": "§C3 v1：立方格点 + 抖动摆放；平衡/采样协议另行报资源后执行",
        }})
    # replace_terms 不携带 box_nm（决策 0023 语义），盒子显式重建
    boxed = ForceFieldIR(
        name=boxed.name, atoms=boxed.atoms, bonds=boxed.bonds, angles=boxed.angles,
        proper_torsions=boxed.proper_torsions, impropers=boxed.impropers,
        exclusions=tuple(exclusions), pairs14=tuple(pairs14),
        total_charge=ir.total_charge * n_molecules,
        forcefield_family=ir.forcefield_family, combination_rule=ir.combination_rule,
        coulomb_14_scale=ir.coulomb_14_scale, vdw_14_scale=ir.vdw_14_scale,
        box_nm=((length, 0.0, 0.0), (0.0, length, 0.0), (0.0, 0.0, length)),
        provenance=boxed.provenance)
    validate(boxed)
    return boxed, positions.reshape(-1, 3)


def to_liquid_system(
    ir: ForceFieldIR,
    *,
    cutoff_nm: float = 1.0,
    dispersion_correction: bool = True,
) -> "openmm.System":
    """PME + LJ 截断 + 色散修正的 System（MD 装配，与导出器 NoCutoff 语义分离）。

    exceptions 与 ``exporters.openmm`` 生产模式一致：1-2/1-3 零强度、1-4 按
    ``coulomb_14_scale``/``vdw_14_scale``。无盒或非正交盒报错（pack_box 产物是立方盒）。
    """
    import openmm

    if ir.box_nm is None:
        raise ValueError("液体系统需要周期盒：先 pack_box（或给 IR 配 box_nm）")
    b = np.asarray(ir.box_nm, dtype=np.float64)
    if abs(b[0, 1]) > 1e-12 or abs(b[0, 2]) > 1e-12 or abs(b[1, 0]) > 1e-12 \
            or abs(b[1, 2]) > 1e-12 or abs(b[2, 0]) > 1e-12 or abs(b[2, 1]) > 1e-12:
        raise ValueError("to_liquid_system v1 只支持正交盒（pack_box 产物）")
    if not 2.0 * cutoff_nm < min(b[0, 0], b[1, 1], b[2, 2]):
        raise ValueError(
            f"截断 {cutoff_nm} nm 不满足最小镜像约定（盒长 {b[0, 0]:.3f} nm 需 > 2×截断）")

    system = openmm.System()
    for atom in ir.atoms:
        system.addParticle(atom.mass)
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(b[0, 0], 0.0, 0.0),
        openmm.Vec3(0.0, b[1, 1], 0.0),
        openmm.Vec3(0.0, 0.0, b[2, 2]),
    )

    bond_force = openmm.HarmonicBondForce()
    for bond in ir.bonds:
        bond_force.addBond(bond.i, bond.j, bond.r0, bond.k)
    system.addForce(bond_force)
    angle_force = openmm.HarmonicAngleForce()
    for a in ir.angles:
        angle_force.addAngle(a.i, a.j, a.k_atom, a.theta0, a.force_constant)
    system.addForce(angle_force)
    proper_force = openmm.PeriodicTorsionForce()
    for t in ir.proper_torsions:
        proper_force.addTorsion(t.i, t.j, t.k, t.l, t.periodicity, t.phase, t.force_constant)
    system.addForce(proper_force)
    improper_force = openmm.PeriodicTorsionForce()
    for t in ir.impropers:
        improper_force.addTorsion(t.i, t.j, t.k, t.l, t.periodicity, t.phase, t.force_constant)
    system.addForce(improper_force)

    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(cutoff_nm)
    nonbonded.setUseDispersionCorrection(bool(dispersion_correction))
    for atom in ir.atoms:
        nonbonded.addParticle(atom.partial_charge, atom.sigma, atom.epsilon)
    exclusions = {tuple(p) for p in ir.exclusions}
    pairs14 = {tuple(p) for p in ir.pairs14}
    for i, j in sorted(exclusions):
        nonbonded.addException(i, j, 0.0, 1.0, 0.0)
    for i, j in sorted(pairs14):
        ai, aj = ir.atoms[i], ir.atoms[j]
        nonbonded.addException(
            i, j,
            ir.coulomb_14_scale * ai.partial_charge * aj.partial_charge,
            0.5 * (ai.sigma + aj.sigma),
            ir.vdw_14_scale * math.sqrt(ai.epsilon * aj.epsilon),
        )
    system.addForce(nonbonded)
    return system


def liquid_dynamics(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    *,
    temperature_k: float = 298.15,
    pressure_bar: float = 1.0,
    dt_fs: float = 2.0,
    friction_per_ps: float = 1.0,
    barostat_interval_steps: int = 25,
    platform_name: str = "CPU",
    cutoff_nm: float = 1.0,
) -> "openmm.app.Simulation":
    """NPT Simulation（LangevinMiddle + MonteCarloBarostat + PME）。

    平台默认 CPU（项目机器无 GPU 约定）；返回的 Simulation 尚未 minimize，
    运行协议（最小化/平衡/采样时长与块平均）由调用方按已批资源执行。
    """
    import openmm
    import openmm.app
    from openmm import unit

    system = to_liquid_system(ir, cutoff_nm=cutoff_nm)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin, friction_per_ps / unit.picosecond,
        dt_fs * unit.femtosecond)
    barostat = openmm.MonteCarloBarostat(
        pressure_bar * unit.bar, temperature_k * unit.kelvin, barostat_interval_steps)
    system.addForce(barostat)
    topology = to_openmm_topology(ir)
    if ir.box_nm is not None:
        box = np.asarray(ir.box_nm, dtype=np.float64)
        topology.setPeriodicBoxVectors((
            openmm.Vec3(box[0, 0], box[0, 1], box[0, 2]),
            openmm.Vec3(box[1, 0], box[1, 1], box[1, 2]),
            openmm.Vec3(box[2, 0], box[2, 1], box[2, 2]),
        ))
    simulation = openmm.app.Simulation(
        topology, system, integrator,
        openmm.Platform.getPlatformByName(platform_name))
    simulation.context.setPositions(
        [openmm.Vec3(*p) for p in np.asarray(positions_nm, dtype=np.float64)]
        * unit.nanometer)
    return simulation
