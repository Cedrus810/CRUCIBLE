"""§C3 液体系统装配（validation/liquid.py）的单元测试：pack_box、to_liquid_system、liquid_dynamics。

用 ethanol 固件（9 原子、q=0、extent ~0.27 nm）；只做装配、单点与最小化+两步积分，
不做采样协议（平衡/时长/块平均属运行计划，另行报资源后执行）。参考平台双精度。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import unit

from ffcompiler.importers import from_amber
from ffcompiler.ir.schema import ForceFieldIR
from ffcompiler.validation.liquid import pack_box, to_liquid_system, liquid_dynamics

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ETHANOL = FIXTURES / "ethanol"
N_MOL = 27
RHO = 0.789  # g/cm³，乙醇 298 K


@pytest.fixture(scope="module")
def ethanol():
    ir, pos = from_amber(*(ETHANOL / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    return ir, pos


@pytest.fixture(scope="module")
def packed(ethanol):
    ir, pos = ethanol
    return pack_box(ir, pos, n_molecules=N_MOL, density_g_cm3=RHO, seed=3)


def test_pack_box_topology_and_density(packed, ethanol):
    boxed, positions = packed
    ir, _ = ethanol
    assert len(boxed.atoms) == N_MOL * len(ir.atoms)
    assert len(boxed.bonds) == N_MOL * len(ir.bonds)
    assert len(boxed.exclusions) == N_MOL * len(ir.exclusions)
    assert len(boxed.pairs14) == N_MOL * len(ir.pairs14)
    assert boxed.total_charge == 0
    box = np.asarray(boxed.box_nm)
    assert np.allclose(box, np.diag([box[0, 0]] * 3))
    # 密度回算（独立于模块常数的物理换算）：m/V
    mass_g = sum(a.mass for a in boxed.atoms) * 1.66053906660e-24
    volume_cm3 = (box[0, 0] * 1e-7) ** 3
    assert mass_g / volume_cm3 == pytest.approx(RHO, rel=1e-6)
    assert positions.shape == (N_MOL * len(ir.atoms), 3)
    assert boxed.provenance["pack"]["n_molecules"] == N_MOL
    assert boxed.provenance["pack"]["density_g_cm3"] == RHO


def test_pack_box_deterministic_and_seeded(ethanol):
    ir, pos = ethanol
    kwargs = dict(n_molecules=N_MOL, density_g_cm3=RHO)
    a_ir, a_pos = pack_box(ir, pos, seed=3, **kwargs)
    b_ir, b_pos = pack_box(ir, pos, seed=3, **kwargs)
    c_pos = pack_box(ir, pos, seed=4, **kwargs)[1]
    assert np.array_equal(a_pos, b_pos)
    assert not np.array_equal(a_pos, c_pos)
    assert a_ir.provenance["pack"] == b_ir.provenance["pack"]


def test_pack_box_intermolecular_gap_guarantee(packed):
    """格点 + 抖动的几何保证：不同分子原子间距 ≥ 2×min_gap（摆放坐标逐轴论证，
    不做逐分子居中——同构象分子居中后会互相重合）。"""
    _, positions = packed
    n_atoms = 9
    flat = positions
    mol_of = np.repeat(np.arange(N_MOL), n_atoms)
    diff = flat[:, None, :] - flat[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    mask = mol_of[:, None] != mol_of[None, :]
    np.fill_diagonal(mask, False)
    assert dist[mask].min() >= 2 * 0.05 - 1e-9


def test_pack_box_rejects_impossible_spacing(ethanol):
    ir, pos = ethanol
    with pytest.raises(ValueError, match="放不下"):
        # 高密度 → 格点间距小于 extent + 2×min_gap
        pack_box(ir, pos, n_molecules=8, density_g_cm3=20.0)


def test_pack_box_rejects_single_molecule(ethanol):
    ir, pos = ethanol
    with pytest.raises(ValueError, match="至少 2 个分子"):
        pack_box(ir, pos, n_molecules=1, density_g_cm3=RHO)


def test_to_liquid_system_pme_and_exceptions(packed):
    boxed, positions = packed
    system = to_liquid_system(boxed, cutoff_nm=0.6)
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    assert nb.getNonbondedMethod() == openmm.NonbondedForce.PME
    assert nb.getCutoffDistance().value_in_unit(unit.nanometer) == pytest.approx(0.6)
    assert nb.getUseDispersionCorrection()
    a, b, c = system.getDefaultPeriodicBoxVectors()
    box = np.asarray(boxed.box_nm)
    assert (a.x, b.y, c.z) == (box[0, 0], box[1, 1], box[2, 2])
    # exceptions 只有分子内 1-2/1-3（零强度）+ 1-4（缩放），无跨分子 exception
    assert nb.getNumExceptions() == len(boxed.exclusions) + len(boxed.pairs14)
    assert system.getNumParticles() == N_MOL * 9
    # 参考平台单点能量有限且为负（液体有吸引相互作用）
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001),
                         openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions([openmm.Vec3(*p) for p in positions] * unit.nanometer)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del ctx
    assert np.isfinite(e)  # 符号不作判据：未最小化的格点初始构型近接触为正是预期


def test_to_liquid_system_requires_box_orthogonal_and_cutoff(ethanol, packed):
    ir, _ = ethanol
    with pytest.raises(ValueError, match="周期盒"):
        to_liquid_system(ir)
    boxed, _ = packed
    box = np.asarray(boxed.box_nm)
    sheared = ForceFieldIR(
        name=boxed.name, atoms=boxed.atoms, bonds=boxed.bonds, angles=boxed.angles,
        proper_torsions=boxed.proper_torsions, impropers=boxed.impropers,
        exclusions=boxed.exclusions, pairs14=boxed.pairs14,
        total_charge=boxed.total_charge, forcefield_family=boxed.forcefield_family,
        combination_rule=boxed.combination_rule,
        coulomb_14_scale=boxed.coulomb_14_scale, vdw_14_scale=boxed.vdw_14_scale,
        box_nm=((box[0, 0], 0.0, 0.0), (0.5, box[1, 1], 0.0), (0.0, 0.0, box[2, 2])),
        provenance=boxed.provenance)
    with pytest.raises(ValueError, match="正交盒"):
        to_liquid_system(sheared)
    with pytest.raises(ValueError, match="最小镜像"):
        to_liquid_system(boxed, cutoff_nm=0.9)  # 2×0.9 > L≈1.378


def test_liquid_dynamics_minimize_and_two_steps(packed):
    """装配 + 最小化 + 两步积分（能量有限、周期盒存在）；不做平衡/采样协议。"""
    boxed, positions = packed
    sim = liquid_dynamics(boxed, positions, platform_name="CPU", cutoff_nm=0.6)
    sim.minimizeEnergy(maxIterations=20)
    sim.step(2)
    state = sim.context.getState(getEnergy=True)
    assert np.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    assert sim.context.getState().getPeriodicBoxVolume() > 0.0 * unit.nanometer**3
