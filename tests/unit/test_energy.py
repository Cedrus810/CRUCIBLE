"""T3：numpy 参考能量的必测用例。"""

from __future__ import annotations

import math

import numpy as np
import openmm
from openmm import unit

from ffcompiler.ir import Angle, Atom, Bond, ForceFieldIR, ProperTorsion
from ffcompiler.ir.energy import COULOMB_KJ_NM, ENERGY_KEYS, dihedral_angle, evaluate
from ffcompiler.ir.topology import expected_exclusions_and_pairs
from ffcompiler.ir.validation import validate
from ffcompiler.validation.engines import reference_context


def _atom(i, element="C", z=6, mass=12.0, q=0.0, sigma=0.3, epsilon=0.4, atype=None):
    return Atom(
        index=i, name=f"A{i}", element=element, atomic_number=z, mass=mass,
        atom_type=atype or f"T{i}", formal_charge=0, partial_charge=q,
        sigma=sigma, epsilon=epsilon,
    )


def test_two_atom_bond_energy_exact() -> None:
    ir = ForceFieldIR(
        name="dimer", atoms=(_atom(0), _atom(1)),
        bonds=(Bond(i=0, j=1, r0=0.1, k=1000.0, source="t"),),
        angles=(), proper_torsions=(), impropers=(),
        exclusions=((0, 1),), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )
    pos = np.array([[0.0, 0.0, 0.0], [0.11, 0.0, 0.0]])
    out = evaluate(ir, pos)
    assert abs(out["bond"] - 0.5 * 1000.0 * 0.01 ** 2) <= 1e-12  # = 0.05 kJ/mol
    assert out["total"] == out["bond"]


def test_three_atom_angle_energy() -> None:
    # θ = 100°, θ0 = 90°, k = 500 → E = ½·500·(π/18)²
    ir = ForceFieldIR(
        name="trimer",
        atoms=(_atom(0), _atom(1), _atom(2)),
        bonds=(),  # evaluate 不依赖 bond 列表
        angles=(Angle(i=0, j=1, k_atom=2, theta0=math.pi / 2, force_constant=500.0, source="t"),),
        proper_torsions=(), impropers=(),
        exclusions=((0, 1), (0, 2), (1, 2)), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )
    theta = math.radians(100.0)
    pos = np.array([
        [math.cos(theta), math.sin(theta), 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    out = evaluate(ir, pos)
    expected = 0.5 * 500.0 * (math.pi / 18.0) ** 2
    assert abs(out["angle"] - expected) <= 1e-12


def test_dihedral_definitions_and_openmm_sign() -> None:
    s3 = math.sqrt(3) / 2
    i = np.array([-0.5, -s3, 0.0])
    j = np.array([0.0, 0.0, 0.0])
    k = np.array([1.0, 0.0, 0.0])
    l_cis = np.array([1.5, -s3, 0.0])
    l_trans = np.array([1.5, s3, 0.0])
    # 顺式 = 0，反式 = π
    assert abs(dihedral_angle(i, j, k, l_cis)) <= 1e-12
    assert abs(abs(dihedral_angle(i, j, k, l_trans)) - math.pi) <= 1e-12

    # 绕 b2（j→k，+x 方向）右手旋转 +60° → φ = +π/3
    c, s = math.cos(math.pi / 3), math.sin(math.pi / 3)
    y, z = l_cis[1], l_cis[2]
    l_60 = np.array([l_cis[0], y * c - z * s, y * s + z * c])
    phi = dihedral_angle(i, j, k, l_60)
    assert abs(phi - math.pi / 3) <= 1e-12

    # 与 OpenMM 符号一致：PeriodicTorsionForce(n=1, δ=0, k=1) → E = 1 + cos φ
    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0)
    torsion = openmm.PeriodicTorsionForce()
    torsion.addTorsion(0, 1, 2, 3, 1, 0.0, 1.0)
    system.addForce(torsion)
    context = reference_context(system)
    for name, geo in (("cis", l_cis), ("trans", l_trans), ("plus60", l_60)):
        pos = np.vstack([i, j, k, geo])
        context.setPositions([openmm.Vec3(*p) for p in pos] * unit.nanometer)
        e_omm = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        e_ref = 1.0 + math.cos(dihedral_angle(pos[0], pos[1], pos[2], pos[3]))
        assert abs(e_omm - e_ref) <= 1e-9, name


def test_lj_coulomb_classification_and_scaling() -> None:
    def two_atom_ir(exclusions, pairs14):
        return ForceFieldIR(
            name="pair",
            atoms=(
                _atom(0, q=0.7, sigma=0.3, epsilon=0.5),
                _atom(1, element="O", z=8, mass=16.0, q=-0.7, sigma=0.3, epsilon=0.5),
            ),
            bonds=(), angles=(), proper_torsions=(), impropers=(),
            exclusions=exclusions, pairs14=pairs14, total_charge=0,
            forcefield_family="amber", combination_rule="lorentz-berthelot",
            coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
        )

    r = 0.25
    pos = np.array([[0.0, 0.0, 0.0], [r, 0.0, 0.0]])
    sr6 = (0.3 / r) ** 6
    lj_hand = 4.0 * 0.5 * (sr6 * sr6 - sr6)
    q_hand = COULOMB_KJ_NM * 0.7 * (-0.7) / r

    free = evaluate(two_atom_ir((), ()), pos)
    assert abs(free["lj"] - lj_hand) <= 1e-9
    assert abs(free["coulomb"] - q_hand) <= 1e-9
    assert free["lj14"] == 0.0 and free["coulomb14"] == 0.0

    scaled = evaluate(two_atom_ir((), ((0, 1),)), pos)
    assert abs(scaled["lj14"] - 0.5 * lj_hand) <= 1e-9
    assert abs(scaled["coulomb14"] - (1.0 / 1.2) * q_hand) <= 1e-9
    assert scaled["lj"] == 0.0 and scaled["coulomb"] == 0.0

    excluded = evaluate(two_atom_ir(((0, 1),), ()), pos)
    assert excluded["lj"] == 0.0 and excluded["coulomb"] == 0.0
    assert excluded["lj14"] == 0.0 and excluded["coulomb14"] == 0.0


def _chain_ir() -> ForceFieldIR:
    atoms = (
        _atom(0, q=-0.3, sigma=0.33, epsilon=0.45),
        _atom(1, q=0.1, sigma=0.25, epsilon=0.3),
        _atom(2, q=-0.1, sigma=0.35, epsilon=0.5),
        _atom(3, q=0.3, sigma=0.2, epsilon=0.25),
    )
    ir = ForceFieldIR(
        name="chain4", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.15, k=250000.0, source="t"),
            Bond(i=1, j=2, r0=0.153, k=300000.0, source="t"),
            Bond(i=2, j=3, r0=0.148, k=280000.0, source="t"),
        ),
        angles=(
            Angle(i=0, j=1, k_atom=2, theta0=1.91, force_constant=520.0, source="t"),
            Angle(i=1, j=2, k_atom=3, theta0=1.88, force_constant=480.0, source="t"),
        ),
        proper_torsions=(
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=1, phase=0.0,
                          force_constant=1.2, source="t"),
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=3, phase=math.pi,
                          force_constant=0.5, source="t"),
        ),
        impropers=(), exclusions=(), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    exclusions, pairs14 = expected_exclusions_and_pairs(ir)
    ir = ForceFieldIR(
        name=ir.name, atoms=ir.atoms, bonds=ir.bonds, angles=ir.angles,
        proper_torsions=ir.proper_torsions, impropers=ir.impropers,
        exclusions=exclusions, pairs14=pairs14, total_charge=ir.total_charge,
        forcefield_family=ir.forcefield_family,
        combination_rule=ir.combination_rule,
        coulomb_14_scale=ir.coulomb_14_scale, vdw_14_scale=ir.vdw_14_scale,
    )
    validate(ir)
    return ir


def _chain_openmm_system(ir: ForceFieldIR) -> openmm.System:
    system = openmm.System()
    for atom in ir.atoms:
        system.addParticle(atom.mass)
    bond_force = openmm.HarmonicBondForce()
    for b in ir.bonds:
        bond_force.addBond(b.i, b.j, b.r0, b.k)
    system.addForce(bond_force)

    angle_force = openmm.HarmonicAngleForce()
    for a in ir.angles:
        angle_force.addAngle(a.i, a.j, a.k_atom, a.theta0, a.force_constant)
    system.addForce(angle_force)

    for t in ir.proper_torsions:
        torsion_force = openmm.PeriodicTorsionForce()
        torsion_force.addTorsion(t.i, t.j, t.k, t.l, t.periodicity, t.phase, t.force_constant)
        system.addForce(torsion_force)

    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for atom in ir.atoms:
        nonbonded.addParticle(atom.partial_charge, atom.sigma, atom.epsilon)
    pairs14 = {tuple(p) for p in ir.pairs14}
    excluded = {tuple(p) for p in ir.exclusions}
    for i in range(len(ir.atoms)):
        for j in range(i + 1, len(ir.atoms)):
            if (i, j) in excluded:
                nonbonded.addException(i, j, 0.0, 0.5, 0.0)
            elif (i, j) in pairs14:
                ai, aj = ir.atoms[i], ir.atoms[j]
                nonbonded.addException(
                    i, j,
                    ir.coulomb_14_scale * ai.partial_charge * aj.partial_charge,
                    0.5 * (ai.sigma + aj.sigma),
                    ir.vdw_14_scale * math.sqrt(ai.epsilon * aj.epsilon),
                )
            else:
                raise AssertionError(f"pair ({i},{j}) unexpectedly full-strength")
    system.addForce(nonbonded)
    return system


def test_evaluate_vs_handbuilt_openmm_10_configs() -> None:
    ir = _chain_ir()
    system = _chain_openmm_system(ir)
    integrator = openmm.VerletIntegrator(0.001)
    context = reference_context(system)

    rng = np.random.default_rng(20261004)
    base = np.array([
        [0.0, 0.0, 0.0],
        [0.15, 0.0, 0.0],
        [0.15 + 0.15 * math.cos(2.0), 0.15 * math.sin(2.0), 0.0],
        [0.30 + 0.15 * math.cos(2.0) + 0.15 * math.cos(4.0),
         0.15 * math.sin(2.0) + 0.15 * math.sin(4.0), 0.05],
    ])
    for n in range(10):
        pos = base + rng.normal(0.0, 0.02, size=(4, 3))
        context.setPositions([openmm.Vec3(*p) for p in pos] * unit.nanometer)
        e_omm = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        e_ref = evaluate(ir, pos)["total"]
        rel = abs(e_omm - e_ref) / max(1.0, abs(e_ref))
        assert rel <= 1e-9, f"config {n}: {e_omm} vs {e_ref}"

    assert set(ENERGY_KEYS) == {
        "bond", "angle", "proper", "improper", "lj", "coulomb", "lj14", "coulomb14",
    }
