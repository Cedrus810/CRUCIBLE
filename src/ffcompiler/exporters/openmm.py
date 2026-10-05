"""OpenMM 导出器（MVP-0 计划 T6）。

decompose=False（生产模式）：一个 NonbondedForce 承担全部非键（1-2/1-3 零强度
exception、1-4 缩放 exception），propers / impropers 各一个 PeriodicTorsionForce。
decompose=True：按计划 §2.3 把八类能量分到八个力组，供跨引擎比较。
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import openmm
import openmm.app
from openmm import unit

from ffcompiler.ir import ForceFieldIR, validate
from ffcompiler.ir.energy import COULOMB_KJ_NM

FORCE_GROUPS: dict[str, int] = {
    "bond": 0,
    "angle": 1,
    "proper": 2,
    "improper": 3,
    "lj": 4,
    "coulomb": 5,
    "lj14": 6,
    "coulomb14": 7,
}


def _add_torsion_terms(force: openmm.PeriodicTorsionForce, terms) -> None:
    for t in terms:
        force.addTorsion(t.i, t.j, t.k, t.l, t.periodicity, t.phase, t.force_constant)


def to_openmm_system(
    ir: ForceFieldIR, *, decompose: bool = False, validate_ir: bool = True
) -> openmm.System:
    # validate_ir=False 供 H_fixed 等裁剪 IR 的数值工具使用（T22）：
    # 剥掉 bond 的 IR 无法通过完整 validate，但力仍可正确构造。
    if validate_ir:
        validate(ir)
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

    proper_force = openmm.PeriodicTorsionForce()
    _add_torsion_terms(proper_force, ir.proper_torsions)
    system.addForce(proper_force)

    improper_force = openmm.PeriodicTorsionForce()
    _add_torsion_terms(improper_force, ir.impropers)
    system.addForce(improper_force)

    exclusions = {tuple(p) for p in ir.exclusions}
    pairs14 = {tuple(p) for p in ir.pairs14}

    if decompose:
        # OpenMM 的 CustomNonbondedForce 不自动组合逐粒子参数：
        # 在表达式里显式写 Lorentz–Berthelot（σ 算术平均、ε 几何平均）
        lj_force = openmm.CustomNonbondedForce(
            "4*sqrt(epsilon1*epsilon2)*"
            "((0.5*(sigma1+sigma2)/r)^12-(0.5*(sigma1+sigma2)/r)^6)"
        )
        lj_force.addPerParticleParameter("sigma")
        lj_force.addPerParticleParameter("epsilon")
        for atom in ir.atoms:
            lj_force.addParticle((atom.sigma, atom.epsilon))
        coulomb_force = openmm.CustomNonbondedForce(
            f"{COULOMB_KJ_NM!r}*charge1*charge2/r"
        )
        coulomb_force.addPerParticleParameter("charge")
        for atom in ir.atoms:
            coulomb_force.addParticle((atom.partial_charge,))
        for i, j in sorted(exclusions | pairs14):
            lj_force.addExclusion(i, j)
            coulomb_force.addExclusion(i, j)

        lj14_force = openmm.CustomBondForce("4*epsilon*((sigma/r)^12-(sigma/r)^6)")
        lj14_force.addPerBondParameter("epsilon")
        lj14_force.addPerBondParameter("sigma")
        coulomb14_force = openmm.CustomBondForce(f"{COULOMB_KJ_NM!r}*qprod/r")
        coulomb14_force.addPerBondParameter("qprod")
        for i, j in sorted(pairs14):
            ai, aj = ir.atoms[i], ir.atoms[j]
            lj14_force.addBond(
                i, j,
                (ir.vdw_14_scale * math.sqrt(ai.epsilon * aj.epsilon),
                 0.5 * (ai.sigma + aj.sigma)),
            )
            coulomb14_force.addBond(
                i, j, (ir.coulomb_14_scale * ai.partial_charge * aj.partial_charge,),
            )

        for force in (lj_force, coulomb_force, lj14_force, coulomb14_force):
            system.addForce(force)
        for name, force in (
            ("bond", bond_force), ("angle", angle_force),
            ("proper", proper_force), ("improper", improper_force),
            ("lj", lj_force), ("coulomb", coulomb_force),
            ("lj14", lj14_force), ("coulomb14", coulomb14_force),
        ):
            force.setForceGroup(FORCE_GROUPS[name])
    else:
        nonbonded = openmm.NonbondedForce()
        nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        for atom in ir.atoms:
            nonbonded.addParticle(atom.partial_charge, atom.sigma, atom.epsilon)
        for i in range(len(ir.atoms)):
            for j in range(i + 1, len(ir.atoms)):
                if (i, j) in exclusions:
                    nonbonded.addException(i, j, 0.0, 1.0, 0.0)
                elif (i, j) in pairs14:
                    ai, aj = ir.atoms[i], ir.atoms[j]
                    nonbonded.addException(
                        i, j,
                        ir.coulomb_14_scale * ai.partial_charge * aj.partial_charge,
                        0.5 * (ai.sigma + aj.sigma),
                        ir.vdw_14_scale * math.sqrt(ai.epsilon * aj.epsilon),
                    )
        system.addForce(nonbonded)

    return system


def to_openmm_topology(ir: ForceFieldIR) -> openmm.app.Topology:
    topology = openmm.app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue(ir.atoms[0].residue_name, chain)
    openmm_atoms = []
    for atom in ir.atoms:
        element = openmm.app.element.Element.getByAtomicNumber(atom.atomic_number)
        openmm_atoms.append(topology.addAtom(atom.name, element, residue))
    for bond in ir.bonds:
        topology.addBond(openmm_atoms[bond.i], openmm_atoms[bond.j])
    return topology


def write_openmm(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    outdir: Path,
) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    system = to_openmm_system(ir, decompose=False)
    xml_path = outdir / "system.xml"
    xml_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")

    topology = to_openmm_topology(ir)
    positions = [openmm.Vec3(*p) for p in np.asarray(positions_nm, dtype=np.float64)]
    pdb_path = outdir / "molecule.pdb"
    with open(pdb_path, "w", encoding="utf-8") as fh:
        openmm.app.PDBFile.writeFile(topology, positions * unit.nanometer, fh)
    return {"system": xml_path, "pdb": pdb_path}
