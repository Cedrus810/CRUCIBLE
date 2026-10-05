"""T6：OpenMM 导出器的必测用例。"""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import unit

from ffcompiler.exporters.openmm import FORCE_GROUPS, to_openmm_system, write_openmm
from ffcompiler.ir import Atom, Bond, ForceFieldIR, ImproperTorsion
from ffcompiler.ir.energy import evaluate, forces_fd
from ffcompiler.ir.validation import validate
from tests.unit.test_energy import _chain_ir
from tests.unit.test_ir_schema import make_4atom_ir
from ffcompiler.validation.engines import reference_context

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
    "acetonitrile",
)


def _fixture_ir(name: str):
    from ffcompiler.importers import from_amber
    mdir = FIXTURES / name
    return from_amber(mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf")


def _eval_system(system: openmm.System, positions: np.ndarray, *, forces: bool = False):
    context = reference_context(system)
    context.setPositions([openmm.Vec3(*p) for p in positions] * unit.nanometer)
    if forces:
        state = context.getState(getEnergy=True, getForces=True)
        e = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        f = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        del context
        return e, f
    del context
    return context, None  # pragma: no cover


def _group_energies(system: openmm.System, positions: np.ndarray) -> tuple[dict[str, float], float]:
    context = reference_context(system)
    context.setPositions([openmm.Vec3(*p) for p in positions] * unit.nanometer)
    out = {}
    for name, group in FORCE_GROUPS.items():
        state = context.getState(getEnergy=True, groups=1 << group)
        out[name] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    del context
    return out, total


TOL_MEM_E = lambda e: 1e-6 + 1e-9 * abs(e)


@pytest.mark.parametrize("name", NAMES)
def test_decomposed_groups_match_evaluate(name: str) -> None:
    ir, positions = _fixture_ir(name)
    system = to_openmm_system(ir, decompose=True)
    groups, _total = _group_energies(system, positions)
    ref = evaluate(ir, positions)
    for key in FORCE_GROUPS:
        assert abs(groups[key] - ref[key]) <= TOL_MEM_E(ref[key]), (
            f"{name}: {key} OpenMM {groups[key]} vs IR {ref[key]}"
        )


def test_decomposed_groups_on_hand_built_irs() -> None:
    from ffcompiler.ir import ProperTorsion

    rng = np.random.default_rng(20261004)
    # 同一四元组 n=1,2,3 三项 proper（Review Focus 4）
    quartet_ir = _chain_ir().replace_terms(
        proper_torsions=tuple(
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=n, phase=0.1 * n,
                          force_constant=1.0 - 0.1 * n, source="t")
            for n in (1, 2, 3)
        )
    )
    for ir in (quartet_ir, _chain_ir()):
        positions = rng.normal(0.0, 0.15, size=(len(ir.atoms), 3))
        system = to_openmm_system(ir, decompose=True)
        groups, _ = _group_energies(system, positions)
        ref = evaluate(ir, positions)
        for key in FORCE_GROUPS:
            assert abs(groups[key] - ref[key]) <= TOL_MEM_E(max(ref[key], 1.0)), key


def test_production_mode_total_matches_decomposed() -> None:
    rng = np.random.default_rng(20261004)
    # 现实构型（骨架 + 0.02 nm 扰动）：纯高斯随机构型会出现近共线键角，
    # arccos 类求值在 θ→π 处的 ~1e-8 rad 误差会放大到能量 ~1e-7 量级
    base = np.array([
        [0.0, 0.0, 0.0],
        [0.15, 0.0, 0.0],
        [0.15 + 0.15 * math.cos(2.0), 0.15 * math.sin(2.0), 0.0],
        [0.30 + 0.15 * math.cos(2.0) + 0.15 * math.cos(4.0),
         0.15 * math.sin(2.0) + 0.15 * math.sin(4.0), 0.05],
    ])
    cases = [
        (_chain_ir(), base + rng.normal(0.0, 0.02, size=(4, 3))),
        _fixture_ir("methyl_benzoate"),
        _fixture_ir("acetate"),
    ]
    for ir, positions in cases:
        e_prod = _group_energies(to_openmm_system(ir, decompose=False), positions)[1]
        e_dec = _group_energies(to_openmm_system(ir, decompose=True), positions)[1]
        ref = evaluate(ir, positions)["total"]
        # e_prod vs e_dec：不同力内核（Nonbonded vs Custom*）求和次序不同，
        # 双精度实测地板 ~3e-9 相对。e_prod vs ref：直接用计划 §2.4 的
        # "内存内"容差（1e-6 + 1e-9|E|）。
        assert abs(e_prod - e_dec) <= 1e-8 * max(1.0, abs(e_dec))
        assert abs(e_prod - ref) <= 1e-6 + 1e-9 * abs(ref)


def _improper_order_ir(imp: ImproperTorsion) -> ForceFieldIR:
    atoms = tuple(
        Atom(index=i, name=f"A{i}", element="C", atomic_number=6, mass=12.0,
             atom_type=f"T{i}", formal_charge=0, partial_charge=0.0,
             sigma=0.0, epsilon=0.0)
        for i in range(4)
    )
    ir = ForceFieldIR(
        name="imp_order", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.15, k=3.0e5, source="t"),
            Bond(i=0, j=2, r0=0.15, k=3.0e5, source="t"),
            Bond(i=0, j=3, r0=0.15, k=3.0e5, source="t"),
        ),
        angles=(), proper_torsions=(), impropers=(imp,),
        exclusions=((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)), pairs14=(),
        total_charge=0, forcefield_family="amber",
        combination_rule="lorentz-berthelot", coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )
    validate(ir)
    return ir


def test_improper_order_not_permuted() -> None:
    s3 = math.sqrt(3) / 2
    positions = np.array([
        [-0.5, -s3, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.5, -0.433, -0.75],
    ])
    ir = _improper_order_ir(
        ImproperTorsion(i=0, j=1, k=2, l=3, periodicity=1, phase=0.5,
                        force_constant=2.0, source="t")
    )
    ir2 = _improper_order_ir(
        ImproperTorsion(i=3, j=1, k=2, l=0, periodicity=1, phase=0.5,
                        force_constant=2.0, source="t")
    )

    e_ir = evaluate(ir, positions)["improper"]
    e_ir_swapped = evaluate(ir2, positions)["improper"]
    assert abs(e_ir - e_ir_swapped) > 1e-6  # 两种顺序能量确实不同

    groups, _ = _group_energies(to_openmm_system(ir, decompose=True), positions)
    assert abs(groups["improper"] - e_ir) <= 1e-9  # OpenMM 按 IR 原顺序求值


def test_system_xml_roundtrip_energy_and_forces(tmp_path: Path) -> None:
    ir, positions = _fixture_ir("ethanol")
    rng = np.random.default_rng(20261004)
    configs = [
        positions,
        positions + rng.normal(0.0, 0.005, size=positions.shape),
        positions + rng.normal(0.0, 0.005, size=positions.shape),
    ]

    out = write_openmm(ir, positions, tmp_path)
    assert set(out) == {"system", "pdb"}
    system_xml = openmm.XmlSerializer.deserialize(out["system"].read_text())
    system_mem = to_openmm_system(ir, decompose=False)

    for pos in configs:
        e_xml, f_xml = _eval_system(system_xml, pos, forces=True)
        e_mem, f_mem = _eval_system(system_mem, pos, forces=True)
        assert abs(e_xml - e_mem) <= 1e-6 + 1e-9 * abs(e_mem)
        assert np.max(np.abs(f_xml - f_mem)) <= 1e-6


@pytest.mark.parametrize("name", NAMES)
def test_forces_match_central_difference(name: str) -> None:
    ir, positions = _fixture_ir(name)
    rng = np.random.default_rng(20261004)
    configs = [positions, positions + rng.normal(0.0, 0.005, size=positions.shape)]
    if name == "acetonitrile":
        # 乙腈的原始几何里 C-C≡N 角恰为 180°，二面角 φ 在线性中心上无定义，
        # 解析力与差分力都退化。比较前把两个构型都微微扰动（MVP-1 起由
        # linear_angle_policy 正式处理，见 MVP-1 计划 Review Focus 2）。
        configs = [pos + rng.normal(0.0, 0.005, size=pos.shape) for pos in configs]
    for pos in configs:
        system = to_openmm_system(ir, decompose=False)
        context = reference_context(system)
        context.setPositions([openmm.Vec3(*p) for p in pos] * unit.nanometer)
        f_omm = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        del context
        f_ref = forces_fd(ir, pos, step_nm=1e-6)
        assert np.max(np.abs(f_omm - f_ref)) <= 1e-4, name
