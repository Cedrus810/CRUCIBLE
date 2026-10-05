"""T2：IR 校验的必测用例——8 个校验项各一个反例，外加一个 3 项同时违反。"""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

from ffcompiler.ir import (
    Angle,
    Atom,
    Bond,
    ForceFieldIR,
    ImproperTorsion,
    ProperTorsion,
)
from ffcompiler.ir.topology import expected_exclusions_and_pairs
from ffcompiler.ir.validation import IRValidationError, validate


def _atoms(n: int, **overrides) -> tuple[Atom, ...]:
    out = []
    for i in range(n):
        kw = dict(
            index=i, name=f"A{i}", element="C", atomic_number=6, mass=12.0 + 0.01 * i,
            atom_type=f"T{i}", formal_charge=0, partial_charge=0.0,
            sigma=0.3 + 0.01 * i, epsilon=0.4 + 0.01 * i,
        )
        kw.update(overrides.get(i, {}))
        out.append(Atom(**kw))
    return tuple(out)


def _chain_ir() -> ForceFieldIR:
    """4 原子链 C0-C1-C2-C3：2 个 angle、1 个 proper，bookkeeping 完整、可过 validate。"""
    atoms = _atoms(4)
    ir = ForceFieldIR(
        name="chain", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.15, k=3.0e5, source="t"),
            Bond(i=1, j=2, r0=0.15, k=3.0e5, source="t"),
            Bond(i=2, j=3, r0=0.15, k=3.0e5, source="t"),
        ),
        angles=(
            Angle(i=0, j=1, k_atom=2, theta0=1.9, force_constant=500.0, source="t"),
            Angle(i=1, j=2, k_atom=3, theta0=1.9, force_constant=500.0, source="t"),
        ),
        proper_torsions=(
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=3, phase=0.0,
                          force_constant=1.0, source="t"),
        ),
        impropers=(), exclusions=(), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    exclusions, pairs14 = expected_exclusions_and_pairs(ir)
    return replace(ir, exclusions=exclusions, pairs14=pairs14)


def _improper_ir() -> ForceFieldIR:
    """C0 中心 + O1/H2/H3：improper (1,2,0,3)，第三个原子 C0 与 i、j、l 都成键。"""
    atoms = _atoms(4)
    ir = ForceFieldIR(
        name="imp", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.12, k=4.0e5, source="t"),
            Bond(i=0, j=2, r0=0.11, k=3.0e5, source="t"),
            Bond(i=0, j=3, r0=0.11, k=3.0e5, source="t"),
        ),
        angles=(
            Angle(i=1, j=0, k_atom=2, theta0=2.0, force_constant=400.0, source="t"),
            Angle(i=1, j=0, k_atom=3, theta0=2.0, force_constant=400.0, source="t"),
            Angle(i=2, j=0, k_atom=3, theta0=2.1, force_constant=300.0, source="t"),
        ),
        proper_torsions=(),
        impropers=(
            ImproperTorsion(i=1, j=2, k=0, l=3, periodicity=2, phase=math.pi,
                            force_constant=20.0, source="t"),
        ),
        exclusions=(), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    exclusions, pairs14 = expected_exclusions_and_pairs(ir)
    return replace(ir, exclusions=exclusions, pairs14=pairs14)


def test_valid_base_irs_pass() -> None:
    validate(_chain_ir())
    validate(_improper_ir())


@pytest.mark.parametrize("field,value", [
    ("partial_charge", float("nan")), ("sigma", float("nan")),
    ("epsilon", -0.5), ("epsilon", float("inf")), ("mass", -1.0),
    ("mass", 0.0), ("sigma", -0.1),
])
def test_invalid_atom_numbers_rejected(field, value):
    ir = _chain_ir()
    bad = replace(ir, atoms=(replace(ir.atoms[0], **{field: value}),) + ir.atoms[1:])
    with pytest.raises(IRValidationError, match=field):
        validate(bad)


def test_zero_lj_remains_supported():
    ir = _chain_ir()
    validate(replace(ir, atoms=(replace(ir.atoms[0], sigma=0.0, epsilon=0.0),) + ir.atoms[1:]))


@pytest.mark.parametrize("field", ["k", "r0", "uncertainty"])
def test_nonfinite_bond_numbers_rejected(field):
    ir = _chain_ir()
    with pytest.raises(IRValidationError, match=field):
        validate(replace(ir, bonds=(replace(ir.bonds[0], **{field: float("inf")}),) + ir.bonds[1:]))


@pytest.mark.parametrize("field", ["phase", "force_constant"])
def test_nonfinite_torsion_numbers_rejected(field):
    ir = _chain_ir()
    with pytest.raises(IRValidationError, match=field):
        validate(replace(ir, proper_torsions=(replace(ir.proper_torsions[0], **{field: float("nan")}),)))


def test_1_index_out_of_order_or_out_of_range() -> None:
    ir = _chain_ir()
    bad_atoms = list(ir.atoms)
    bad_atoms[2] = replace(bad_atoms[2], index=9)
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, atoms=tuple(bad_atoms)))
    assert any("index" in s for s in err.value.issues)

    bad_bond = Bond(i=0, j=9, r0=0.15, k=3.0e5, source="t")
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, bonds=ir.bonds + (bad_bond,)))
    assert any("out-of-range" in s for s in err.value.issues)


def test_2_missing_connectivity() -> None:
    ir = _chain_ir()
    ghost = Angle(i=0, j=1, k_atom=3, theta0=1.9, force_constant=500.0, source="t")
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, angles=ir.angles + (ghost,)))
    assert any("angle" in s and "lacks bond (1,3)" in s for s in err.value.issues)

    imp = _improper_ir()
    broken = replace(imp, bonds=imp.bonds[:2])  # 删掉 (0,3)
    with pytest.raises(IRValidationError) as err:
        validate(broken)
    # 决策 0002：中心原子可以是任何位置，但必须恰有一个原子与其余三个成键
    assert any("improper (1,2,0,3)" in s and "exactly one atom" in s
               for s in err.value.issues)


def test_3_duplicates() -> None:
    ir = _chain_ir()
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, bonds=ir.bonds + (Bond(i=2, j=3, r0=0.16, k=2.9e5, source="t"),)))
    assert any("duplicate bond" in s for s in err.value.issues)

    with pytest.raises(IRValidationError) as err:
        validate(replace(
            ir,
            proper_torsions=ir.proper_torsions + (
                ProperTorsion(i=0, j=1, k=2, l=3, periodicity=3, phase=0.5,
                              force_constant=2.0, source="t"),
            ),
        ))
    assert any("duplicate proper" in s and "periodicity 3" in s for s in err.value.issues)

    # 反向四元组视为同一个：原 proper 是 (0,1,2,3) n=3，反向 (3,2,1,0) n=3 构成重复
    with pytest.raises(IRValidationError) as err:
        validate(replace(
            ir,
            proper_torsions=ir.proper_torsions + (
                ProperTorsion(i=3, j=2, k=1, l=0, periodicity=3, phase=0.0,
                              force_constant=2.0, source="t"),
            ),
        ))
    assert any("duplicate proper" in s and "periodicity 3" in s for s in err.value.issues)


def test_4_parameter_ranges() -> None:
    ir = _chain_ir()
    bad_bonds = list(ir.bonds)
    bad_bonds[0] = replace(bad_bonds[0], k=-1.0)
    bad_bonds[1] = replace(bad_bonds[1], r0=0.0)
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, bonds=tuple(bad_bonds)))
    assert any("bond (0,1) has k=-1.0" in s for s in err.value.issues)
    assert any("r0=0.0" in s for s in err.value.issues)

    bad_angles = list(ir.angles)
    bad_angles[0] = replace(bad_angles[0], force_constant=0.0)
    bad_angles[1] = replace(bad_angles[1], theta0=math.pi + 0.1)
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, angles=tuple(bad_angles)))
    assert any("force_constant=0.0" in s for s in err.value.issues)
    assert any("theta0" in s for s in err.value.issues)

    bad_torsion = ProperTorsion(i=0, j=1, k=2, l=3, periodicity=7, phase=0.0,
                                force_constant=1.0, source="t")
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, proper_torsions=(bad_torsion,)))
    assert any("periodicity=7" in s for s in err.value.issues)


def test_5_exclusions_and_pairs14_must_match_expected() -> None:
    ir = _chain_ir()
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, exclusions=()))
    assert any("exclusions" in s for s in err.value.issues)

    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, pairs14=((0, 2),)))
    assert any("pairs14" in s for s in err.value.issues)

    # 同一个对不能同时落在两处
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, exclusions=ir.exclusions + ((0, 3),)))
    assert any("both exclusions and pairs14" in s for s in err.value.issues)


def test_6_charge_conservation() -> None:
    ir = _chain_ir()
    atoms = list(ir.atoms)
    atoms[0] = replace(atoms[0], partial_charge=0.5)
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, atoms=tuple(atoms)))
    assert any("charge not conserved" in s for s in err.value.issues)


def test_7_same_atom_type_same_params() -> None:
    ir = _chain_ir()
    atoms = list(ir.atoms)
    atoms[1] = replace(atoms[1], atom_type="T0")  # 与原子 0 同类型但 sigma/mass 不同
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, atoms=tuple(atoms)))
    assert any("atom_type 'T0' has inconsistent" in s for s in err.value.issues)


def test_8_family_and_scales() -> None:
    ir = _chain_ir()
    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, forcefield_family="gaff"))
    assert any("forcefield_family" in s for s in err.value.issues)

    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, combination_rule="geometric"))
    assert any("combination_rule" in s for s in err.value.issues)

    with pytest.raises(IRValidationError) as err:
        validate(replace(ir, coulomb_14_scale=1.5))
    assert any("coulomb_14_scale" in s for s in err.value.issues)


def test_multiple_issues_all_reported_at_once() -> None:
    ir = _chain_ir()
    broken = replace(
        ir,
        bonds=ir.bonds + (Bond(i=2, j=3, r0=0.15, k=3.0e5, source="t"),),  # 重复 bond
        proper_torsions=(
            ProperTorsion(i=0, j=1, k=2, l=3, periodicity=9, phase=0.0,
                          force_constant=1.0, source="t"),  # periodicity 非法
        ),
        forcefield_family="gaff",  # family 错
    )
    with pytest.raises(IRValidationError) as err:
        validate(broken)
    issues = err.value.issues
    assert any("duplicate bond" in s for s in issues)
    assert any("periodicity=9" in s for s in issues)
    assert any("forcefield_family" in s for s in issues)
    assert len(issues) >= 3
