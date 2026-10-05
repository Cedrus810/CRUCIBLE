"""T19：DRIH 的必测用例。

oracle 为真实断言：Wilson B 以 Bohr 为长度单位构造后，全部 oracle case
与 Sobtop DRIH 一致（决策 0012，取代 0010）。SF6 诊断只打印不断言。
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.drih import drih
from ffcompiler.bonded.internal import (
    LinearAngleError,
    linear_angles,
    proper_quartets,
    wilson_b,
)
from ffcompiler.bonded.sobtop_ref import read_sobtop_params
from ffcompiler.hessian.base import CartesianHessian, HessianQualityError
from ffcompiler.hessian.gaussian import read_fchk

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from sobtop_oracle import (  # noqa: E402
    GATE_CASES,
    REFS,
    load_case,
    oracle_mismatches,
    quartets_from,
)

_quartets_from = quartets_from


def test_oracle_full_matrix() -> None:
    """用例 1：除 sf6_fchk 外全部 case 与 param_m4.txt 一致（§2.8 容差，含平衡值）。"""
    failures = {}
    for case in GATE_CASES:
        h, bonds, angles = load_case(case)
        terms = drih(h, bonds, angles, dihedrals=quartets_from(bonds))
        bad = oracle_mismatches(terms, read_sobtop_params(REFS / case / "param_m4.txt"))
        if bad:
            failures[case] = bad[:5]
    assert not failures, failures


def test_linear_angle_in_input_raises() -> None:
    """§2.4：所有方法遇到线性键角都抛 LinearAngleError，不得静默丢弃。"""
    h, bonds, angles = load_case("sf6_fchk")
    with pytest.raises(LinearAngleError):
        drih(h, bonds, angles, dihedrals=[])


def test_dihedrals_through_linear_angle_excluded() -> None:
    """library 策略的调用方式：angle 列表已去掉线性角，经过线性角的四元组由
    drih 按几何自行剔除（任一端三元组的顺序都要识别）。"""
    # 乙腈骨架 H3C-C≡N：N(0)-C(1)-C(2)，C(2) 上 3 个 H
    pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1157, 0.0, 0.0],
            [0.2615, 0.0, 0.0],
            [0.2980, 0.1028, 0.0],
            [0.2980, -0.0514, 0.0890],
            [0.2980, -0.0514, -0.0890],
        ]
    )
    bonds = [(0, 1), (1, 2), (2, 3), (2, 4), (2, 5)]
    all_angles = [(0, 1, 2), (1, 2, 3), (1, 2, 4), (1, 2, 5), (3, 2, 4), (3, 2, 5), (4, 2, 5)]
    lin = linear_angles(pos, all_angles)
    assert lin
    angles = [a for a in all_angles if a not in lin]
    quartets = quartets_from(bonds)
    # (3,2,1,0) 这类四元组的 (2,1,0) 端三元组与线性角 (0,1,2) 顺序相反
    assert any(q[1:] == (2, 1, 0) or q[:3] == (0, 1, 2) for q in quartets)
    coords = [("bond", b) for b in bonds] + [("angle", a) for a in angles]
    b_matrix = wilson_b(pos, coords)
    hessian = sum(1000.0 * np.outer(b, b) for b in b_matrix)
    h = CartesianHessian(
        atomic_numbers=np.array([7, 6, 6, 1, 1, 1]), positions_nm=pos, hessian=hessian,
        gradient=np.zeros(18), energy_kj=None, provider="test", method="synthetic",
    )
    terms = drih(h, bonds, angles, dihedrals=quartets, allow_nonstationary=True)
    assert len(terms) == len(bonds) + len(angles)


def test_sf6_diagnostic(capfd: pytest.CaptureFixture) -> None:
    """SF6（180° 键角）：bond 与 Sobtop 一致；90° 角差异打印不断言。"""
    h, bonds, angles = load_case("sf6_fchk")
    params = read_sobtop_params(REFS / "sf6_fchk" / "param_m4.txt")
    # linear_angle_policy 的等效做法：先去掉 180° 键角
    linear = set(linear_angles(h.positions_nm, angles))
    terms = drih(h, bonds, [a for a in angles if a not in linear], dihedrals=[])
    for t in terms:
        if len(t.atoms) == 2:
            _, k_ref = params.bonds[tuple(sorted(t.atoms))]
            ours = t.k / 100.0
            assert abs(ours - k_ref) <= 1e-4 * abs(k_ref) + 5e-3
    for kind, atoms, k in (
        ("a", t.atoms, t.k) for t in terms if len(t.atoms) == 3
    ):
        i, j, k_at = atoms
        _, k_ref = params.angles[(min(i, k_at), j, max(i, k_at))]
        print(f"SF6 90°角 {atoms}: ours={k:.3f} sobtop={k_ref:.3f} 差={k - k_ref:.3f}")


def _synthetic_hessian(pos, coords, ks):
    b_matrix = wilson_b(pos, coords)
    return sum(k * np.outer(b, b) for k, b in zip(ks, b_matrix))


def _synthetic_h(pos, coords, ks) -> CartesianHessian:
    n3 = 3 * pos.shape[0]
    return CartesianHessian(
        atomic_numbers=np.ones(pos.shape[0], dtype=np.int64),
        positions_nm=pos,
        hessian=_synthetic_hessian(pos, coords, ks),
        gradient=np.zeros(n3),
        energy_kj=None, provider="test", method="synthetic",
    )


def test_exact_recovery_water() -> None:
    pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.09572, 0.0, 0.0],
            [0.09572 * math.cos(math.radians(104.52)), 0.09572 * math.sin(math.radians(104.52)), 0.0],
        ]
    )
    coords = [("bond", (0, 1)), ("bond", (0, 2)), ("angle", (1, 0, 2))]
    ks = [4.0e5, 5.0e5, 300.0]
    h = _synthetic_h(pos, coords, ks)
    terms = drih(h, [(0, 1), (0, 2)], [(1, 0, 2)], dihedrals=[])
    for t, k_ref in zip(terms, ks):
        assert t.k == pytest.approx(k_ref, rel=1e-8)


def test_exact_recovery_h2o2() -> None:
    # 非平面 H2O2：二面角约 110°
    phi = math.radians(110.0)
    r_oo, r_oh, theta = 0.145, 0.097, math.radians(100.0)
    pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [r_oo, 0.0, 0.0],
            [r_oh * math.cos(theta), r_oh * math.sin(theta), 0.0],
            [
                r_oo - r_oh * math.cos(theta),
                r_oh * math.sin(theta) * math.cos(phi),
                r_oh * math.sin(theta) * math.sin(phi),
            ],
        ]
    )
    coords = [
        ("bond", (0, 1)), ("bond", (0, 2)), ("bond", (1, 3)),
        ("angle", (2, 0, 1)), ("angle", (0, 1, 3)), ("dihedral", (2, 0, 1, 3)),
    ]
    ks = [3.5e5, 4.5e5, 4.5e5, 350.0, 350.0, 120.0]
    h = _synthetic_h(pos, coords, ks)
    terms = drih(
        h, [(0, 1), (0, 2), (1, 3)], [(2, 0, 1), (0, 1, 3)],
        dihedrals=[(2, 0, 1, 3)],
    )
    # 输出只有 bond+angle（dihedral 只进坐标集）
    assert len(terms) == 5
    for t, k_ref in zip(terms, ks[:5]):
        assert t.k == pytest.approx(k_ref, rel=1e-8)


def test_redundant_ethanol_recorded(capsys: pytest.CaptureFixture) -> None:
    """sp³ 冗余下不能精确恢复：打印最大相对偏差，不断言（DRIH 固有性质）。"""
    h, bonds, angles = load_case("ethanol_gfn2")
    coords = [("bond", b) for b in bonds] + [("angle", a) for a in angles]
    from ffcompiler.bonded.internal import wilson_b as wb

    b_matrix = wb(h.positions_nm, coords)
    ks = np.full(len(coords), 1000.0)
    hessian = sum(k * np.outer(b, b) for k, b in zip(ks, b_matrix))
    h2 = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=h.positions_nm,
        hessian=hessian, gradient=np.zeros(hessian.shape[0]),
        energy_kj=None, provider="test", method="synthetic",
    )
    terms = drih(h2, bonds, angles, dihedrals=_quartets_from(bonds))
    worst = max(abs(t.k - 1000.0) / 1000.0 for t in terms)
    print(f"ethanol 冗余恢复最大相对偏差 = {worst:.3e}（>0 即为 DRIH 固有混合）")
    assert worst > 0.0


def test_nonstationary_rejected() -> None:
    """用例 5：fchk 的梯度人为放大到 max|g| = 100 kJ/mol/nm → 默认报错（Review Focus 3）。"""
    h, bonds, angles = load_case("methyl_benzoate_fchk")
    assert h.gradient is not None
    scale = 100.0 / float(np.max(np.abs(h.gradient)))
    h_big = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=h.positions_nm,
        hessian=h.hessian, gradient=h.gradient * scale,
        energy_kj=None, provider="test", method=h.method,
    )
    quartets = quartets_from(bonds)
    with pytest.raises(HessianQualityError, match="stationary"):
        drih(h_big, bonds, angles, dihedrals=quartets)
    drih(h_big, bonds, angles, dihedrals=quartets, allow_nonstationary=True)  # 不抛错
