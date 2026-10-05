"""T11：units、CartesianHessian、质量检查与谐振频率的必测用例。

不对称度阈值按 docs/decisions/0007：>1e-3 报错、>1e-6 警告
（计划测试用例原文"1e-4 → 报错"与 §2.3 表矛盾，用户确认以表为准）。
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.constants as sc
from scipy.linalg import null_space
from rdkit import Chem

from ffcompiler.hessian.base import (
    CartesianHessian,
    HessianQualityError,
    _tr_basis,
    average_masses,
    check_hessian,
    harmonic_frequencies,
    symmetrized,
)
from ffcompiler.units import BOHR_NM, C_CM_S, HARTREE_KJ_MOL

K_BOND = 5e5  # kJ/mol/nm²，双原子用例的键力常数


def _diatomic_hessian_object(
    *, k: float = K_BOND, gradient: np.ndarray | None = None
) -> CartesianHessian:
    """两原子相距 0.1 nm、沿 x 的解析键 Hessian（§2.2 用例）。"""
    z = np.array([1, 17])
    pos = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    h = np.zeros((6, 6))
    h[0, 0] = k
    h[3, 3] = k
    h[0, 3] = -k
    h[3, 0] = -k
    return CartesianHessian(
        atomic_numbers=z,
        positions_nm=pos,
        hessian=h,
        gradient=gradient,
        energy_kj=None,
        provider="test",
        method="analytic",
    )


def test_units_match_scipy_constants() -> None:
    hartree_kj = sc.value("Hartree energy") * sc.value("Avogadro constant") / 1000.0
    assert abs(HARTREE_KJ_MOL - hartree_kj) / hartree_kj <= 1e-9
    bohr_nm = sc.value("Bohr radius") / 1e-9
    assert abs(BOHR_NM - bohr_nm) / bohr_nm <= 1e-9
    assert abs(C_CM_S - sc.c * 100.0) / (sc.c * 100.0) <= 1e-9


def test_diatomic_frequency_analytic() -> None:
    h = _diatomic_hessian_object()
    masses = np.array([1.0, 19.0])
    freqs = harmonic_frequencies(h.hessian, h.positions_nm, masses)
    assert freqs.shape == (1,)  # 线性分子 3N−5 = 1
    mu = 1.0 * 19.0 / (1.0 + 19.0)
    expected = math.sqrt(K_BOND / mu * 1e24) / (2.0 * math.pi * C_CM_S)
    assert abs(freqs[0] - expected) <= 1e-6

    diag = check_hessian(h, masses)
    assert diag.frequencies_cm.shape == (1,)
    assert diag.frequencies_cm[0] > 0.0
    assert diag.warnings == ()


def test_diatomic_mass_fallback_rdkit() -> None:
    h = _diatomic_hessian_object()
    diag = check_hessian(h)  # masses=None → RDKit 平均原子量
    table = Chem.GetPeriodicTable()
    m1 = table.GetAtomicWeight(1)
    m2 = table.GetAtomicWeight(17)
    mu = m1 * m2 / (m1 + m2)
    expected = math.sqrt(K_BOND / mu * 1e24) / (2.0 * math.pi * C_CM_S)
    assert abs(diag.frequencies_cm[0] - expected) <= 1e-6
    assert np.allclose(average_masses(np.array([1, 17])), [m1, m2])


def test_translation_violation_raises() -> None:
    h = _diatomic_hessian_object()
    bad = h.hessian.copy()
    bad[0, 0] += 1e-2 * np.abs(h.hessian).max()
    broken = CartesianHessian(
        atomic_numbers=h.atomic_numbers,
        positions_nm=h.positions_nm,
        hessian=bad,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="analytic",
    )
    with pytest.raises(HessianQualityError, match="translation"):
        check_hessian(broken, np.array([1.0, 19.0]))


def test_asymmetry_thresholds_per_adr0007() -> None:
    # >1e-3 → 报错
    h = _diatomic_hessian_object()
    bad = h.hessian.copy()
    bad[0, 3] += 1e-2 * np.abs(h.hessian).max()
    broken = CartesianHessian(
        atomic_numbers=h.atomic_numbers,
        positions_nm=h.positions_nm,
        hessian=bad,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="analytic",
    )
    with pytest.raises(HessianQualityError, match="asymmetry"):
        check_hessian(broken, np.array([1.0, 19.0]))

    # 1e-4 → 警告（计划原文写报错，决策 0007 改为警告）
    warn = h.hessian.copy()
    warn[0, 3] += 0.5e-4 * np.abs(h.hessian).max()
    warn[3, 0] -= 0.5e-4 * np.abs(h.hessian).max()
    warned = CartesianHessian(
        atomic_numbers=h.atomic_numbers,
        positions_nm=h.positions_nm,
        hessian=warn,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="analytic",
    )
    diag = check_hessian(warned, np.array([1.0, 19.0]))
    # 扰动后 max|H| = k+δ/2，故比值略低于 1e-4
    assert diag.asymmetry_rel == pytest.approx(1e-4, rel=1e-3)
    assert any("asymmetry" in w for w in diag.warnings)
    # 反对称扰动经 ½(H+Hᵀ) 后消失：平移残差为 0，警告只有不对称一条
    assert diag.translation_residual_rel == 0.0
    assert len(diag.warnings) == 1

    # 1e-8 → 通过，无警告；symmetrized 后严格对称
    small = h.hessian.copy()
    small[0, 3] += 1e-8 * np.abs(h.hessian).max()
    minor = CartesianHessian(
        atomic_numbers=h.atomic_numbers,
        positions_nm=h.positions_nm,
        hessian=small,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="analytic",
    )
    diag = check_hessian(minor, np.array([1.0, 19.0]))
    assert diag.warnings == ()
    s = symmetrized(minor)
    assert np.array_equal(s.hessian, s.hessian.T)
    assert np.array_equal(s.atomic_numbers, minor.atomic_numbers)
    assert np.array_equal(s.positions_nm, minor.positions_nm)


def test_nonstationary_gradient() -> None:
    gradient = np.zeros(6)
    gradient[0] = 100.0  # kJ/mol/nm，远超 25
    h = _diatomic_hessian_object(gradient=gradient)
    with pytest.raises(HessianQualityError, match="stationary"):
        check_hessian(h, np.array([1.0, 19.0]))

    diag = check_hessian(
        h, np.array([1.0, 19.0]), allow_nonstationary=True
    )  # 只记录
    assert diag.max_gradient == pytest.approx(100.0)
    assert any("stationary" in w for w in diag.warnings)


def _water_with_vibrations(wavenumbers_cm: tuple[float, float, float]) -> CartesianHessian:
    """水分子几何上构造 H = M^{1/2}·Qc·diag(λ)·Qcᵀ·M^{1/2}。

    λ 由给定波数反推：振动特征值精确可控，用于虚频报错/警告用例。
    """
    masses = np.array([15.999, 1.008, 1.008])
    theta = math.radians(104.52)
    r_oh = 0.09572
    pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [r_oh, 0.0, 0.0],
            [r_oh * math.cos(theta), r_oh * math.sin(theta), 0.0],
        ]
    )
    q_tr = np.array(_tr_basis(pos, masses)).T            # (9,6)
    q_vib = null_space(q_tr.T)
    assert q_vib.shape == (9, 3)
    lam = np.array(
        [
            math.copysign((2.0 * math.pi * C_CM_S * abs(nu)) ** 2 / 1e24, nu)
            for nu in wavenumbers_cm
        ]
    )
    hw = q_vib @ np.diag(lam) @ q_vib.T
    sqrt_m3 = np.repeat(np.sqrt(masses), 3)
    h = sqrt_m3[:, None] * hw * sqrt_m3[None, :]
    return CartesianHessian(
        atomic_numbers=np.array([8, 1, 1]),
        positions_nm=pos,
        hessian=h,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="synthetic",
    )


def test_imaginary_frequency_error_and_warning() -> None:
    saddle = _water_with_vibrations((1600.0, 3800.0, -50.0))
    with pytest.raises(HessianQualityError, match="imaginary"):
        check_hessian(saddle, np.array([15.999, 1.008, 1.008]))

    shallow = _water_with_vibrations((1600.0, 3800.0, -5.0))
    diag = check_hessian(shallow, np.array([15.999, 1.008, 1.008]))
    assert diag.n_imaginary == 1
    assert diag.frequencies_cm[0] == pytest.approx(-5.0, abs=1e-6)
    assert any("imaginary" in w for w in diag.warnings)
