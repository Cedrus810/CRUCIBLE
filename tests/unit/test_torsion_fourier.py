"""T24：Fourier 拟合的必测用例（MVP-2 计划 §2.4）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffcompiler.torsion.fourier import fit_fourier, fourier_energy, to_ir_terms

K_TRUE = (1.3, 0.0, 4.2, 0.7, 0.0, 0.25)
PHASE_TRUE = (0.4, 0.0, math.pi, -1.1, 0.0, 2.5)
C0_TRUE = -12.0


def _synthetic(phi: np.ndarray) -> np.ndarray:
    e = np.full_like(phi, C0_TRUE)
    for n, (k, d) in enumerate(zip(K_TRUE, PHASE_TRUE), start=1):
        e += k * np.cos(n * phi - d)
    return e


def _wrap(x: float) -> float:
    return math.atan2(math.sin(x), math.cos(x))


def test_recovers_known_coefficients_on_nonuniform_grid() -> None:
    rng = np.random.default_rng(20261004)
    phi = np.sort(rng.uniform(-math.pi, math.pi, 36))
    fit = fit_fourier(phi, _synthetic(phi))
    assert len(fit.k) == 6 and len(fit.phase) == 6
    for n in range(6):
        assert fit.k[n] == pytest.approx(K_TRUE[n], abs=1e-9)
        if K_TRUE[n] > 0:
            assert abs(_wrap(fit.phase[n] - PHASE_TRUE[n])) <= 1e-9
    assert fit.rmse <= 1e-9
    assert fit.max_abs_residual <= 1e-9


def test_ir_terms_match_fourier_energy_up_to_constant() -> None:
    rng = np.random.default_rng(1)
    phi = np.linspace(-math.pi, math.pi, 36, endpoint=False)
    fit = fit_fourier(phi, _synthetic(phi))
    terms = to_ir_terms(fit, (0, 1, 2, 3), source="mlp", meta={"model": "test"})
    assert [t.periodicity for t in terms] == [1, 2, 3, 4, 5, 6]  # 六项全写，不按大小删
    assert all(t.source == "mlp" and t.meta["model"] == "test" for t in terms)
    assert all((t.i, t.j, t.k, t.l) == (0, 1, 2, 3) for t in terms)
    probe = rng.uniform(-math.pi, math.pi, 100)
    e_ir = np.array([
        sum(t.force_constant * (1.0 + math.cos(t.periodicity * p - t.phase)) for t in terms)
        for p in probe
    ])
    diff = e_ir - fourier_energy(fit, probe)
    assert float(np.std(diff)) <= 1e-9


def test_too_few_points_raises() -> None:
    phi = np.linspace(-math.pi, math.pi, 12, endpoint=False)
    with pytest.raises(ValueError):
        fit_fourier(phi, np.zeros_like(phi))
