"""Fourier 拟合与 IR 周期项换算（MVP-2 计划 §2.4，T24）。

y ≈ c0 + Σ_{n=1}^{n_max} (a_n cos nφ + b_n sin nφ)，均匀权重线性最小二乘；
a_n cos nφ + b_n sin nφ = k_n cos(nφ − δ_n)，k_n = √(a_n²+b_n²)，δ_n = atan2(b_n, a_n)。
IR 形式 Σ k_n [1 + cos(nφ − δ_n)] 比上式多一个常数 Σ k_n，并入 c0。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from ..ir.schema import ProperTorsion


@dataclass(frozen=True)
class FourierFit:
    k: tuple[float, ...]        # n = 1..n_max，kJ/mol，≥ 0
    phase: tuple[float, ...]    # rad
    c0: float                   # kJ/mol
    rmse: float
    max_abs_residual: float


def _design(phi: np.ndarray, n_max: int) -> np.ndarray:
    cols = [np.ones_like(phi)]
    for n in range(1, n_max + 1):
        cols.append(np.cos(n * phi))
        cols.append(np.sin(n * phi))
    return np.stack(cols, axis=1)


def fit_fourier(phi_rad: np.ndarray, energy: np.ndarray, n_max: int = 6) -> FourierFit:
    phi = np.asarray(phi_rad, dtype=np.float64).ravel()
    y = np.asarray(energy, dtype=np.float64).ravel()
    if phi.shape != y.shape:
        raise ValueError(f"phi {phi.shape} 与 energy {y.shape} 长度不一致")
    n_par = 2 * n_max + 1
    if phi.size < n_par:
        raise ValueError(f"至少需要 {n_par} 个点才能拟合 n_max={n_max}，只有 {phi.size} 个")
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(y))):
        raise ValueError("phi 或 energy 含非有限值")
    a = _design(phi, n_max)
    coef, *_ = np.linalg.lstsq(a, y, rcond=None)
    resid = y - a @ coef
    ks, phases = [], []
    for n in range(1, n_max + 1):
        an, bn = coef[2 * n - 1], coef[2 * n]
        ks.append(float(np.hypot(an, bn)))
        phases.append(float(np.arctan2(bn, an)))
    return FourierFit(
        k=tuple(ks),
        phase=tuple(phases),
        c0=float(coef[0]),
        rmse=float(np.sqrt(np.mean(resid**2))),
        max_abs_residual=float(np.max(np.abs(resid))),
    )


def fourier_energy(fit: FourierFit, phi_rad) -> np.ndarray:
    """拟合曲线（含 c0），kJ/mol。"""
    phi = np.asarray(phi_rad, dtype=np.float64)
    e = np.full_like(phi, fit.c0)
    for n, (k, d) in enumerate(zip(fit.k, fit.phase), start=1):
        e = e + k * np.cos(n * phi - d)
    return e


def to_ir_terms(
    fit: FourierFit,
    quartet: tuple[int, int, int, int],
    *,
    source: str,
    meta: Mapping[str, str | int | float | bool],
) -> list[ProperTorsion]:
    """IR 周期项，n = 1..n_max 全写（不设"太小就删"的阈值，计划 T24 裁定）。"""
    i, j, k, l = (int(a) for a in quartet)
    return [
        ProperTorsion(i=i, j=j, k=k, l=l, periodicity=n, phase=d, force_constant=kn,
                      source=source, meta=dict(meta))
        for n, (kn, d) in enumerate(zip(fit.k, fit.phase), start=1)
    ]
