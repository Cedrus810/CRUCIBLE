"""CartesianHessian、质量检查与谐振频率（MVP-1 计划 T11）。

内部单位（计划 §2.1）：坐标 nm、Hessian kJ/mol/nm²、梯度 kJ/mol/nm、能量 kJ/mol。
检查阈值与频率公式按计划 §2.2、§2.3 锁定；不对称度报错阈值 1e-3 的取值见
docs/decisions/0007（计划测试用例与 §2.3 表矛盾，用户确认以表为准）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Mapping

import numpy as np
from rdkit import Chem

from ffcompiler.units import C_CM_S

# §2.3 检查阈值
_ASYMMETRY_ERROR = 1e-3
_ASYMMETRY_WARN = 1e-6
_TRANSLATION_ERROR = 1e-3
_TRANSLATION_WARN = 1e-5
_STATIONARY_MAX_GRADIENT = 25.0   # kJ/mol/nm（≈5×10⁻⁴ Ha/Bohr）
_IMAGINARY_ERROR_CM = -20.0       # 最小频率低于此值报错

# Gram–Schmidt 投影基向量的相对舍弃阈值：严格线性分子的绕轴转动向量恰为 0，
# 近线性分子（偏差 ≫1e-8）仍保留 6 个基向量
_PROJECTION_DROP_TOL = 1e-8


class HessianQualityError(ValueError):
    """Hessian 质量检查失败：不对称、破坏平移不变、非驻点或存在虚频。"""


@dataclass(frozen=True)
class CartesianHessian:
    """统一单位、统一原子顺序的 Cartesian Hessian（MVP-1 计划 §2.1）。"""

    atomic_numbers: np.ndarray          # (N,) int
    positions_nm: np.ndarray            # (N,3)
    hessian: np.ndarray                 # (3N,3N) kJ/mol/nm²
    gradient: np.ndarray | None         # (3N,) kJ/mol/nm
    energy_kj: float | None
    provider: str                       # "gaussian-fchk" | "xtb"
    method: str                         # 例如 "gfnff"、"gfn2"、fchk 里的方法名
    source_file: str | None = None
    source_sha256: str | None = None
    program_version: str | None = None
    # provider 的补充说明（如 xtb --ohess 的 "stationary_by"，计划 T13）
    meta: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        z = np.asarray(self.atomic_numbers)
        pos = np.asarray(self.positions_nm, dtype=np.float64)
        hess = np.asarray(self.hessian, dtype=np.float64)
        n = pos.shape[0]
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"positions_nm 形状应为 (N,3)，得到 {pos.shape}")
        if z.shape != (n,):
            raise ValueError(f"atomic_numbers 形状应为 ({n},)，得到 {z.shape}")
        if hess.shape != (3 * n, 3 * n):
            raise ValueError(f"hessian 形状应为 ({3*n},{3*n})，得到 {hess.shape}")
        if self.gradient is not None:
            g = np.asarray(self.gradient, dtype=np.float64)
            if g.shape != (3 * n,):
                raise ValueError(f"gradient 形状应为 ({3*n},)，得到 {g.shape}")
            object.__setattr__(self, "gradient", g)
        object.__setattr__(self, "atomic_numbers", z.astype(np.int64))
        object.__setattr__(self, "positions_nm", pos)
        object.__setattr__(self, "hessian", hess)


@dataclass(frozen=True)
class HessianDiagnostics:
    """check_hessian 的检查结果；只在全部报错条件未触发时返回。"""

    asymmetry_rel: float
    translation_residual_rel: float
    max_gradient: float | None
    frequencies_cm: np.ndarray
    n_imaginary: int
    warnings: tuple[str, ...]


def symmetrized(h: CartesianHessian) -> CartesianHessian:
    """返回 ½(H+Hᵀ) 的新 CartesianHessian；其余字段不变。

    IEEE 加法交换律保证 (H[i,j]+H[j,i])/2 与 (H[j,i]+H[i,j])/2 逐位相同，
    结果严格对称。
    """
    return replace(h, hessian=0.5 * (h.hessian + h.hessian.T))


def average_masses(atomic_numbers: np.ndarray) -> np.ndarray:
    """RDKit 周期表平均原子量（amu），计划 §2.2 的无 IR 兜底。"""
    table = Chem.GetPeriodicTable()
    return np.array(
        [table.GetAtomicWeight(int(z)) for z in atomic_numbers], dtype=np.float64
    )


def _tr_basis(positions_nm: np.ndarray, masses: np.ndarray) -> list[np.ndarray]:
    """质量加权平移+转动基向量（以质心为原点），Gram–Schmidt 正交化。

    线性分子绕分子轴的转动向量恰为 0（或近 0），按相对阈值舍弃，
    基向量个数即 6（非线性）或 5（线性）。
    """
    n = positions_nm.shape[0]
    com = masses @ positions_nm / masses.sum()
    sqrt_m = np.sqrt(masses)

    vectors: list[np.ndarray] = []
    for a in range(3):
        v = np.zeros((n, 3))
        v[:, a] = sqrt_m
        vectors.append(v.ravel())
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        disp = np.cross(e, positions_nm - com)
        vectors.append((sqrt_m[:, None] * disp).ravel())

    scale = float(np.linalg.norm(vectors[0]))
    basis: list[np.ndarray] = []
    for v in vectors:
        w = v.astype(np.float64).copy()
        for b in basis:
            w = w - (b @ w) * b
        norm = float(np.linalg.norm(w))
        if norm > _PROJECTION_DROP_TOL * scale:
            basis.append(w / norm)
    return basis


def harmonic_frequencies(
    hessian: np.ndarray, positions_nm: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    """谐振频率（cm⁻¹），升序；负值表示虚频（计划 §2.2）。

    质量加权 H̃ = M^{-1/2} H M^{-1/2}，投影掉以质心为原点的平移和转动
    （6 个基向量，线性分子 5 个，Gram–Schmidt 正交化），再对角化。
    特征值 λ（kJ mol⁻¹ nm⁻² amu⁻¹）→ ν̃ = sign(λ)·√(|λ|×10²⁴)/(2πc)。
    返回 3N−6（线性 3N−5）个频率。
    """
    h = np.asarray(hessian, dtype=np.float64)
    pos = np.asarray(positions_nm, dtype=np.float64)
    m = np.asarray(masses, dtype=np.float64)
    n = pos.shape[0]
    if m.shape != (n,):
        raise ValueError(f"masses 形状应为 ({n},)，得到 {m.shape}")
    if h.shape != (3 * n, 3 * n):
        raise ValueError(f"hessian 形状应为 ({3*n},{3*n})，得到 {h.shape}")

    basis = _tr_basis(pos, m)
    sqrt_m3 = np.repeat(np.sqrt(m), 3)
    hw = h / np.outer(sqrt_m3, sqrt_m3)
    hw = 0.5 * (hw + hw.T)
    if basis:
        q_tr = np.array(basis).T                      # (3N,k)
        p = np.eye(3 * n) - q_tr @ q_tr.T
        hw = p @ hw @ p

    eigenvalues = np.linalg.eigvalsh(hw)
    # 投影子空间留下 |λ|≈0 的特征值；按 |λ| 选出 3N−k 个振动特征值，
    # 这样负的振动特征值（虚频）不会被误删
    n_vib = 3 * n - len(basis)
    order = np.argsort(np.abs(eigenvalues))
    vib = np.sort(eigenvalues[order[len(basis):]])
    signs = np.where(vib < 0.0, -1.0, 1.0)
    return signs * np.sqrt(np.abs(vib) * 1e24) / (2.0 * math.pi * C_CM_S)


def check_hessian(
    h: CartesianHessian,
    masses: np.ndarray | None = None,
    *,
    allow_nonstationary: bool = False,
) -> HessianDiagnostics:
    """按计划 §2.3 检查 Hessian 质量；通过返回 HessianDiagnostics，失败抛错。

    - 不对称度 >1e-3 或平移残差 >1e-3 → HessianQualityError；
    - 不对称度 >1e-6 → 警告并用 ½(H+Hᵀ) 继续；平移残差 >1e-5 → 警告；
    - 有梯度且 max|g| > 25 kJ/mol/nm → HessianQualityError，
      除非 allow_nonstationary=True（此时只记录进警告）；
    - 最小频率 < −20 cm⁻¹ → HessianQualityError（鞍点）；−20…0 → 警告。

    masses 不给时用 RDKit 平均原子量；有 IR 时应传 IR 的 mass。
    """
    warnings: list[str] = []
    hmax = float(np.abs(h.hessian).max())
    if hmax == 0.0:
        raise HessianQualityError("Hessian 全零，无法检查")

    asymmetry = float(np.abs(h.hessian - h.hessian.T).max()) / hmax
    if asymmetry > _ASYMMETRY_ERROR:
        raise HessianQualityError(
            f"Hessian asymmetry（不对称度）{asymmetry:.3e} 超过报错阈值 1e-3"
        )
    work = 0.5 * (h.hessian + h.hessian.T)
    if asymmetry > _ASYMMETRY_WARN:
        warnings.append(
            f"Hessian asymmetry（不对称度）{asymmetry:.3e} 超过警告阈值 1e-6，"
            "已用 ½(H+Hᵀ) 对称化"
        )

    n = h.positions_nm.shape[0]
    # Σ_j H[3i+a, 3j+b] 对每个 (i,a,b)：平移不变性要求全部为 0
    block_sums = work.reshape(n, 3, n, 3).sum(axis=2)
    translation = float(np.abs(block_sums).max()) / hmax
    if translation > _TRANSLATION_ERROR:
        raise HessianQualityError(
            f"Hessian translation residual（平移残差）{translation:.3e} "
            "超过报错阈值 1e-3：破坏平移不变性"
        )
    if translation > _TRANSLATION_WARN:
        warnings.append(
            f"Hessian translation residual（平移残差）{translation:.3e} "
            "超过警告阈值 1e-5"
        )

    max_gradient: float | None = None
    if h.gradient is not None:
        max_gradient = float(np.abs(h.gradient).max())
        if max_gradient > _STATIONARY_MAX_GRADIENT:
            message = (
                f"non-stationary geometry（非驻点）：max|gradient| = "
                f"{max_gradient:.3e} kJ/mol/nm 超过 25"
            )
            if allow_nonstationary:
                warnings.append(message + "（allow_nonstationary=True，仅记录）")
            else:
                raise HessianQualityError(
                    message + "；如确需继续，传 allow_nonstationary=True"
                )

    if masses is None:
        masses = average_masses(h.atomic_numbers)
    frequencies = harmonic_frequencies(work, h.positions_nm, np.asarray(masses, dtype=np.float64))
    n_imaginary = int((frequencies < 0.0).sum())
    if frequencies.size and frequencies[0] < _IMAGINARY_ERROR_CM:
        raise HessianQualityError(
            f"{n_imaginary} 个 imaginary frequencies（虚频）低于 −20 cm⁻¹"
            f"（最小 {frequencies[0]:.2f} cm⁻¹）：参考几何可能是鞍点"
        )
    if n_imaginary:
        warnings.append(
            f"{n_imaginary} 个 imaginary frequencies（虚频），最小 "
            f"{frequencies[0]:.2f} cm⁻¹，介于 −20…0 cm⁻¹"
        )

    return HessianDiagnostics(
        asymmetry_rel=asymmetry,
        translation_residual_rel=translation,
        max_gradient=max_gradient,
        frequencies_cm=frequencies,
        n_imaginary=n_imaginary,
        warnings=tuple(warnings),
    )
