"""FF Hessian、几何优化与 Kabsch RMSD（MVP-1 计划 T22，validation 侧）。

FF Hessian 用 OpenMM Reference 的解析力做中心差分（步长默认 1e-5 nm），
结果对称化。全部 float64、只用 Reference 平台（MVP-0 §1.3）。
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import openmm
import openmm.app
from scipy.optimize import minimize
from openmm import unit

from ffcompiler.bonded.internal import linear_angles
from ffcompiler.exporters.openmm import to_openmm_system
from ffcompiler.ir import ForceFieldIR


def _reference_context(ir: ForceFieldIR) -> openmm.Context:
    # 数值工具不走完整校验：裁剪 IR（H_fixed）不满足完整 validate（T22）
    system = to_openmm_system(ir, validate_ir=False)
    platform = openmm.Platform.getPlatformByName("Reference")
    return openmm.Context(system, openmm.VerletIntegrator(0.001), platform)


def _forces_at(context: openmm.Context, positions_nm: np.ndarray) -> np.ndarray:
    context.setPositions(
        [openmm.Vec3(*p) for p in positions_nm] * unit.nanometer
    )
    state = context.getState(getForces=True)
    forces = state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    return np.asarray(forces, dtype=np.float64).ravel()


def ff_hessian(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    *,
    exclude_terms: Literal["none", "bonds_angles"] = "none",
    step_nm: float = 1e-5,
) -> np.ndarray:
    """IR 的 FF Hessian（中心差分解析力），形状 (3N,3N)，已对称化。

    exclude_terms="bonds_angles"：只去掉待推导的 bond/angle（global_fit 的固定
    部分）；在 positions_nm 上线性（> 175°）的键角不推导，留在体系里（T22）。
    """
    pos = np.asarray(positions_nm, dtype=np.float64)
    if exclude_terms == "bonds_angles":
        linear = set(linear_angles(pos, [(a.i, a.j, a.k_atom) for a in ir.angles]))
        ir = ir.replace_terms(
            bonds=(), angles=tuple(a for a in ir.angles if (a.i, a.j, a.k_atom) in linear)
        )
    elif exclude_terms != "none":
        raise ValueError(f"未知 exclude_terms {exclude_terms!r}")
    context = _reference_context(ir)
    n3 = 3 * pos.shape[0]
    hessian = np.zeros((n3, n3))
    for a in range(n3):
        plus = pos.copy()
        plus[a // 3, a % 3] += step_nm
        minus = pos.copy()
        minus[a // 3, a % 3] -= step_nm
        f_plus = _forces_at(context, plus)
        f_minus = _forces_at(context, minus)
        # H[a,b] = −∂F_b/∂x_a
        hessian[a] = -(f_plus - f_minus) / (2.0 * step_nm)
    return 0.5 * (hessian + hessian.T)


def minimize_ff(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    *,
    tolerance_kj_nm: float = 1e-3,
) -> np.ndarray:
    """FF 能量极小化（L-BFGS-B，解析力作梯度），返回优化后坐标 (N,3)。

    收敛判据是最大力分量 ≤ tolerance_kj_nm；达不到就抛 RuntimeError，
    不静默返回半收敛的几何。L-BFGS-B 的能量相对下降判据（ftol）关闭，
    否则会在力远未达标时提前停止。
    """
    context = _reference_context(ir)

    def energy(x: np.ndarray) -> float:
        context.setPositions(
            [openmm.Vec3(*p) for p in x.reshape(-1, 3)] * unit.nanometer
        )
        return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

    def gradient(x: np.ndarray) -> np.ndarray:
        return -_forces_at(context, x.reshape(-1, 3))

    x0 = np.asarray(positions_nm, dtype=np.float64).ravel()
    result = minimize(
        energy, x0, jac=gradient, method="L-BFGS-B",
        options={"gtol": tolerance_kj_nm, "ftol": 0.0, "maxiter": 20000},
    )
    max_force = float(np.abs(gradient(result.x)).max())
    if max_force > tolerance_kj_nm:
        raise RuntimeError(
            f"minimize_ff 未收敛：max|F| = {max_force:.3e} kJ/mol/nm > {tolerance_kj_nm:.1e}"
            f"（{result.message}，{result.nit} 步）"
        )
    return result.x.reshape(-1, 3)


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """全原子 Kabsch 对齐后的 RMSD（nm）。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a_centered = a - a.mean(axis=0)
    b_centered = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a_centered.T @ b_centered)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    aligned = a_centered @ rotation
    diff = aligned - b_centered
    return float(np.sqrt((diff**2).sum() / len(a)))
