"""QMDFF 式全局 Hessian 拟合（MVP-1b，MVP-1 计划 T22，§2.9 锁定）。

模型：H_ref − H_fixed = Σ_c s_c·k_c^prior·(Σ_{t∈c} b_tᵀb_t)，对 s ≥ 0 做
带 λ‖s−1‖² 正则的线性最小二乘（scipy.optimize.lsq_linear，下界 0）。
数学定义与"为什么 bbᵀ 只在 q⁰=q(x₀) 成立"见 docs/theory/hessian_methods.md §5。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import lsq_linear

from ffcompiler.bonded.assign import _derive_terms, check_hessian_topology
from ffcompiler.bonded.internal import (
    LinearAngleError,
    angle_value,
    bond_value,
    linear_angles,
    wilson_b,
)
from ffcompiler.bonded.symmetry import term_classes
from ffcompiler.hessian.base import CartesianHessian, check_hessian
from ffcompiler.ir import Angle, Bond, ForceFieldIR
from ffcompiler.validation.vibrations import _forces_at, _reference_context, ff_hessian

_BONDED_METHODS = ("mseminario", "m2seminario", "drih")


@dataclass(frozen=True)
class GlobalFitReport:
    alpha: float
    lam: float
    cond: float
    rel_residual_fit: float
    rel_residual_prior: float
    n_at_bound: int
    scale_by_class: dict[str, float]
    max_force_at_x0: float  # kJ/mol/nm
    subtract_fixed: bool


def _term_label(atoms: tuple[int, ...]) -> str:
    if len(atoms) == 2:
        return f"bond({atoms[0]},{atoms[1]})"
    return f"angle({atoms[0]},{atoms[1]},{atoms[2]})"


def _prior_ks(
    ir: ForceFieldIR,
    h: CartesianHessian,
    *,
    prior_method: str,
    prior_ir: ForceFieldIR | None,
    symmetrize: bool,
    bonds: list[tuple[int, int]],
    angles: list[tuple[int, int, int]],
    allow_nonstationary: bool,
) -> np.ndarray:
    """每个待推导项的先验 k（½k 约定），与项顺序对齐。"""
    if prior_ir is not None:
        bond_k = {
            (min(b.i, b.j), max(b.i, b.j)): b.k for b in prior_ir.bonds
        }
        angle_k = {
            (min(a.i, a.k_atom), a.j, max(a.i, a.k_atom)): a.force_constant
            for a in prior_ir.angles
        }
        values = []
        for i, j in bonds:
            values.append(bond_k[(min(i, j), max(i, j))])
        for i, j, k_at in angles:
            values.append(angle_k[(min(i, k_at), j, max(i, k_at))])
        return np.array(values, dtype=np.float64)
    if prior_method not in _BONDED_METHODS:
        raise ValueError(f"未知先验方法 {prior_method!r}")
    terms = _derive_terms(h, ir, prior_method, bonds, angles, allow_nonstationary)
    if symmetrize:
        classes = term_classes(ir, [t.atoms for t in terms])
        terms, _ = _symmetrize(terms, classes)
    return np.array([t.k for t in terms], dtype=np.float64)


def _symmetrize(terms, classes):
    from ffcompiler.bonded.symmetry import symmetrize as sym

    return sym(terms, classes)


def global_fit(
    ir: ForceFieldIR,
    h: CartesianHessian,
    *,
    prior_method: str = "mseminario",
    prior_ir: ForceFieldIR | None = None,
    alpha: float = 1e-3,
    symmetrize: bool = True,
    linear_angle_policy: str = "error",
    allow_nonstationary: bool = False,
    fd_step_nm: float = 1e-5,
    input_positions_nm: np.ndarray | None = None,
    subtract_fixed: bool = True,
) -> tuple[ForceFieldIR, np.ndarray, GlobalFitReport]:
    """全局 Hessian 拟合（MVP-1b，§2.9）。返回 (新 IR, Hessian 几何, 报告)。

    拟合目标 y = vec(H_ref − H_fixed)，H_ref 取 h.hessian（对称化），H_fixed 为
    IR 去掉待推导 bond/angle 后在 x₀ 的 FF Hessian（subtract_fixed=False 时为 0，
    仅供诊断）。symmetrize=False 时每项自成一类。线性键角按策略处理：
    "error" 抛 LinearAngleError；"library" 保留原参数并留在固定部分。
    不在此处调用 validate：下界命中的类会得到 k=0（计划 T22 用例 4）。
    """
    if linear_angle_policy not in ("error", "library"):
        raise ValueError(f"未知 linear_angle_policy {linear_angle_policy!r}")
    check_hessian_topology(ir, h, input_positions_nm)
    check_hessian(
        h, np.array([a.mass for a in ir.atoms]), allow_nonstationary=allow_nonstationary
    )
    pos = h.positions_nm
    n3 = 3 * pos.shape[0]

    bonds = [(b.i, b.j) for b in ir.bonds]
    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]
    linear = {(min(i, k), j, max(i, k)) for i, j, k in linear_angles(pos, angles)}
    if linear and linear_angle_policy == "error":
        raise LinearAngleError(sorted(linear))
    angles_for_fit = [
        a for a in angles if (min(a[0], a[2]), a[1], max(a[0], a[2])) not in linear
    ]

    # 先验（逐项，与 coords 顺序对齐）
    k_prior = _prior_ks(
        ir, h, prior_method=prior_method, prior_ir=prior_ir, symmetrize=symmetrize,
        bonds=bonds, angles=angles_for_fit, allow_nonstationary=allow_nonstationary,
    )

    h_ref = 0.5 * (h.hessian + h.hessian.T)
    h_fixed = (
        ff_hessian(ir, pos, exclude_terms="bonds_angles", step_nm=fd_step_nm)
        if subtract_fixed
        else np.zeros((n3, n3))
    )

    # 等价类：symmetrize=False 时每项自成一类（§2.9）
    coords = [("bond", b) for b in bonds] + [("angle", a) for a in angles_for_fit]
    if symmetrize:
        classes = list(term_classes(ir, [c[1] for c in coords]))
    else:
        classes = list(range(len(coords)))
    class_ids = sorted(set(classes))
    class_index = {cid: i for i, cid in enumerate(class_ids)}
    members: dict[int, list[int]] = {cid: [] for cid in class_ids}
    for row, cls in enumerate(classes):
        members[cls].append(row)

    # 设计矩阵：D[:,c] = k_prior_c · vec(Σ_{t∈c} b_t b_tᵀ)
    b_matrix = wilson_b(pos, coords)
    prior_by_class = np.array([np.mean(k_prior[members[cid]]) for cid in class_ids])
    columns = []
    for cid, k_c in zip(class_ids, prior_by_class):
        rows = b_matrix[members[cid]]
        columns.append(k_c * (rows.T @ rows).ravel())
    d_matrix = np.column_stack(columns)
    y = (h_ref - h_fixed).ravel()

    n_c = len(class_ids)
    lam = alpha * float(np.trace(d_matrix.T @ d_matrix)) / n_c
    if alpha > 0.0:
        a_aug = np.vstack([d_matrix, np.sqrt(lam) * np.eye(n_c)])
        b_aug = np.concatenate([y, np.sqrt(lam) * np.ones(n_c)])
    else:
        a_aug, b_aug = d_matrix, y
    scales = lsq_linear(a_aug, b_aug, bounds=(0.0, np.inf), lsq_solver="exact").x

    rel_prior = float(np.linalg.norm(y - d_matrix @ np.ones(n_c)) / np.linalg.norm(y))
    rel_fit = float(np.linalg.norm(y - d_matrix @ scales) / np.linalg.norm(y))
    n_at_bound = int(np.sum(scales <= 1e-10))

    # 写回 IR：k = s_c·k_c^prior；平衡值 q⁰ = q(x₀)，对称化时取类内平均
    values = np.array([
        bond_value(pos, *c[1]) if c[0] == "bond" else angle_value(pos, *c[1]) for c in coords
    ])
    value_by_class = {cid: float(np.mean(values[members[cid]])) for cid in class_ids}
    fitted: dict[tuple[int, ...], tuple[float, float, float]] = {}
    for row, (kind, atoms) in enumerate(coords):
        ci = class_index[classes[row]]
        key = (min(atoms), max(atoms)) if kind == "bond" else (
            min(atoms[0], atoms[2]), atoms[1], max(atoms[0], atoms[2])
        )
        fitted[key] = (
            float(scales[ci] * prior_by_class[ci]),
            value_by_class[classes[row]],
            float(scales[ci]),
        )
    scale_by_class = {
        _term_label(coords[members[cid][0]][1]): float(scales[class_index[cid]])
        for cid in class_ids
    }

    meta_common = {
        "prior_method": "prior_ir" if prior_ir is not None else prior_method,
        "alpha": alpha,
        "symmetrized": symmetrize,
        "subtract_fixed": subtract_fixed,
    }
    new_bonds = []
    for bond in ir.bonds:
        k, r0, scale = fitted[(min(bond.i, bond.j), max(bond.i, bond.j))]
        new_bonds.append(
            Bond(i=bond.i, j=bond.j, r0=r0, k=k, source="global_fit",
                 meta=dict(meta_common, scale=scale))
        )
    new_angles = []
    for angle in ir.angles:
        key = (min(angle.i, angle.k_atom), angle.j, max(angle.i, angle.k_atom))
        if key in fitted:
            k, theta0, scale = fitted[key]
            new_angles.append(
                Angle(i=angle.i, j=angle.j, k_atom=angle.k_atom, theta0=theta0,
                      force_constant=k, source="global_fit",
                      meta=dict(meta_common, scale=scale))
            )
        else:  # library 线性角：原参数保留
            new_angles.append(
                Angle(i=angle.i, j=angle.j, k_atom=angle.k_atom,
                      theta0=angle.theta0, force_constant=angle.force_constant,
                      source=angle.source, meta=dict(angle.meta, linear_angle=True))
            )
    new_ir = ir.replace_terms(bonds=tuple(new_bonds), angles=tuple(new_angles))

    # 拟合后 FF 在 x₀ 的最大力（"平衡位置不一致"的量化）
    max_force = float(np.abs(_forces_at(_reference_context(new_ir), pos)).max())

    report = GlobalFitReport(
        alpha=alpha,
        lam=lam,
        cond=float(np.linalg.cond(d_matrix)),
        rel_residual_fit=rel_fit,
        rel_residual_prior=rel_prior,
        n_at_bound=n_at_bound,
        scale_by_class=scale_by_class,
        max_force_at_x0=max_force,
        subtract_fixed=subtract_fixed,
    )
    return new_ir, pos, report
