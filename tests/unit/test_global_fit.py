"""T22：全局 Hessian 拟合（MVP-1b）的必测用例。

自洽构造按计划 T22 用例 2：ethanol GAFF2 IR 的 bond/angle 平衡值设成 x₀ 上的值，
H_ref = 该 IR 的完整 FF Hessian（含 torsion、非键），以 CartesianHessian 传入。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.bonded.global_fit import global_fit
from ffcompiler.bonded.internal import angle_value, bond_value, wilson_b
from ffcompiler.hessian.base import CartesianHessian, HessianQualityError
from ffcompiler.hessian.xtb import read_xtb_hessian
from ffcompiler.importers import from_amber
from ffcompiler.ir import Angle, Atom, Bond, ForceFieldIR
from ffcompiler.ir.topology import bond_graph
from ffcompiler.validation.vibrations import ff_hessian

REPO = Path(__file__).resolve().parents[2]
GAFF2 = REPO / "tests" / "fixtures" / "gaff2"
XTB_FIXTURE = REPO / "tests" / "fixtures" / "hessian" / "xtb_ethanol_gfn2"
SEED = 20261004


def _load_ir(name: str) -> tuple[ForceFieldIR, np.ndarray]:
    mdir = GAFF2 / name
    return from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )


def _at_x0(ir: ForceFieldIR, x0: np.ndarray) -> ForceFieldIR:
    """把全部 bond/angle 平衡值设成 x₀ 上的值（q⁰ = q(x₀)）。"""
    return ir.replace_terms(
        bonds=tuple(
            Bond(i=b.i, j=b.j, r0=bond_value(x0, b.i, b.j), k=b.k,
                 source=b.source, meta=dict(b.meta))
            for b in ir.bonds
        ),
        angles=tuple(
            Angle(i=a.i, j=a.j, k_atom=a.k_atom,
                  theta0=angle_value(x0, a.i, a.j, a.k_atom),
                  force_constant=a.force_constant, source=a.source, meta=dict(a.meta))
            for a in ir.angles
        ),
    )


def _ff_reference(ir: ForceFieldIR, x0: np.ndarray) -> CartesianHessian:
    return CartesianHessian(
        atomic_numbers=np.array([a.atomic_number for a in ir.atoms]),
        positions_nm=x0, hessian=ff_hessian(ir, x0), gradient=None, energy_kj=None,
        provider="test", method="ff",
    )


def _ethanol_fixture() -> tuple[ForceFieldIR, np.ndarray, CartesianHessian]:
    ir, _ = _load_ir("ethanol")
    x0 = read_xtb_hessian(
        XTB_FIXTURE / "hessian", XTB_FIXTURE / "xtbopt.xyz", method="gfn2"
    ).positions_nm
    ir0 = _at_x0(ir, x0)
    return ir0, x0, _ff_reference(ir0, x0)


def _rng_priors(ir0: ForceFieldIR) -> ForceFieldIR:
    rng = np.random.default_rng(SEED)
    new_bonds = tuple(
        Bond(i=b.i, j=b.j, r0=b.r0, k=b.k * rng.uniform(0.5, 1.5),
             source=b.source, meta=dict(b.meta))
        for b in ir0.bonds
    )
    new_angles = tuple(
        Angle(i=a.i, j=a.j, k_atom=a.k_atom, theta0=a.theta0,
              force_constant=a.force_constant * rng.uniform(0.5, 1.5),
              source=a.source, meta=dict(a.meta))
        for a in ir0.angles
    )
    return ir0.replace_terms(bonds=new_bonds, angles=new_angles)


def _worst_rel_error(fitted: ForceFieldIR, truth_ir: ForceFieldIR) -> float:
    truth = {("b", b.i, b.j): b.k for b in truth_ir.bonds}
    truth.update({("a", a.i, a.j, a.k_atom): a.force_constant for a in truth_ir.angles})
    worst = 0.0
    for b in fitted.bonds:
        ref = truth[("b", b.i, b.j)]
        worst = max(worst, abs(b.k - ref) / ref)
    for a in fitted.angles:
        ref = truth[("a", a.i, a.j, a.k_atom)]
        worst = max(worst, abs(a.force_constant - ref) / ref)
    return worst


def test_global_fit_self_consistent_recovery() -> None:
    """用例 2：α=0 时 symmetrize 关/开都精确恢复真值（≤1e-5）。"""
    ir0, x0, h = _ethanol_fixture()
    prior_ir = _rng_priors(ir0)
    for sym in (False, True):
        new_ir, _, report = global_fit(ir0, h, prior_ir=prior_ir, alpha=0.0, symmetrize=sym)
        worst = _worst_rel_error(new_ir, ir0)
        print(f"symmetrize={sym}: 恢复最大相对误差 {worst:.3e}")
        assert worst <= 1e-5
        assert report.n_at_bound == 0


def test_fit_uses_reference_hessian() -> None:
    """拟合目标必须是 h.hessian：不扣固定部分、α=0 时 H_ref×2 ⇒ 全部 k×2。"""
    ir0, x0, h = _ethanol_fixture()
    h2 = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=x0, hessian=2.0 * h.hessian,
        gradient=None, energy_kj=None, provider="test", method="ff",
    )
    kw = dict(prior_ir=_rng_priors(ir0), alpha=0.0, symmetrize=False, subtract_fixed=False)
    fit1, _, _ = global_fit(ir0, h, **kw)
    fit2, _, _ = global_fit(ir0, h2, **kw)
    for b1, b2 in zip(fit1.bonds, fit2.bonds):
        assert b2.k == pytest.approx(2.0 * b1.k, rel=1e-9)


def test_symmetrize_false_fits_each_term_separately() -> None:
    """§2.9：symmetrize=False 时每项自成一类——等价项真值不同也要各自恢复。"""
    ir, _ = _load_ir("ethanol")
    _, x0, _ = _ethanol_fixture()
    graph = bond_graph(ir)
    carbon0 = next(i for i in graph if len(graph[i]) == 4
                   and sum(ir.atoms[j].atomic_number == 1 for j in graph[i]) == 3)
    target = next(b for b in ir.bonds if carbon0 in (b.i, b.j)
                  and ir.atoms[b.j if b.i == carbon0 else b.i].atomic_number == 1)
    bumped = ir.replace_terms(bonds=tuple(
        Bond(i=b.i, j=b.j, r0=b.r0, k=b.k * (1.3 if b is target else 1.0),
             source=b.source, meta=dict(b.meta))
        for b in ir.bonds
    ))
    ir0 = _at_x0(bumped, x0)
    new_ir, _, _ = global_fit(
        ir0, _ff_reference(ir0, x0), prior_ir=_rng_priors(ir0), alpha=0.0, symmetrize=False
    )
    assert _worst_rel_error(new_ir, ir0) <= 1e-5


def test_subtract_fixed_necessity() -> None:
    """用例 3：同一构造下不扣除固定部分，恢复误差 > 1e-2（全部 k，与用例 2 同口径）。"""
    ir0, x0, h = _ethanol_fixture()
    prior_ir = _rng_priors(ir0)
    new_no, _, report_no = global_fit(
        ir0, h, prior_ir=prior_ir, alpha=0.0, symmetrize=False, subtract_fixed=False,
    )
    worst_no = _worst_rel_error(new_no, ir0)
    print(f"不扣除固定部分：恢复最大相对误差 {worst_no:.3e}")
    assert worst_no > 1e-2
    assert report_no.subtract_fixed is False
    assert all(b.meta["subtract_fixed"] is False for b in new_no.bonds)


def _diatomic(k_ir: float, k_hessian: float) -> tuple[ForceFieldIR, CartesianHessian]:
    pos = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    b_row = wilson_b(pos, [("bond", (0, 1))])[0]
    atoms = tuple(
        Atom(index=i, name=f"A{i}", element="H", atomic_number=1, mass=1.008,
             atom_type="T", formal_charge=0, partial_charge=0.0, sigma=0.3, epsilon=0.0)
        for i in range(2)
    )
    ir = ForceFieldIR(
        name="diatomic", atoms=atoms,
        bonds=(Bond(i=0, j=1, r0=0.1, k=k_ir, source="t"),),
        angles=(), proper_torsions=(), impropers=(),
        exclusions=((0, 1),), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )
    h = CartesianHessian(
        atomic_numbers=np.array([1, 1]), positions_nm=pos,
        hessian=k_hessian * np.outer(b_row, b_row),
        gradient=np.zeros(6), energy_kj=None, provider="test", method="synthetic",
    )
    return ir, h


def test_lower_bound_hit() -> None:
    """用例 4：H_ref 给负曲率（−5，虚频 ≈ −17 cm⁻¹ 只警告），先验 +5 ⇒ 无约束解 s=−1，拟合后 s=0。"""
    ir, h = _diatomic(k_ir=5.0, k_hessian=-5.0)
    new_ir, _, report = global_fit(ir, h, prior_ir=ir, alpha=0.0, symmetrize=False)
    assert report.n_at_bound == 1
    assert new_ir.bonds[0].k == pytest.approx(0.0, abs=1e-12)
    assert np.isfinite(report.max_force_at_x0)


def test_report_fields_and_meta() -> None:
    """用例 5 + T22 规格：报告字段有限，meta 含 prior_method、alpha、scale、subtract_fixed。"""
    ir0, x0, h = _ethanol_fixture()
    new_ir, _, report = global_fit(ir0, h, prior_ir=_rng_priors(ir0), alpha=1e-3, symmetrize=True)
    assert np.isfinite(report.cond) and report.cond > 0
    assert np.isfinite(report.rel_residual_fit) and report.rel_residual_fit >= 0
    assert np.isfinite(report.rel_residual_prior) and report.rel_residual_prior >= 0
    assert report.rel_residual_fit <= report.rel_residual_prior + 1e-12
    assert np.isfinite(report.max_force_at_x0) and report.max_force_at_x0 > 0
    assert report.lam > 0 and report.alpha > 0
    assert report.subtract_fixed is True
    for term in (*new_ir.bonds, *new_ir.angles):
        assert term.source == "global_fit"
        assert term.meta["prior_method"] == "prior_ir"
        assert term.meta["alpha"] == 1e-3
        assert term.meta["subtract_fixed"] is True
        assert np.isfinite(term.meta["scale"])


def test_symmetrize_averages_equilibrium_values() -> None:
    """§2.9 第 1 条：对称化时平衡值取类内平均（CH3 三个 C-H 的 r0 相同）。"""
    ir, positions = _load_ir("ethanol")
    h = read_xtb_hessian(XTB_FIXTURE / "hessian", XTB_FIXTURE / "xtbopt.xyz", method="gfn2")
    new_ir, _, _ = global_fit(ir, h, prior_method="m2seminario", symmetrize=True,
                              input_positions_nm=positions)
    graph = bond_graph(ir)
    carbon0 = next(i for i in graph if len(graph[i]) == 4
                   and sum(ir.atoms[j].atomic_number == 1 for j in graph[i]) == 3)
    r0s = {b.r0 for b in new_ir.bonds if carbon0 in (b.i, b.j)
           and ir.atoms[b.j if b.i == carbon0 else b.i].atomic_number == 1}
    assert len(r0s) == 1


def test_topology_checks() -> None:
    """Review Focus 1：原子顺序打乱、或键长与输入几何不符 → TopologyMismatchError。"""
    ir, positions = _load_ir("ethanol")
    h = read_xtb_hessian(XTB_FIXTURE / "hessian", XTB_FIXTURE / "xtbopt.xyz", method="gfn2")
    perm = np.arange(len(h.atomic_numbers))[::-1]
    order = (3 * perm[:, None] + np.arange(3)).ravel()
    shuffled = CartesianHessian(
        atomic_numbers=h.atomic_numbers[perm], positions_nm=h.positions_nm[perm],
        hessian=h.hessian[np.ix_(order, order)], gradient=None, energy_kj=None,
        provider="xtb", method="gfn2",
    )
    with pytest.raises(TopologyMismatchError):
        global_fit(ir, shuffled, prior_method="m2seminario")
    stretched_input = positions.copy()
    stretched_input[2] = positions[1] + 1.5 * (positions[2] - positions[1])
    with pytest.raises(TopologyMismatchError, match=r"1, 2"):
        global_fit(ir, h, prior_method="m2seminario", input_positions_nm=stretched_input)


def test_invalid_linear_angle_policy_rejected() -> None:
    ir0, _, h = _ethanol_fixture()
    with pytest.raises(ValueError, match="linear_angle_policy"):
        global_fit(ir0, h, prior_ir=ir0, linear_angle_policy="ignore")


def test_nonstationary_rejected() -> None:
    """用例 6（Review Focus 3）。"""
    ir0, x0, h = _ethanol_fixture()
    h_bad = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=x0,
        hessian=h.hessian, gradient=np.full(h.hessian.shape[0], 100.0),
        energy_kj=None, provider="test", method="ff",
    )
    with pytest.raises(HessianQualityError):
        global_fit(ir0, h_bad, prior_method="m2seminario")


@pytest.mark.parametrize("prior_method", ["mseminario", "m2seminario", "drih"])
def test_prior_methods(prior_method: str) -> None:
    """prior_method 接受全部 BondedMethod（计划 T22 接口）。"""
    ir, positions = _load_ir("ethanol")
    h = read_xtb_hessian(XTB_FIXTURE / "hessian", XTB_FIXTURE / "xtbopt.xyz", method="gfn2")
    new_ir, _, report = global_fit(ir, h, prior_method=prior_method,
                                   input_positions_nm=positions)
    assert report.rel_residual_fit <= report.rel_residual_prior + 1e-12
    assert all(b.meta["prior_method"] == prior_method for b in new_ir.bonds)
