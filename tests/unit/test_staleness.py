"""staleness（bonded/staleness.py，C3/C4）的必测用例：快照、JSON 往返、
global_fit 依赖记录与过期判定、MLP MM 扣除依赖。纯内存，无外部程序。"""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pytest

from ffcompiler.bonded.global_fit import global_fit
from ffcompiler.bonded.internal import wilson_b
from ffcompiler.bonded.staleness import (
    h_fixed_snapshot,
    mm_subtraction_snapshot,
    staleness,
)
from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.ir import Bond, ForceFieldIR, ProperTorsion
from ffcompiler.ir.schema import Atom


def _diatomic(k_ir: float = 5.0, k_hessian: float = 5.0) -> tuple[ForceFieldIR, CartesianHessian]:
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


def test_snapshot_json_roundtrip_stable() -> None:
    ir, _ = _diatomic()
    snap = h_fixed_snapshot(ir)
    assert json.loads(json.dumps(snap)) == snap  # JSON 可往返（float 逐位一致）
    assert mm_subtraction_snapshot(ir)["bonds"] == [[0, 1, 0.1, 5.0]]


def test_global_fit_records_dependency_and_is_fresh() -> None:
    ir, h = _diatomic()
    new_ir, _, report = global_fit(ir, h, prior_ir=ir, alpha=1e-3, symmetrize=False)
    assert report.n_at_bound == 0
    dep = new_ir.provenance["bonded"]["h_fixed_dependency"]
    assert dep["inputs"] == h_fixed_snapshot(ir)
    assert dep["hessian"]["provider"] == "test"
    assert dep["fd_step_nm"] == 1e-5
    assert new_ir.provenance["bonded"]["method"] == "global_fit"
    assert staleness(new_ir) == []  # 刚拟合完不过期


def test_charge_change_makes_global_fit_stale() -> None:
    ir, h = _diatomic()
    new_ir, _, _ = global_fit(ir, h, prior_ir=ir, alpha=1e-3, symmetrize=False)
    q = [0.05, -0.05]
    charged = new_ir.replace_terms(
        atoms=tuple(replace(a, partial_charge=qi) for a, qi in zip(new_ir.atoms, q)),
        provenance_update={"charge": {"method": "loaded"}},
    )
    items = staleness(charged)
    assert any("电荷" in it for it in items)
    assert all(it.startswith("global_fit:") for it in items)


def test_torsion_change_makes_global_fit_stale() -> None:
    ir, h = _diatomic()
    new_ir, _, _ = global_fit(ir, h, prior_ir=ir, alpha=1e-3, symmetrize=False)
    with_torsion = new_ir.replace_terms(
        proper_torsions=(ProperTorsion(0, 0, 1, 0, 1, 0.0, 1.0, "mlp"),)
    )
    assert any("proper 二面角" in it for it in staleness(with_torsion))


def test_non_global_fit_ir_has_no_staleness_items() -> None:
    ir, _ = _diatomic()
    assert staleness(ir) == []


def test_missing_dependency_recorded_as_stale() -> None:
    ir, h = _diatomic()
    new_ir, _, _ = global_fit(ir, h, prior_ir=ir, alpha=1e-3, symmetrize=False)
    prov = dict(new_ir.provenance)
    bonded = dict(prov["bonded"])
    bonded.pop("h_fixed_dependency")
    legacy = new_ir.replace_terms(provenance_update={"bonded": bonded})
    assert any("依赖未记录" in it for it in staleness(legacy))


def test_mlp_torsion_mm_dependency_staleness() -> None:
    ir, _ = _diatomic()
    torsion = ProperTorsion(0, 0, 1, 0, 1, 0.0, 1.0, "mlp")
    recorded = ir.replace_terms(
        proper_torsions=(torsion,),
        provenance_update={"torsion": {"provider": "mlp",
                                       "mm_dependency": mm_subtraction_snapshot(ir)}},
    )
    assert staleness(recorded) == []
    softened = recorded.replace_terms(
        bonds=(Bond(i=0, j=1, r0=0.1, k=2.5, source="t"),)
    )
    items = staleness(softened)
    assert any("MM 扣除的 bond" in it for it in items)

    no_dep = ir.replace_terms(
        proper_torsions=(torsion,), provenance_update={"torsion": {"provider": "mlp"}}
    )
    assert any("MM 扣除依赖未记录" in it for it in staleness(no_dep))


def test_improper_change_makes_mlp_torsion_stale() -> None:
    # improper 在 mm_energy_without_bond 的扣除能量里（2026-10-17 验收补进快照）
    from ffcompiler.ir import ImproperTorsion

    ir, _ = _diatomic()
    torsion = ProperTorsion(0, 0, 1, 0, 1, 0.0, 1.0, "mlp")
    recorded = ir.replace_terms(
        proper_torsions=(torsion,),
        provenance_update={"torsion": {"provider": "mlp",
                                       "mm_dependency": mm_subtraction_snapshot(ir)}},
    )
    with_improper = recorded.replace_terms(
        impropers=(ImproperTorsion(0, 1, 0, 1, 2, 3.141592653589793, 4.6, "gaff2"),)
    )
    assert any("MM 扣除的 improper" in it for it in staleness(with_improper))


def test_fresh_global_fit_survives_ir_json_roundtrip(tmp_path) -> None:
    from ffcompiler.ir import io as ir_io

    ir, h = _diatomic()
    new_ir, _, _ = global_fit(ir, h, prior_ir=ir, alpha=1e-3, symmetrize=False)
    ir_io.to_json(new_ir, tmp_path / "ir.json")
    assert staleness(ir_io.from_json(tmp_path / "ir.json")) == []
