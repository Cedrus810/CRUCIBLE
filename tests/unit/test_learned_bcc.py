"""learned BCC 骨架（charge/learned_bcc.py，MVP-4 管线）的必测用例。

只测管线：守恒、符号约定、错误路径、provenance。不涉及任何训练。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeError, ChargeProvider
from ffcompiler.charge.learned_bcc import IdentityCorrections, LearnedBCCCharges
from ffcompiler.ir import Bond, ForceFieldIR
from ffcompiler.ir.schema import Atom


class _FakeModel:
    def __init__(self, deltas):
        self._deltas = np.asarray(deltas, dtype=np.float64)

    def predict(self, ir, positions_nm, mol):
        return self._deltas

    def identity(self):
        return {"name": "fake", "trained": False, "n_params": 0}


def _diatomic_ir(partial: tuple[float, float] = (0.0, 0.0)) -> ForceFieldIR:
    atoms = tuple(
        Atom(index=i, name=f"A{i}", element="H", atomic_number=1, mass=1.008,
             atom_type="T", formal_charge=0, partial_charge=q, sigma=0.3, epsilon=0.0)
        for i, q in enumerate(partial)
    )
    return ForceFieldIR(
        name="diatomic", atoms=atoms,
        bonds=(Bond(i=0, j=1, r0=0.1, k=5.0, source="t"),),
        angles=(), proper_torsions=(), impropers=(),
        exclusions=((0, 1),), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )


class _MolStub:
    pass


class _FixedBase(ChargeProvider):
    """固定的基础电荷 q⁽⁰⁾（与 IR 当前电荷无关）。"""

    method = "fixed-base"

    def __init__(self, q):
        self.q = np.asarray(q, dtype=np.float64)

    def compute(self, ir, positions_nm, mol):
        return self.q.copy()

    def source_info(self):
        return {"kind": "test"}


_BASE0 = _FixedBase([0.0, 0.0])
_BASE2 = _FixedBase([0.02, -0.02])


def test_identity_model_returns_base_charges() -> None:
    ir = _diatomic_ir((0.02, -0.02))
    provider = LearnedBCCCharges(IdentityCorrections(), base_provider=_BASE2)
    out = provider.compute(ir, np.zeros((2, 3)), _MolStub())
    assert out.tolist() == [0.02, -0.02]
    assert provider.source_info()["model"]["name"] == "identity-untrained"


def test_bond_antisymmetry_and_conservation() -> None:
    ir = _diatomic_ir((0.0, 0.0))
    provider = LearnedBCCCharges(_FakeModel([0.25]), base_provider=_BASE0)
    out = provider.compute(ir, np.zeros((2, 3)), _MolStub())
    # 符号约定：键 (0,1) 书写顺序 → i=0 端 +Δq，j=1 端 −Δq
    assert out.tolist() == pytest.approx([0.25, -0.25])
    assert out.sum() == pytest.approx(0.0, abs=1e-15)


def test_wrong_shape_rejected() -> None:
    ir = _diatomic_ir()
    provider = LearnedBCCCharges(_FakeModel([0.1, 0.2]), base_provider=_BASE0)  # 2 个 Δq，只有 1 条键
    with pytest.raises(ChargeError, match="每键"):
        provider.compute(ir, np.zeros((2, 3)), _MolStub())


def test_nonfinite_rejected() -> None:
    ir = _diatomic_ir()
    provider = LearnedBCCCharges(_FakeModel([float("nan")]), base_provider=_BASE0)
    with pytest.raises(ChargeError, match="非有限"):
        provider.compute(ir, np.zeros((2, 3)), _MolStub())


def test_missing_artifact_rejected_at_construction(tmp_path: Path) -> None:
    with pytest.raises(ChargeError, match="未立项|工件不存在"):
        LearnedBCCCharges(IdentityCorrections(), base_provider=_BASE0,
                          model_path=tmp_path / "nope.pt")


def test_base_provider_is_required() -> None:
    # S1（base charge）未定：不设"取 IR 当前电荷"的默认，否则重复赋电荷会叠加修正
    with pytest.raises(TypeError):
        LearnedBCCCharges(IdentityCorrections())


def test_assign_charges_records_learned_bcc_provenance() -> None:
    ir = _diatomic_ir((0.02, -0.02))
    provider = LearnedBCCCharges(_FakeModel([0.1]), base_provider=_BASE2)
    new_ir, report = assign_charges(ir, provider, np.zeros((2, 3)), _MolStub())
    assert report.method == "learned_bcc"
    assert new_ir.provenance["charge"]["method"] == "learned_bcc"
    assert new_ir.provenance["charge"]["previous"] is None
    assert [a.partial_charge for a in new_ir.atoms] == pytest.approx([0.12, -0.12])
    assert new_ir.provenance["charge"]["source"]["base"]["method"] == "fixed-base"
    # 上一来源在再次赋电荷时保留（T40 行为一致）；base 固定，修正不叠加
    again, _ = assign_charges(
        new_ir, LearnedBCCCharges(IdentityCorrections(), base_provider=_BASE2),
        np.zeros((2, 3)), _MolStub())
    assert again.provenance["charge"]["previous"]["method"] == "learned_bcc"
    assert [a.partial_charge for a in again.atoms] == pytest.approx([0.02, -0.02])
