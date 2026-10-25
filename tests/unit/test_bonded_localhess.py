"""localhess（Sobtop 隐藏 k_method=5 复刻，决策 0023，experimental）单元测试。

无 oracle（Sobtop 菜单无该方法入口），测试策略为内部一致性：
B 行解析梯度数值校验、i/j 交换不变性（方法固有对称性）、乙醇 GFN2 上的
量级窗口与正性、线性角/退化块错误路径、derive_bonded 全链路。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.assign import derive_bonded
from ffcompiler.bonded.internal import LinearAngleError
from ffcompiler.bonded.localhess import (
    NonFiniteComplianceError,
    _pinv_cutoff,
    angle_b_row,
    bond_b_row,
    local_compliance_k,
    localhess,
)
from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.hessian.xtb import read_xtb_hessian
from ffcompiler.importers import from_amber

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ETH = FIXTURES / "gaff2" / "ethanol"
XTB = FIXTURES / "hessian" / "xtb_ethanol_gfn2"


def _ethanol():
    ir, pos = from_amber(ETH / "ethanol.prmtop", ETH / "ethanol.inpcrd", ETH / "ethanol.sdf")
    h = read_xtb_hessian(XTB / "hessian", XTB / "xtbopt.xyz", method="gfn2")
    bonds = [(b.i, b.j) for b in ir.bonds]
    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]
    return ir, h, bonds, angles


def _fake_hessian(hess: np.ndarray, n_atoms: int) -> CartesianHessian:
    return CartesianHessian(
        atomic_numbers=np.full(n_atoms, 6),
        positions_nm=np.zeros((n_atoms, 3)),
        hessian=hess,
        gradient=None,
        energy_kj=None,
        provider="test",
        method="test",
    )


def test_bond_b_row_numeric_gradient():
    # r = |x_j - x_i| 的解析 B 行 vs 中心差分
    pos = np.array([[0.0, 0.0, 0.0], [0.15, 0.03, -0.02]])
    b = bond_b_row(pos, 0, 1)
    eps = 1e-7
    num = np.zeros(6)
    for c in range(6):
        p = pos.copy().reshape(-1)
        p[c] += eps
        m = pos.copy().reshape(-1)
        m[c] -= eps
        num[c] = (
            np.linalg.norm(p.reshape(2, 3)[1] - p.reshape(2, 3)[0])
            - np.linalg.norm(m.reshape(2, 3)[1] - m.reshape(2, 3)[0])
        ) / (2 * eps)
    assert np.allclose(b, num, atol=1e-5)


def test_angle_b_row_numeric_gradient():
    pos = np.array(
        [[0.1, 0.02, 0.0], [0.0, 0.0, 0.0], [-0.08, 0.12, 0.05]]
    )
    b = angle_b_row(pos, 0, 1, 2)
    eps = 1e-7

    def theta(p):
        u = p[0] - p[1]
        v = p[2] - p[1]
        return np.arccos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

    num = np.zeros(9)
    for c in range(9):
        p = pos.copy().reshape(-1); p[c] += eps
        m = pos.copy().reshape(-1); m[c] -= eps
        num[c] = (theta(p.reshape(3, 3)) - theta(m.reshape(3, 3))) / (2 * eps)
    assert np.allclose(b, num, atol=1e-5)


def test_pinv_cutoff_zeroes_tiny_singulars():
    m = np.diag([2.0, 1.0, 3e-12])
    p = _pinv_cutoff(m)
    assert np.allclose(np.diag(p), [0.5, 1.0, 0.0])


def test_local_compliance_symmetric_under_swap():
    ir, h, bonds, angles = _ethanol()
    i, j = bonds[0]
    k1 = local_compliance_k(h, (i, j), bond_b_row(h.positions_nm, i, j))
    k2 = local_compliance_k(h, (j, i), bond_b_row(h.positions_nm, j, i))
    # 方法固有对称：H 块与 B 行在 (i,j)→(j,i) 下只是行列重排
    assert k1 == pytest.approx(k2, rel=1e-12)


def test_ethanol_localhess_terms_finite_positive():
    ir, h, bonds, angles = _ethanol()
    terms = localhess(h, bonds, angles)
    assert len(terms) == len(bonds) + len(angles)
    for t in terms[: len(bonds)]:
        assert t.k > 0 and 1e5 < t.k < 1e6  # 键 kJ/mol/nm² 量级窗口
    for t in terms[len(bonds):]:
        assert t.k > 0 and 50 < t.k < 5e3   # 角 kJ/mol/rad² 量级窗口


def test_zero_block_raises_nonfinite():
    h = _fake_hessian(np.zeros((6, 6)), 2)
    b = np.array([0.0, 0.0, 1.0, 0.0, 0.0, -1.0])
    with pytest.raises(NonFiniteComplianceError):
        local_compliance_k(h, (0, 1), b)


def test_linear_angle_rejected():
    # 共线三原子：sinθ→0，B 行奇异 → LinearAngleError
    pos = np.zeros((3, 3))
    pos[1] = [0.0, 0.0, 0.0]
    pos[0] = [0.1, 0.0, 0.0]
    pos[2] = [-0.1, 0.0, 0.0]
    with pytest.raises(LinearAngleError):
        angle_b_row(pos, 0, 1, 2)


def test_derive_bonded_pipeline_localhess():
    ir, h, _, _ = _ethanol()
    new_ir, geom, report = derive_bonded(ir, h, method="localhess", symmetrize=True)
    assert report.method == "localhess"
    # 对称化后甲基 C-H 三条键 k 相同
    ch = [b for b in new_ir.bonds if {b.i, b.j} in ({0, 3}, {0, 4}, {0, 5})]
    ks = {round(b.k, 6) for b in ch}
    assert len(ks) == 1
    # 换回的键长来自 Hessian 几何
    assert all(b.source == "localhess" for b in new_ir.bonds)
