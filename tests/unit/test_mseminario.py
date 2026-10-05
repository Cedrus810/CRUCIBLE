"""T17：mSeminario 的必测用例（oracle 为真实断言，决策 0011）。"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.mseminario import mseminario
from ffcompiler.bonded.sobtop_ref import read_sobtop_params
from ffcompiler.hessian.base import CartesianHessian

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from sobtop_oracle import GATE_CASES, REFS, load_case, oracle_mismatches  # noqa: E402

K_SPRING = 5e5  # kJ/mol/nm²


def _diatomic(k: float = K_SPRING) -> CartesianHessian:
    pos = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    h = np.zeros((6, 6))
    h[0, 0] = h[3, 3] = k
    h[0, 3] = h[3, 0] = -k
    return CartesianHessian(
        atomic_numbers=np.array([6, 6]), positions_nm=pos, hessian=h,
        gradient=None, energy_kj=None, provider="test", method="analytic",
    )


def test_diatomic_spring_exact() -> None:
    h = _diatomic()
    terms = mseminario(h, [(0, 1)], [])
    assert len(terms) == 1
    assert terms[0].k == pytest.approx(K_SPRING, rel=1e-10)
    assert terms[0].value == pytest.approx(0.1)


def test_atom_reordering_invariance() -> None:
    h, bonds, angles = load_case("ethanol_gfn2")
    terms = mseminario(h, bonds, angles)
    rng = np.random.default_rng(20261004)
    perm = rng.permutation(len(h.atomic_numbers))
    inverse = np.argsort(perm)
    order = (3 * perm[:, None] + np.arange(3)).ravel()
    h2 = CartesianHessian(
        atomic_numbers=h.atomic_numbers[perm],
        positions_nm=h.positions_nm[perm],
        hessian=h.hessian[np.ix_(order, order)],
        gradient=None, energy_kj=None, provider="test", method="gfn2",
    )
    # 只换下标、不改键方向：Sobtop mSeminario 的键 k 取 K_ij 一侧（决策 0011）
    bonds2 = [tuple(int(inverse[a]) for a in b) for b in bonds]
    angles2 = []
    for i, j, k_at in angles:
        angles2.append((int(inverse[i]), int(inverse[j]), int(inverse[k_at])))
    terms2 = mseminario(h2, bonds2, angles2)

    bond_map = {
        tuple(sorted(t.atoms)): idx
        for idx, t in enumerate(terms)
        if len(t.atoms) == 2
    }
    angle_map = {}
    for idx, t in enumerate(terms):
        if len(t.atoms) == 3:
            i, j, k = t.atoms
            angle_map[(min(i, k), j, max(i, k))] = idx
    for t in terms2:
        if len(t.atoms) == 2:
            old_atoms = tuple(sorted(int(perm[a]) for a in t.atoms))
            ref = terms[bond_map[old_atoms]]
        else:
            oi, oj, ok = (int(perm[a]) for a in t.atoms)
            ref = terms[angle_map[(min(oi, ok), oj, max(oi, ok))]]
        assert t.k == pytest.approx(ref.k, rel=1e-9)
        assert t.value == pytest.approx(ref.value, rel=1e-9)


def test_oracle_sobtop_mseminario() -> None:
    """用例 1：除 sf6_fchk 外全部 case 与 param_m2.txt 一致（§2.8 容差，含平衡值）。"""
    failures = {}
    for case in GATE_CASES:
        h, bonds, angles = load_case(case)
        bad = oracle_mismatches(
            mseminario(h, bonds, angles), read_sobtop_params(REFS / case / "param_m2.txt")
        )
        if bad:
            failures[case] = bad[:5]
    assert not failures, failures
