"""T18：m2seminario（C2 候选）的必测用例。

按计划 §2.6 命名规则：C2 在全部 oracle case 上与 Sobtop m2Seminario 一致，
公开方法名为 `m2seminario`（决策 0011，取代 0009）。
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.internal import proper_quartets, wilson_b
from ffcompiler.bonded.sobtop_ref import read_sobtop_params
from ffcompiler.bonded.m2seminario import m2seminario
from ffcompiler.hessian.base import CartesianHessian

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from sobtop_oracle import GATE_CASES, REFS, load_case, oracle_mismatches  # noqa: E402


def test_bond_direction_symmetry() -> None:
    h, bonds, angles = load_case("ethanol_gfn2")
    terms_ij = m2seminario(h, bonds, [])
    swapped = [(j, i) for i, j in bonds]
    terms_ji = m2seminario(h, swapped, [])
    for t_ij, t_ji in zip(terms_ij, terms_ji):
        assert t_ij.k == pytest.approx(t_ji.k, rel=1e-12)


def test_angle_endpoint_order_symmetry() -> None:
    """用例 2（angle 部分）：端点交换 (i,j,k)→(k,j,i)，k 完全相同。"""
    h, bonds, angles = load_case("ethanol_gfn2")
    forward = [t for t in m2seminario(h, bonds, angles) if len(t.atoms) == 3]
    flipped = [(k, j, i) for i, j, k in angles]
    backward = [t for t in m2seminario(h, bonds, flipped) if len(t.atoms) == 3]
    for t_f, t_b in zip(forward, backward):
        assert t_b.atoms == tuple(reversed(t_f.atoms))
        assert t_f.k == pytest.approx(t_b.k, rel=1e-12)


def test_equivalent_terms_same_on_synthetic_d6h_benzene() -> None:
    # 严格 D6h 苯：C-C 0.139 nm、C-H 0.108 nm、全部键角 120°
    r_cc, r_ch = 0.139, 0.108
    pos = []
    for n in range(6):
        angle = 2.0 * math.pi * n / 6.0
        pos.append([r_cc * math.cos(angle), r_cc * math.sin(angle), 0.0])
    for n in range(6):
        angle = 2.0 * math.pi * n / 6.0
        pos.append(
            [
                (r_cc + r_ch) * math.cos(angle),
                (r_cc + r_ch) * math.sin(angle),
                0.0,
            ]
        )
    pos = np.array(pos)
    bonds = [(n, (n + 1) % 6) for n in range(6)] + [(n, 6 + n) for n in range(6)]
    angles = []
    for v in range(6):
        nbrs = sorted({(v - 1) % 6, (v + 1) % 6, 6 + v})
        for x in range(len(nbrs)):
            for y in range(x + 1, len(nbrs)):
                angles.append((nbrs[x], v, nbrs[y]))
    # 合成 Hessian：全部 C-C 键 k1、C-H 键 k2、C-C-C 角 k3、C-C-H 角 k4
    k1, k2, k3, k4 = 3000.0, 3500.0, 400.0, 320.0
    coords = [("bond", b) for b in bonds] + [("angle", a) for a in angles]
    b_matrix = wilson_b(pos, coords)
    ks = [k1] * 6 + [k2] * 6 + [
        k3 if all(a < 6 for a in ang) else k4 for ang in angles
    ]
    hessian = sum(k * np.outer(b, b) for k, b in zip(ks, b_matrix))
    h = CartesianHessian(
        atomic_numbers=np.array([6] * 6 + [1] * 6), positions_nm=pos,
        hessian=hessian, gradient=None, energy_kj=None,
        provider="test", method="synthetic",
    )
    terms = m2seminario(h, bonds, angles)

    def spread(values):
        mean = sum(values) / len(values)
        return (max(values) - min(values)) / abs(mean)

    cc = [t.k for t in terms if len(t.atoms) == 2 and all(a < 6 for a in t.atoms)]
    ch = [t.k for t in terms if len(t.atoms) == 2 and any(a >= 6 for a in t.atoms)]
    ccc = [t.k for t in terms if len(t.atoms) == 3 and all(a < 6 for a in t.atoms)]
    cch = [t.k for t in terms if len(t.atoms) == 3 and any(a >= 6 for a in t.atoms)]
    assert len(cc) == 6 and len(ch) == 6 and len(ccc) == 6 and len(cch) == 12
    assert spread(cc) <= 1e-8
    assert spread(ch) <= 1e-8
    assert spread(ccc) <= 1e-8
    assert spread(cch) <= 1e-8


def test_atom_reordering_invariance() -> None:
    h, bonds, angles = load_case("ethanol_gfn2")
    terms = m2seminario(h, bonds, angles)
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
    bonds2 = [tuple(int(inverse[a]) for a in b) for b in bonds]
    angles2 = [
        (int(inverse[i]), int(inverse[j]), int(inverse[k_at]))
        for i, j, k_at in angles
    ]
    terms2 = m2seminario(h2, bonds2, angles2)
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


def test_oracle_sobtop_m2seminario() -> None:
    """用例 1：除 sf6_fchk 外全部 case 与 param_m3.txt 一致（§2.8 容差，含平衡值）。"""
    failures = {}
    for case in GATE_CASES:
        h, bonds, angles = load_case(case)
        bad = oracle_mismatches(
            m2seminario(h, bonds, angles), read_sobtop_params(REFS / case / "param_m3.txt")
        )
        if bad:
            failures[case] = bad[:5]
    assert not failures, failures


def test_known_limitation_cyclopropane_gfnff_matches_sobtop() -> None:
    """决策 0014：GFN-FF 环丙烷的环角 k 为负，与 Sobtop 在同一 Hessian 上的结果一致。"""
    from sobtop_oracle import KNOWN_LIMITATIONS

    h, bonds, angles = load_case("cyclopropane_gfnff", KNOWN_LIMITATIONS)
    terms = m2seminario(h, bonds, angles)
    params = read_sobtop_params(KNOWN_LIMITATIONS / "cyclopropane_gfnff" / "param_m3.txt")
    assert not oracle_mismatches(terms, params)
    ring = [t.k for t in terms if len(t.atoms) == 3 and set(t.atoms) <= {0, 1, 2}]
    assert len(ring) == 3 and all(k < 0 for k in ring)
