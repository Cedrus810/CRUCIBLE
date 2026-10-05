"""T2：拓扑工具的必测用例（手工构造的 IR，只需要原子和 bond）。"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import replace

from ffcompiler.ir import Atom, Bond, ForceFieldIR
from ffcompiler.ir.topology import (
    bond_graph,
    expected_exclusions_and_pairs,
    topological_pairs,
)

_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "Cl": 17}


def make_ir(elements: list[str], bonds: list[tuple[int, int]], name: str = "t") -> ForceFieldIR:
    atoms = tuple(
        Atom(
            index=i, name=f"{e}{i}", element=e, atomic_number=_Z[e], mass=12.0,
            atom_type=f"T{i}", formal_charge=0, partial_charge=0.0,
            sigma=0.3, epsilon=0.1,
        )
        for i, e in enumerate(elements)
    )
    bond_terms = tuple(Bond(i=i, j=j, r0=0.15, k=3.0e5, source="test") for i, j in bonds)
    skeleton = ForceFieldIR(
        name=name, atoms=atoms, bonds=bond_terms, angles=(), proper_torsions=(),
        impropers=(), exclusions=(), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    exclusions, pairs14 = expected_exclusions_and_pairs(skeleton)
    return replace(skeleton, exclusions=exclusions, pairs14=pairs14)


def _ref_pairs_by_distance(
    n: int, edges: list[tuple[int, int]], want: int
) -> set[tuple[int, int]]:
    """独立的 BFS 参照实现，用于大分子的交叉核对。"""
    adj: dict[int, set[int]] = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    out: set[tuple[int, int]] = set()
    for s in range(n):
        dist = {s: 0}
        queue = deque([s])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        for t, d in dist.items():
            if t > s and d == want:
                out.add((s, t))
    return out


def test_butane_one_14_pair() -> None:
    ir = make_ir(["C"] * 4, [(0, 1), (1, 2), (2, 3)], name="butane")
    pairs12, pairs13, pairs14 = topological_pairs(ir)
    assert pairs12 == frozenset({(0, 1), (1, 2), (2, 3)})
    assert pairs13 == frozenset({(0, 2), (1, 3)})
    assert pairs14 == frozenset({(0, 3)})


def test_cyclopropane_no_13_14() -> None:
    ir = make_ir(["C"] * 3, [(0, 1), (1, 2), (2, 0)], name="cyclopropane")
    pairs12, pairs13, pairs14 = topological_pairs(ir)
    assert pairs12 == frozenset({(0, 1), (0, 2), (1, 2)})
    assert pairs13 == frozenset()
    assert pairs14 == frozenset()


def test_cyclobutane_13_diagonals_no_14() -> None:
    ir = make_ir(["C"] * 4, [(0, 1), (1, 2), (2, 3), (3, 0)], name="cyclobutane")
    pairs12, pairs13, pairs14 = topological_pairs(ir)
    assert pairs13 == frozenset({(0, 2), (1, 3)})
    assert pairs14 == frozenset()


def test_cyclopentane_ring_no_14_but_h_yes() -> None:
    ring = [(i, (i + 1) % 5) for i in range(5)]
    ir = make_ir(["C"] * 5, ring, name="cyclopentane")
    _, _, pairs14 = topological_pairs(ir)
    assert pairs14 == frozenset()  # 环上任意两点最短路径 ≤ 2

    # 带 H：每个 C 两个 H（idx 5+2i, 6+2i）
    elements = ["C"] * 5 + ["H"] * 10
    edges = list(ring)
    for c in range(5):
        edges.append((c, 5 + 2 * c))
        edges.append((c, 6 + 2 * c))
    irh = make_ir(elements, edges, name="cyclopentane_h")
    pairs12h, pairs13h, pairs14h = topological_pairs(irh)
    assert pairs14h == frozenset(_ref_pairs_by_distance(15, edges, 3))
    assert len(pairs14h) == 40  # 10H×2 个 (H,C) + 相邻碳的 H–H 20 个，且每对只出现一次
    assert (2, 5) in pairs14h   # H(C0) 与隔两根键的 C2
    assert (5, 7) in pairs14h   # 相邻碳上的 H–H
    assert not (pairs14h & pairs13h) and not (pairs14h & pairs12h)


def test_norbornane_14_once_and_disjoint() -> None:
    # C1CC2CCC1C2 骨架：8 条键，桥头 C2、C5
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6), (2, 6)]
    ir = make_ir(["C"] * 7, bonds, name="norbornane")
    pairs12, pairs13, pairs14 = topological_pairs(ir)

    assert pairs14 == frozenset({(0, 3), (1, 4)})  # 手工推导：只有这两个距离 3 的对
    assert pairs14 == frozenset(_ref_pairs_by_distance(7, bonds, 3))
    assert len(pairs14) == len(_ref_pairs_by_distance(7, bonds, 3))  # 每对恰好一次
    assert not (pairs14 & pairs13) and not (pairs14 & pairs12)  # 没有对同时落在两类

    exclusions, pairs = expected_exclusions_and_pairs(ir)
    assert exclusions == tuple(sorted(pairs12 | pairs13))
    assert pairs == tuple(sorted(pairs14))
    graph = bond_graph(ir)
    assert graph[2] == frozenset({1, 3, 6})  # 桥头原子
