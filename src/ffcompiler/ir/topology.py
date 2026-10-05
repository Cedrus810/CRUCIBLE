"""键图与按最短路径键数分类的原子对（MVP-0 计划 T2）。

约定（计划 §2.1）：``exclusions`` = 拓扑 1-2 ∪ 1-3；``pairs14`` = 最短路径
恰为 3 条键的原子对（AMBER 家族约定）。一个对不能同时出现在两处。
"""

from __future__ import annotations

from collections import deque

from .schema import ForceFieldIR


def bond_graph(ir: ForceFieldIR) -> dict[int, frozenset[int]]:
    """邻接表，只由 ``ir.bonds`` 决定。

    对索引不在 ``ir.atoms`` 里的端点也容忍（setdefault），校验器据此统一报问题。
    """
    graph: dict[int, set[int]] = {atom.index: set() for atom in ir.atoms}
    for bond in ir.bonds:
        graph.setdefault(bond.i, set()).add(bond.j)
        graph.setdefault(bond.j, set()).add(bond.i)
    return {idx: frozenset(nbrs) for idx, nbrs in graph.items()}


def _bfs_distances(graph: dict[int, frozenset[int]], src: int) -> dict[int, int]:
    dist = {src: 0}
    queue = deque([src])
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if nbr not in dist:
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def topological_pairs(
    ir: ForceFieldIR,
) -> tuple[frozenset[tuple[int, int]], frozenset[tuple[int, int]], frozenset[tuple[int, int]]]:
    """按最短路径键数分出 (1-2, 1-3, 1-4) 三个集合，元素 i<j，三者互不相交。"""
    graph = bond_graph(ir)
    pairs12: set[tuple[int, int]] = set()
    pairs13: set[tuple[int, int]] = set()
    pairs14: set[tuple[int, int]] = set()
    for src in sorted(graph):
        for dst, d in _bfs_distances(graph, src).items():
            if dst <= src:
                continue
            pair = (src, dst)
            if d == 1:
                pairs12.add(pair)
            elif d == 2:
                pairs13.add(pair)
            elif d == 3:
                pairs14.add(pair)
    return frozenset(pairs12), frozenset(pairs13), frozenset(pairs14)


def expected_exclusions_and_pairs(
    ir: ForceFieldIR,
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """按 §2.1 约定生成 (exclusions, pairs14)，各自升序。"""
    pairs12, pairs13, pairs14 = topological_pairs(ir)
    exclusions = tuple(sorted(pairs12 | pairs13))
    pairs = tuple(sorted(pairs14))
    return exclusions, pairs
