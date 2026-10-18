"""对称等价类与对称化（MVP-1 计划 T20）。

原子等价类 = 带颜色键图的自同构轨道（颜色 = (atomic_number, atom_type,
formal_charge)；只用图，不用 3D 几何）。项的等价类在原子轨道之上再要求
存在一个图自同构把一个项映到另一个（不能只比较类标签，计划 §2.7）。

不枚举整个自同构群（两个叔丁基的分子 |Aut| 就有 ~10⁷）。做法：
1. 颜色细化（1-WL，等价于 RDKit CanonicalRankAtoms(breakTies=False)）得到
   候选划分；它是自同构不变的，所以轨道一定落在同一个细化类里。
2. "钉住若干原子后是否存在自同构" 用个体化-细化回溯搜索判定
   （两份图联合细化，标签可直接比较）。
3. 沿稳定化子链（逐个钉住基点）求一组强生成元：每层只对细化类中
   尚未被已有生成元轨道覆盖的候选做一次搜索。
原子轨道 / 项轨道 = 生成元作用下的连通分量（项轨道在全体图项上做闭包，
所以传入的项集合不必在自同构下封闭）。
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Hashable, Sequence

from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.ir import ForceFieldIR


def _colored_graph(ir: ForceFieldIR) -> tuple[list[int], list[tuple[int, ...]]]:
    """带颜色键图：返回 (每个原子的颜色号, 邻接表)。颜色 = (Z, atom_type, formal_charge)。"""
    n = len(ir.atoms)
    color_index: dict[tuple[int, str, int], int] = {}
    colors: list[int] = []
    for atom in ir.atoms:
        color = (int(atom.atomic_number), str(atom.atom_type), int(atom.formal_charge))
        colors.append(color_index.setdefault(color, len(color_index)))
    neighbors: list[set[int]] = [set() for _ in range(n)]
    for bond in ir.bonds:
        i, j = int(bond.i), int(bond.j)
        if i == j or not (0 <= i < n and 0 <= j < n) or j in neighbors[i]:
            raise ValueError(f"非法键 ({i}, {j})")
        neighbors[i].add(j)
        neighbors[j].add(i)
    return colors, [tuple(sorted(nb)) for nb in neighbors]


def _refine(adj: Sequence[Sequence[int]], init: Sequence[Hashable]) -> list[int]:
    """颜色细化到稳定（等价划分）。标签由排序后的签名决定，与原子编号无关，
    所以对两份图的不交并细化后，两侧标签可以直接比较。"""

    def relabel(keys: Sequence) -> list[int]:
        table = {key: idx for idx, key in enumerate(sorted(set(keys)))}
        return [table[key] for key in keys]

    labels = relabel(list(init))
    n_classes = len(set(labels))
    while True:
        keys = [
            (labels[x], tuple(sorted(labels[y] for y in adj[x])))
            for x in range(len(adj))
        ]
        new = relabel(keys)
        n_new = len(set(new))
        if n_new == n_classes:
            return labels
        labels, n_classes = new, n_new


class _Graph:
    """带颜色图 + 钉住原子的自同构存在性判定。"""

    def __init__(self, colors: list[int], adj: list[tuple[int, ...]]) -> None:
        self.n = len(colors)
        self.colors = colors
        self.adj = adj
        self.edges = {(i, j) for i in range(self.n) for j in adj[i]}
        # 两份图的不交并：右侧原子编号 + n
        n = self.n
        self._union_adj = list(adj) + [tuple(y + n for y in nb) for nb in adj]

    def find_automorphism(
        self, pins: Sequence[tuple[int, int]]
    ) -> tuple[int, ...] | None:
        """找一个颜色保持的自同构 σ，满足每个 (u, v) ∈ pins 都有 σ(u) = v；没有则 None。"""
        n = self.n
        mapping: dict[int, int] = {}
        for u, v in pins:
            if mapping.setdefault(u, v) != v or self.colors[u] != self.colors[v]:
                return None
        if len(set(mapping.values())) != len(mapping):
            return None  # 两个原子钉到同一个像上，不是双射
        pin_left = {u: tag for tag, u in enumerate(mapping)}
        pin_right = {v: tag for tag, v in enumerate(mapping.values())}
        init = [(self.colors[x], pin_left.get(x, -1)) for x in range(n)] + [
            (self.colors[x], pin_right.get(x, -1)) for x in range(n)
        ]
        return self._search(_refine(self._union_adj, init))

    def _search(self, labels: list[int]) -> tuple[int, ...] | None:
        n = self.n
        left, right = labels[:n], labels[n:]
        if Counter(left) != Counter(right):
            return None
        cells_l: dict[int, list[int]] = defaultdict(list)
        cells_r: dict[int, list[int]] = defaultdict(list)
        for x in range(n):
            cells_l[left[x]].append(x)
            cells_r[right[x]].append(x)
        open_cells = [lab for lab, members in cells_l.items() if len(members) > 1]
        if not open_cells:
            sigma = tuple(cells_r[left[x]][0] for x in range(n))
            if all((sigma[i], sigma[j]) in self.edges for i, j in self.edges):
                return sigma
            return None
        # 个体化：取最小的非单元素类，左侧固定一个原子，依次试右侧候选
        target = min(open_cells, key=lambda lab: (len(cells_l[lab]), lab))
        u = cells_l[target][0]
        fresh = max(labels) + 1
        for v in cells_r[target]:
            init = list(labels)
            init[u] = fresh
            init[n + v] = fresh
            found = self._search(_refine(self._union_adj, init))
            if found is not None:
                return found
        return None


class _UnionFind:
    def __init__(self, items: Sequence[Hashable] = ()) -> None:
        self._parent: dict = {item: item for item in items}

    def add(self, x: Hashable) -> bool:
        if x in self._parent:
            return False
        self._parent[x] = x
        return True

    def find(self, x: Hashable) -> Hashable:
        root = x
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[x] != root:
            self._parent[x], x = root, self._parent[x]
        return root

    def union(self, a: Hashable, b: Hashable) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self._parent[rb] = ra


def _generators(ir: ForceFieldIR) -> list[tuple[int, ...]]:
    """颜色保持自同构群的一组（强）生成元，不枚举群元素。

    稳定化子链：依次钉住基点 b1, b2, …；第 i 层在"钉住 b1..b_{i-1} 后的细化类"
    中逐个检查 b_i 的候选像 v：已被本层生成元轨道覆盖就跳过，否则搜一个
    σ（固定 b1..b_{i-1}，σ(b_i)=v），找到即加入。各层轨道都被覆盖，故所得
    生成元生成整个群（Schreier–Sims 论证）。钉住后细化为离散划分时停止。
    """
    colors, adj = _colored_graph(ir)
    graph = _Graph(colors, adj)
    n = graph.n
    gens: list[tuple[int, ...]] = []
    base: list[int] = []
    while True:
        init = [(colors[x], base.index(x) if x in base else -1) for x in range(n)]
        labels = _refine(adj, init)
        cells: dict[int, list[int]] = defaultdict(list)
        for x in range(n):
            cells[labels[x]].append(x)
        open_cells = [members for members in cells.values() if len(members) > 1]
        if not open_cells:
            return gens
        cell = min(open_cells, key=lambda m: (len(m), m[0]))
        b = cell[0]
        orbit = _UnionFind(range(n))
        fixed = [(x, x) for x in base]
        for v in cell[1:]:
            if orbit.find(v) == orbit.find(b):
                continue
            sigma = graph.find_automorphism(fixed + [(b, v)])
            if sigma is None:
                continue
            gens.append(sigma)
            for x in range(n):
                orbit.union(x, sigma[x])
        base.append(b)


def _orbit_labels(items: Sequence[Hashable], uf: _UnionFind) -> tuple[int, ...]:
    labels: dict = {}
    classes = []
    for item in items:
        root = uf.find(item)
        if root not in labels:
            labels[root] = len(labels)
        classes.append(labels[root])
    return tuple(classes)


def atom_classes(ir: ForceFieldIR) -> tuple[int, ...]:
    """每个原子的类号（自同构轨道，按首次出现顺序编号 0..k−1）。"""
    n = len(ir.atoms)
    uf = _UnionFind(range(n))
    for sigma in _generators(ir):
        for i in range(n):
            uf.union(i, sigma[i])
    return _orbit_labels(range(n), uf)


def _normalize_term(term: tuple[int, ...]) -> tuple[int, ...]:
    """bond 无序端点；angle 保持顶点在中间、两端升序；
    dihedral (i,j,k,l) 与 (l,k,j,i) 取字典序小者（MVP-2 计划 T25）。"""
    if len(term) == 2:
        i, j = term
        return (min(i, j), max(i, j))
    if len(term) == 3:
        i, j, k = term
        return (i, j, k) if i <= k else (k, j, i)
    if len(term) == 4:
        return min(term, term[::-1])
    raise ValueError(f"term 只支持 bond/angle/dihedral，得到 {term!r}")


def term_classes(
    ir: ForceFieldIR, terms: Sequence[tuple[int, ...]]
) -> tuple[int, ...]:
    """bond (i,j)、angle (i,j,k)（j 为顶点）或 dihedral (i,j,k,l) 的类号。

    两个项同类 ⟺ 存在颜色保持的图自同构把一个映到另一个
    （bond 端点无序；angle 顶点固定、两端无序）。
    实现：在生成元作用下对规范化后的项做轨道闭包（可能经过不在 terms
    里的中间项），同一轨道即同类。
    """
    normalized = [_normalize_term(tuple(int(a) for a in t)) for t in terms]
    gens = _generators(ir)
    uf = _UnionFind()
    for term in normalized:
        if not uf.add(term):
            continue  # 已在某个先前项的轨道里
        stack = [term]
        while stack:
            current = stack.pop()
            for sigma in gens:
                image = _normalize_term(tuple(sigma[a] for a in current))
                if uf.add(image):
                    stack.append(image)
                uf.union(current, image)
    return _orbit_labels(normalized, uf)


def symmetrize(
    terms: Sequence[DerivedTerm], classes: Sequence[int]
) -> tuple[list[DerivedTerm], dict[int, float]]:
    """类内 k 与平衡值取算术平均（决策 D5，计划 §2.7）。

    返回按输入顺序排列的对称化项，以及每类 k 的相对离散度
    (max−min)/|mean|（mean 为 0 且 max≠min 时记 inf）。
    """
    if len(terms) != len(classes):
        raise ValueError("terms 与 classes 长度不一致")
    groups: dict[int, list[int]] = defaultdict(list)
    for idx, cls in enumerate(classes):
        groups[cls].append(idx)

    spread: dict[int, float] = {}
    averaged: dict[int, tuple[float, float]] = {}
    for cls, members in groups.items():
        ks = [terms[i].k for i in members]
        values = [terms[i].value for i in members]
        k_mean = sum(ks) / len(ks)
        value_mean = sum(values) / len(values)
        averaged[cls] = (value_mean, k_mean)
        denom = abs(k_mean)
        if denom > 0.0:
            spread[cls] = (max(ks) - min(ks)) / denom
        else:
            spread[cls] = 0.0 if max(ks) == min(ks) else math.inf

    result: list[DerivedTerm] = []
    for term, cls in zip(terms, classes):
        value_mean, k_mean = averaged[cls]
        result.append(DerivedTerm(atoms=tuple(term.atoms), value=value_mean, k=k_mean))
    return result, spread
