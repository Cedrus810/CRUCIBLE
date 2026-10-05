"""T20：对称等价类的必测用例。"""

from __future__ import annotations

import math
import time
from itertools import combinations
from pathlib import Path

import pytest
from rdkit import Chem

from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.bonded.symmetry import atom_classes, symmetrize, term_classes
from ffcompiler.importers import from_amber
from ffcompiler.ir import Atom, Bond, ForceFieldIR
from ffcompiler.ir.topology import bond_graph

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _load(name: str):
    mdir = FIXTURES / name
    return from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )


def _class_groups(classes) -> dict[int, list[int]]:
    groups: dict[int, list[int]] = {}
    for idx, cls in enumerate(classes):
        groups.setdefault(cls, []).append(idx)
    return groups


def test_benzene_classes() -> None:
    ir, _ = _load("benzene")
    atom_cls = atom_classes(ir)
    carbon = [a.index for a in ir.atoms if a.element == "C"]
    hydrogen = [a.index for a in ir.atoms if a.element == "H"]
    assert len({atom_cls[i] for i in carbon}) == 1
    assert len({atom_cls[i] for i in hydrogen}) == 1
    assert atom_cls[carbon[0]] != atom_cls[hydrogen[0]]

    bonds = [(b.i, b.j) for b in ir.bonds]
    bond_cls = term_classes(ir, bonds)
    cc = [(b.i, b.j) for b in ir.bonds if {ir.atoms[b.i].element, ir.atoms[b.j].element} == {"C"}]
    ch = [(b.i, b.j) for b in ir.bonds if {ir.atoms[b.i].element, ir.atoms[b.j].element} == {"C", "H"}]
    assert len(cc) == 6 and len(ch) == 6
    assert len({bond_cls[bonds.index(t)] for t in cc}) == 1
    assert len({bond_cls[bonds.index(t)] for t in ch}) == 1

    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]
    angle_cls = term_classes(ir, angles)
    ccc = [t for t in angles if ir.atoms[t[0]].element == "C" and ir.atoms[t[2]].element == "C"]
    cch = [t for t in angles if {ir.atoms[t[0]].element, ir.atoms[t[2]].element} == {"C", "H"}]
    assert len(ccc) == 6 and len(cch) == 12
    assert len({angle_cls[angles.index(t)] for t in ccc}) == 1
    assert len({angle_cls[angles.index(t)] for t in cch}) == 1


def test_ethanol_hydrogen_classes() -> None:
    ir, _ = _load("ethanol")
    atom_cls = atom_classes(ir)
    groups = _class_groups(atom_cls)
    graph = bond_graph(ir)
    hydrogens = [a.index for a in ir.atoms if a.element == "H"]
    assert len(hydrogens) == 6  # CH3 ×3 + CH2 ×2 + OH ×1
    h_classes = [atom_cls[i] for i in hydrogens]
    assert len(set(h_classes)) == 3  # 三种氢

    def h_on(heavy: int) -> list[int]:
        return sorted(j for j in graph[heavy] if ir.atoms[j].element == "H")

    oxygen = next(a.index for a in ir.atoms if a.element == "O")
    ch2 = next(j for j in graph[oxygen] if ir.atoms[j].element == "C")
    ch3 = next(
        a.index for a in ir.atoms if a.element == "C" and a.index != ch2
    )
    h_ch3, h_ch2, h_oh = h_on(ch3), h_on(ch2), h_on(oxygen)
    assert (len(h_ch3), len(h_ch2), len(h_oh)) == (3, 2, 1)
    # CH3 的 3 个 H 恰为一类，CH2 的 2 个 H 恰为一类，OH 的 H 单独一类
    assert sorted(groups[atom_cls[h_ch3[0]]]) == h_ch3
    assert sorted(groups[atom_cls[h_ch2[0]]]) == h_ch2
    assert groups[atom_cls[h_oh[0]]] == h_oh


def test_norbornane_bridgeheads() -> None:
    ir, _ = _load("norbornane")
    atom_cls = atom_classes(ir)
    graph = bond_graph(ir)
    # 桥头碳 = 有 3 个碳邻居的碳（全图含 H，sp³ 碳度都是 4）；降冰片烷恰好 2 个
    bridgeheads = [
        i for i, nbrs in graph.items()
        if ir.atoms[i].element == "C"
        and sum(1 for j in nbrs if ir.atoms[j].element == "C") == 3
    ]
    assert len(bridgeheads) == 2
    assert atom_cls[bridgeheads[0]] == atom_cls[bridgeheads[1]]


def _pentane_ir() -> ForceFieldIR:
    atoms = tuple(
        Atom(
            index=i, name=f"C{i}", element="C", atomic_number=6, mass=12.011,
            atom_type="C", formal_charge=0, partial_charge=0.0,
            sigma=0.3, epsilon=0.4,
        )
        for i in range(5)
    )
    bonds = tuple(
        Bond(i=i, j=i + 1, r0=0.154, k=2000.0, source="t") for i in range(4)
    )
    return ForceFieldIR(
        name="pentane-skeleton", atoms=atoms, bonds=bonds, angles=(),
        proper_torsions=(), impropers=(), exclusions=(), pairs14=(),
        total_charge=0, forcefield_family="amber",
        combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )


def test_pentane_same_color_but_distinct_classes() -> None:
    ir = _pentane_ir()
    atom_cls = atom_classes(ir)
    assert atom_cls[0] == atom_cls[4]  # C1/C5
    assert atom_cls[1] == atom_cls[3]  # C2/C4
    assert len({atom_cls[0], atom_cls[1], atom_cls[2]}) == 3  # C3 单独

    bond_cls = term_classes(ir, [(b.i, b.j) for b in ir.bonds])
    # C1-C2 与 C4-C5 同类；C2-C3 单独
    assert bond_cls[0] == bond_cls[3]
    assert len({bond_cls[0], bond_cls[1], bond_cls[2]}) == 2


def test_symmetrize_averaging_and_spread() -> None:
    terms = [
        DerivedTerm(atoms=(0, 1), value=0.15, k=1.0),
        DerivedTerm(atoms=(1, 2), value=0.15, k=2.0),
        DerivedTerm(atoms=(2, 3), value=0.15, k=3.0),
        DerivedTerm(atoms=(4, 5), value=0.10, k=7.0),
    ]
    classes = (0, 0, 0, 1)
    result, spread = symmetrize(terms, classes)
    assert all(t.k == pytest.approx(2.0) for t in result[:3])
    assert all(t.value == pytest.approx(0.15) for t in result[:3])
    assert result[3].k == pytest.approx(7.0)
    assert spread[0] == pytest.approx(1.0)  # (3−1)/2
    assert spread[1] == 0.0
    # 输入顺序保持
    assert [t.atoms for t in result] == [t.atoms for t in terms]


def test_symmetrize_zero_mean_spread() -> None:
    terms = [
        DerivedTerm(atoms=(0, 1), value=0.1, k=1.0),
        DerivedTerm(atoms=(1, 2), value=0.1, k=-1.0),
    ]
    _, spread = symmetrize(terms, (0, 0))
    assert math.isinf(spread[0])


# ---------------------------------------------------------------------------
# 高对称分子（自同构群远超 1e5）：不能靠枚举全部自同构
# ---------------------------------------------------------------------------


def _ir_from_smiles(smiles: str) -> ForceFieldIR:
    """SMILES（加全氢）→ 只含原子与键的 IR；atom_type 取元素符号。"""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    atoms = tuple(
        Atom(
            index=a.GetIdx(), name=f"{a.GetSymbol()}{a.GetIdx()}",
            element=a.GetSymbol(), atomic_number=a.GetAtomicNum(),
            mass=a.GetMass(), atom_type=a.GetSymbol(),
            formal_charge=a.GetFormalCharge(), partial_charge=0.0,
            sigma=0.3, epsilon=0.4,
        )
        for a in mol.GetAtoms()
    )
    bonds = tuple(
        Bond(i=b.GetBeginAtomIdx(), j=b.GetEndAtomIdx(), r0=0.15, k=1000.0, source="t")
        for b in mol.GetBonds()
    )
    return ForceFieldIR(
        name=smiles, atoms=atoms, bonds=bonds, angles=(),
        proper_torsions=(), impropers=(), exclusions=(), pairs14=(),
        total_charge=sum(a.formal_charge for a in atoms), forcefield_family="amber",
        combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )


def _all_angles(ir: ForceFieldIR) -> list[tuple[int, int, int]]:
    graph = bond_graph(ir)
    return [
        (i, j, k)
        for j in range(len(ir.atoms))
        for i, k in combinations(sorted(graph[j]), 2)
    ]


def _partition(items, classes) -> set[frozenset]:
    groups: dict[int, set] = {}
    for item, cls in zip(items, classes):
        groups.setdefault(cls, set()).add(item)
    return {frozenset(g) for g in groups.values()}


def _atom_kinds(ir: ForceFieldIR) -> dict[str, list[int]]:
    """按 (元素, 重原子邻居数, 氢数) 粗分原子，便于写期望。"""
    graph = bond_graph(ir)
    kinds: dict[str, list[int]] = {}
    for a in ir.atoms:
        heavy = sum(1 for j in graph[a.index] if ir.atoms[j].element != "H")
        nh = sum(1 for j in graph[a.index] if ir.atoms[j].element == "H")
        if a.element == "H":
            parent = next(iter(graph[a.index]))
            key = "H@" + _kind_of(ir, graph, parent)
        else:
            key = f"{a.element}{heavy}h{nh}"
        kinds.setdefault(key, []).append(a.index)
    return kinds


def _kind_of(ir, graph, i) -> str:
    heavy = sum(1 for j in graph[i] if ir.atoms[j].element != "H")
    nh = sum(1 for j in graph[i] if ir.atoms[j].element == "H")
    return f"{ir.atoms[i].element}{heavy}h{nh}"


def test_tetramethylpentane_two_tert_butyls() -> None:
    # 2,2,4,4-四甲基戊烷：|Aut| = 6^6·6^2·2·2 ≈ 6.7e6，旧实现会撞上 1e5 上限
    ir = _ir_from_smiles("CC(C)(C)CC(C)(C)C")
    kinds = _atom_kinds(ir)
    assert len(kinds["C1h3"]) == 6 and len(kinds["H@C1h3"]) == 18
    assert len(kinds["C4h0"]) == 2 and len(kinds["C2h2"]) == 1
    assert len(kinds["H@C2h2"]) == 2

    atom_cls = atom_classes(ir)
    expected_atoms = {frozenset(v) for v in kinds.values()}
    assert _partition(range(len(ir.atoms)), atom_cls) == expected_atoms

    graph = bond_graph(ir)
    bonds = [(b.i, b.j) for b in ir.bonds]
    bond_cls = term_classes(ir, bonds)
    by_kind: dict[frozenset, set] = {}
    for t in bonds:
        by_kind.setdefault(frozenset(_kind_of(ir, graph, a) if ir.atoms[a].element != "H" else "H" for a in t), set()).add(t)
    assert _partition(bonds, bond_cls) == {frozenset(v) for v in by_kind.values()}
    assert len(set(bond_cls)) == 4  # CH3-H, CH2-H, Cq-CH3, Cq-CH2

    angles = _all_angles(ir)
    angle_cls = term_classes(ir, angles)
    # 7 类：H-C-H(甲基)、Cq-C-H(甲基)、CH3-Cq-CH3、CH3-Cq-CH2、Cq-CH2-Cq、H-CH2-H、Cq-CH2-H
    assert len(set(angle_cls)) == 7
    methyl_hch = [t for t in angles if _kind_of(ir, graph, t[1]) == "C1h3"
                  and ir.atoms[t[0]].element == ir.atoms[t[2]].element == "H"]
    assert len(methyl_hch) == 18
    assert len({angle_cls[angles.index(t)] for t in methyl_hch}) == 1


def test_di_tert_butylbenzene() -> None:
    ir = _ir_from_smiles("CC(C)(C)c1ccc(cc1)C(C)(C)C")
    kinds = _atom_kinds(ir)
    atom_cls = atom_classes(ir)
    assert _partition(range(len(ir.atoms)), atom_cls) == {
        frozenset(v) for v in kinds.values()
    }
    assert len(kinds["H@C1h3"]) == 18 and len(kinds["H@C2h1"]) == 4

    graph = bond_graph(ir)
    bonds = [(b.i, b.j) for b in ir.bonds]
    bond_cls = term_classes(ir, bonds)
    ring_cc = [t for t in bonds if {_kind_of(ir, graph, a) for a in t} == {"C2h1"}]
    ipso_cc = [t for t in bonds if {_kind_of(ir, graph, a) for a in t} == {"C2h1", "C3h0"}]
    assert len(ring_cc) == 2 and len(ipso_cc) == 4
    assert len({bond_cls[bonds.index(t)] for t in ring_cc}) == 1
    assert len({bond_cls[bonds.index(t)] for t in ipso_cc}) == 1
    assert bond_cls[bonds.index(ring_cc[0])] != bond_cls[bonds.index(ipso_cc[0])]
    # 甲基 C-H(18)、芳 C-H(4)、Cq-CH3(6)、Cq-C芳(2)、芳 CH-CH(2)、芳 C-CH(4)
    assert len(set(bond_cls)) == 6


def test_neopentane_fast() -> None:
    # |Aut| = 6^4·4! = 31104：旧实现逐个枚举很慢
    ir = _ir_from_smiles("CC(C)(C)C")
    start = time.perf_counter()
    atom_cls = atom_classes(ir)
    angles = _all_angles(ir)
    angle_cls = term_classes(ir, angles)
    bond_cls = term_classes(ir, [(b.i, b.j) for b in ir.bonds])
    elapsed = time.perf_counter() - start
    assert elapsed < 0.5, f"neopentane 用时 {elapsed:.2f} s"

    kinds = _atom_kinds(ir)
    assert len(kinds["C1h3"]) == 4 and len(kinds["H@C1h3"]) == 12
    assert len({atom_cls[i] for i in kinds["C1h3"]}) == 1
    assert len({atom_cls[i] for i in kinds["H@C1h3"]}) == 1
    assert len(set(atom_cls)) == 3
    assert len(set(bond_cls)) == 2  # C-C、C-H
    assert len(set(angle_cls)) == 3  # C-C-C、C-C-H、H-C-H


def _brute_force_classes(ir: ForceFieldIR, terms) -> tuple[int, ...]:
    """小分子上用 RDKit 枚举全部颜色保持自同构，作为对照。"""
    mol = Chem.RWMol()
    colors: dict = {}
    for a in ir.atoms:
        c = (a.atomic_number, a.atom_type, a.formal_charge)
        atom = Chem.Atom(a.atomic_number)
        atom.SetIsotope(colors.setdefault(c, len(colors) + 1))
        mol.AddAtom(atom)
    for b in ir.bonds:
        mol.AddBond(b.i, b.j, Chem.BondType.SINGLE)
    autos = mol.GetSubstructMatches(mol, uniquify=False, maxMatches=1_000_000)

    def norm(t):
        if len(t) == 1:
            return t
        if len(t) == 2:
            return tuple(sorted(t))
        return t if t[0] <= t[2] else (t[2], t[1], t[0])

    orbit_id: dict[tuple, int] = {}
    out = []
    for t in terms:
        key = min(norm(tuple(s[a] for a in t)) for s in autos)
        out.append(orbit_id.setdefault(key, len(orbit_id)))
    return tuple(out)


@pytest.mark.parametrize(
    "smiles",
    [
        "CCO", "CC(C)C", "Cc1ccccc1", "OC(=O)CC(=O)O", "C1CCC(CC1)C",
        "CC(=O)[O-]", "C[NH3+]",
        "C12C3C4C1C5C2C3C45",  # 立方烷
        "C1C2CC3CC1CC(C2)C3",  # 金刚烷
        "C1CC1.C1CCCCC1",  # 颜色细化（1-WL）分不开两种 CH2，必须靠搜索区分
    ],
)
def test_matches_brute_force(smiles: str) -> None:
    ir = _ir_from_smiles(smiles)
    atoms = [(i,) for i in range(len(ir.atoms))]
    bonds = [(b.i, b.j) for b in ir.bonds]
    angles = _all_angles(ir)
    assert _partition(range(len(atoms)), atom_classes(ir)) == _partition(
        range(len(atoms)), _brute_force_classes(ir, atoms)
    )
    assert _partition(bonds, term_classes(ir, bonds)) == _partition(
        bonds, _brute_force_classes(ir, bonds)
    )
    assert _partition(angles, term_classes(ir, angles)) == _partition(
        angles, _brute_force_classes(ir, angles)
    )


def test_term_subset_still_grouped() -> None:
    # 只传部分项（项集合不在自同构下封闭）也必须按"存在自同构"判同类
    ir = _ir_from_smiles("CC(C)C")
    graph = bond_graph(ir)
    methyls = [i for i in range(len(ir.atoms)) if _kind_of(ir, graph, i) == "C1h3"]
    h_a = min(j for j in graph[methyls[0]] if ir.atoms[j].element == "H")
    h_b = max(j for j in graph[methyls[2]] if ir.atoms[j].element == "H")
    cls = term_classes(ir, [(methyls[0], h_a), (h_b, methyls[2])])
    assert cls[0] == cls[1]


def test_angle_endpoints_unordered_vertex_fixed() -> None:
    ir = _pentane_ir()
    # (0,1,2) 与 (2,3,4)（顶点不同但对称）同类；(1,2,3) 单独；端点顺序无关
    cls = term_classes(ir, [(0, 1, 2), (4, 3, 2), (1, 2, 3), (3, 2, 1)])
    assert cls[0] == cls[1]
    assert cls[2] == cls[3]
    assert cls[0] != cls[2]
