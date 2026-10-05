"""AMBER 导出器（MVP-0 计划 T7）：原生写 prmtop + inpcrd（决策 0015，不依赖 ParmEd）。

½k 约定：IR 按 E = ½kΔq² 储存，AMBER 是 E = KΔq²，导出时 K = k/2（再换 kcal、Å）。
1-4 规则：``pairs14`` 中每个对恰好挂在一个"计算 1-4"的 dihedral 条目上；同一对出现在
多个二面角时其余条目 ignore_end；没有任何 proper 项的 1-4 对补一条 k=0、n=1、δ=0 的
载体条目。不写 frcmod（路线图 §3.4）。

数组布局与 tleap/ParmEd 相同：含氢（任一原子为 H）与不含氢的项分开；二面角第 3 个
原子取负表示不计 1-4、第 4 个取负表示 improper；原子 0 落在第 3/4 位时整体反序
（φ 在反序下不变）。参数类型按取值去重。LJ 按 atom_type 分类型，A/B 系数用
Lorentz–Berthelot 组合。
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ffcompiler.formats.amber import POINTER_NAMES, PrmtopFile, write_prmtop, write_rst7
from ffcompiler.ir import ForceFieldIR, validate

KCAL_KJ = 4.184
_K_BOND = 2.0 * KCAL_KJ * 100.0   # IR kJ/mol/nm²（½k）→ AMBER kcal/mol/Å²（无 ½）
_K_ANGLE = 2.0 * KCAL_KJ          # IR kJ/mol/rad²（½k）→ AMBER kcal/mol/rad²
_K_TORSION = KCAL_KJ
_CHARGE_SCALE = 18.2223
_VERSION = "%VERSION  VERSION_STAMP = V0001.000  WRITER = ffcompiler.exporters.amber"


def _pair(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


class _TypeTable:
    """按参数取值去重的类型表（1-based 下标）。"""

    def __init__(self) -> None:
        self.index: dict[tuple, int] = {}
        self.params: list[tuple] = []

    def get(self, params: tuple) -> int:
        if params not in self.index:
            self.params.append(params)
            self.index[params] = len(self.params)
        return self.index[params]

    def column(self, k: int) -> list:
        return [p[k] for p in self.params]


def _carrier_quartets(ir: ForceFieldIR) -> list[tuple[int, int, int, int]]:
    """没有 proper 项的 1-4 对：沿真实键路径补 k=0 载体四元组。"""
    adjacency: dict[int, list[int]] = {}
    for bond in ir.bonds:
        adjacency.setdefault(bond.i, []).append(bond.j)
        adjacency.setdefault(bond.j, []).append(bond.i)
    with_proper = {_pair(t.i, t.l) for t in ir.proper_torsions}
    quartets = []
    for i, j in sorted(_pair(p[0], p[1]) for p in ir.pairs14):
        if (i, j) in with_proper:
            continue
        path = next(
            ((i, n1, n2, j) for n1 in adjacency[i] for n2 in adjacency[n1]
             if n2 != i and j in adjacency[n2] and j != n1),
            None,
        )
        if path is None:
            raise ValueError(f"1-4 pair ({i},{j}) has no 3-bond path in the bond graph")
        quartets.append(path)
    return quartets


def build_prmtop(ir: ForceFieldIR) -> PrmtopFile:
    validate(ir)
    natom = len(ir.atoms)
    is_h = [a.atomic_number == 1 for a in ir.atoms]

    # ---- LJ 类型（按 atom_type，validate 保证同类型 LJ 相同）-------------------
    lj_types: dict[str, int] = {}
    half_rmin: list[float] = []
    eps: list[float] = []
    for atom in ir.atoms:
        if atom.atom_type not in lj_types:
            lj_types[atom.atom_type] = len(lj_types) + 1
            half_rmin.append((atom.sigma / 0.1) * 2.0 ** (1.0 / 6.0) / 2.0)
            eps.append(atom.epsilon / KCAL_KJ)
    ntypes = len(lj_types)
    ico = [0] * (ntypes * ntypes)
    n_pairs = ntypes * (ntypes + 1) // 2
    acoef = [0.0] * n_pairs
    bcoef = [0.0] * n_pairs
    for t1 in range(ntypes):
        for t2 in range(ntypes):
            lo, hi = min(t1, t2) + 1, max(t1, t2) + 1
            idx = hi * (hi - 1) // 2 + lo
            ico[ntypes * t1 + t2] = idx
            rmin = half_rmin[t1] + half_rmin[t2]
            e = math.sqrt(eps[t1] * eps[t2])
            acoef[idx - 1] = e * rmin**12
            bcoef[idx - 1] = 2.0 * e * rmin**6

    # ---- bond / angle --------------------------------------------------------
    bond_types, angle_types, dihedral_types = _TypeTable(), _TypeTable(), _TypeTable()
    bonds_h, bonds_a = [], []
    for bond in ir.bonds:
        t = bond_types.get((bond.k / _K_BOND, bond.r0 * 10.0))
        (bonds_h if is_h[bond.i] or is_h[bond.j] else bonds_a).extend(
            [bond.i * 3, bond.j * 3, t]
        )
    angles_h, angles_a = [], []
    for angle in ir.angles:
        t = angle_types.get((angle.force_constant / _K_ANGLE, angle.theta0))
        members = (angle.i, angle.j, angle.k_atom)
        (angles_h if any(is_h[m] for m in members) else angles_a).extend(
            [angle.i * 3, angle.j * 3, angle.k_atom * 3, t]
        )

    # ---- dihedral（1-4 归属：每对只允许一个"计算 1-4"的条目）------------------
    scee, scnb = 1.0 / ir.coulomb_14_scale, 1.0 / ir.vdw_14_scale
    pairs14 = {_pair(p[0], p[1]) for p in ir.pairs14}
    first_carrier: dict[tuple[int, int], int] = {}
    for n, term in enumerate(ir.proper_torsions):
        first_carrier.setdefault(_pair(term.i, term.l), n)
    entries: list[tuple[tuple[int, int, int, int], tuple, bool, bool]] = []
    for n, term in enumerate(ir.proper_torsions):
        pair = _pair(term.i, term.l)
        carries = pair in pairs14 and first_carrier[pair] == n
        params = (term.force_constant / _K_TORSION, float(term.periodicity), term.phase, scee, scnb)
        entries.append(((term.i, term.j, term.k, term.l), params, not carries, False))
    for term in ir.impropers:
        params = (term.force_constant / _K_TORSION, float(term.periodicity), term.phase, scee, scnb)
        entries.append(((term.i, term.j, term.k, term.l), params, True, True))
    for quartet in _carrier_quartets(ir):
        entries.append((quartet, (0.0, 1.0, 0.0, scee, scnb), False, False))

    dihed_h, dihed_a = [], []
    for quartet, params, ignore_end, improper in entries:
        t = dihedral_types.get(params)
        a, b, c, d = quartet
        if c == 0 or d == 0:
            a, b, c, d = d, c, b, a
        row = [a * 3, b * 3, c * 3 * (-1 if ignore_end else 1),
               d * 3 * (-1 if improper else 1), t]
        (dihed_h if any(is_h[m] for m in quartet) else dihed_a).extend(row)

    # ---- 排除表：1-2 ∪ 1-3 ∪ 1-4，前向引用（1-based），没有则 0 哨兵 ---------------
    forward: dict[int, list[int]] = {i: [] for i in range(natom)}
    for i, j in {_pair(p[0], p[1]) for p in (*ir.exclusions, *ir.pairs14)}:
        forward[i].append(j)
    n_excluded, excluded = [], []
    for i in range(natom):
        js = sorted(forward[i])
        n_excluded.append(max(len(js), 1))
        excluded.extend([j + 1 for j in js] if js else [0])

    pointers = dict.fromkeys(POINTER_NAMES[:31], 0)
    pointers.update(
        NATOM=natom, NTYPES=ntypes, NBONH=len(bonds_h) // 3, MBONA=len(bonds_a) // 3,
        NTHETH=len(angles_h) // 4, MTHETA=len(angles_a) // 4,
        NPHIH=len(dihed_h) // 5, MPHIA=len(dihed_a) // 5, NNB=len(excluded), NRES=1,
        NBONA=len(bonds_a) // 3, NTHETA=len(angles_a) // 4, NPHIA=len(dihed_a) // 5,
        NUMBND=len(bond_types.params), NUMANG=len(angle_types.params),
        NPTRA=len(dihedral_types.params), NATYP=ntypes, NMXRS=natom,
    )

    prm = PrmtopFile(version=_VERSION)
    e, i8, a4 = "5E16.8", "10I8", "20a4"
    title = ir.name[:80]
    prm.set("TITLE", a4, [title[k:k + 4] for k in range(0, len(title), 4)])
    prm.set("POINTERS", i8, [pointers[name] for name in POINTER_NAMES[:31]])
    prm.set("ATOM_NAME", a4, [a.name for a in ir.atoms])
    prm.set("CHARGE", e, [a.partial_charge * _CHARGE_SCALE for a in ir.atoms])
    prm.set("ATOMIC_NUMBER", i8, [a.atomic_number for a in ir.atoms])
    prm.set("MASS", e, [a.mass for a in ir.atoms])
    prm.set("ATOM_TYPE_INDEX", i8, [lj_types[a.atom_type] for a in ir.atoms])
    prm.set("NUMBER_EXCLUDED_ATOMS", i8, n_excluded)
    prm.set("NONBONDED_PARM_INDEX", i8, ico)
    prm.set("RESIDUE_LABEL", a4, [ir.atoms[0].residue_name])
    prm.set("RESIDUE_POINTER", i8, [1])
    prm.set("BOND_FORCE_CONSTANT", e, bond_types.column(0))
    prm.set("BOND_EQUIL_VALUE", e, bond_types.column(1))
    prm.set("ANGLE_FORCE_CONSTANT", e, angle_types.column(0))
    prm.set("ANGLE_EQUIL_VALUE", e, angle_types.column(1))
    prm.set("DIHEDRAL_FORCE_CONSTANT", e, dihedral_types.column(0))
    prm.set("DIHEDRAL_PERIODICITY", e, dihedral_types.column(1))
    prm.set("DIHEDRAL_PHASE", e, dihedral_types.column(2))
    prm.set("SCEE_SCALE_FACTOR", e, dihedral_types.column(3))
    prm.set("SCNB_SCALE_FACTOR", e, dihedral_types.column(4))
    prm.set("SOLTY", e, [0.0] * ntypes)
    prm.set("LENNARD_JONES_ACOEF", e, acoef)
    prm.set("LENNARD_JONES_BCOEF", e, bcoef)
    prm.set("BONDS_INC_HYDROGEN", i8, bonds_h)
    prm.set("BONDS_WITHOUT_HYDROGEN", i8, bonds_a)
    prm.set("ANGLES_INC_HYDROGEN", i8, angles_h)
    prm.set("ANGLES_WITHOUT_HYDROGEN", i8, angles_a)
    prm.set("DIHEDRALS_INC_HYDROGEN", i8, dihed_h)
    prm.set("DIHEDRALS_WITHOUT_HYDROGEN", i8, dihed_a)
    prm.set("EXCLUDED_ATOMS_LIST", i8, excluded)
    prm.set("HBOND_ACOEF", e, [])
    prm.set("HBOND_BCOEF", e, [])
    prm.set("HBCUT", e, [])
    prm.set("AMBER_ATOM_TYPE", a4, [a.atom_type for a in ir.atoms])
    prm.set("TREE_CHAIN_CLASSIFICATION", a4, ["BLA"] * natom)
    prm.set("JOIN_ARRAY", i8, [0] * natom)
    prm.set("IROTAT", i8, [0] * natom)
    prm.set("RADIUS_SET", "1a80", ["none (gas-phase export, radii not set)"])
    prm.set("RADII", e, [0.0] * natom)
    prm.set("SCREEN", e, [0.0] * natom)
    prm.set("IPOL", "1I8", [0])
    return prm


def write_amber(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    outdir: Path,
) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prmtop = outdir / "molecule.prmtop"
    inpcrd = outdir / "molecule.inpcrd"
    write_prmtop(build_prmtop(ir), prmtop)
    write_rst7(np.asarray(positions_nm, dtype=np.float64) * 10.0, inpcrd, title=ir.name)
    return {"prmtop": prmtop, "inpcrd": inpcrd}
