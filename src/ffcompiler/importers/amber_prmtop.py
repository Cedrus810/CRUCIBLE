"""AMBER prmtop 导入器（MVP-0 计划 T5）。

把 GAFF2 + AM1-BCC 的 prmtop 按计划 §2.2 的换算读进引擎无关的 ForceFieldIR。
导入器只做格式转换，不做任何参数推导。文件解析用原生 ``formats.amber``
（决策 0015，不依赖 ParmEd）；项的顺序与方向与 prmtop 数组一致（先不含氢、后含氢）。
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem

from ffcompiler.formats.amber import PrmtopFile, read_prmtop, read_rst7

from ffcompiler.ir import (
    Angle,
    Atom,
    Bond,
    ForceFieldIR,
    ImproperTorsion,
    ProperTorsion,
    expected_exclusions_and_pairs,
    validate,
)

KCAL_KJ = 4.184
_K_BOND = 2.0 * KCAL_KJ * 100.0   # kcal/mol/Å² → kJ/mol/nm²（AMBER 无 ½，IR 有）
_K_ANGLE = 2.0 * KCAL_KJ          # kcal/mol/rad² → kJ/mol/rad²
_K_TORSION = KCAL_KJ              # kcal/mol → kJ/mol（prmtop 已除过 IDIVF）
_LJ_EPS = KCAL_KJ                 # kcal/mol → kJ/mol
_CHARGE_SCALE = 18.2223           # prmtop CHARGE = q(e) × 18.2223
_LB_REL_TOL = 1e-6                # 非对角 LJ 系数与 LB 组合的相对容差（16.8E 写出精度之上）


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _element(z: int) -> str:
    return Chem.GetPeriodicTable().GetElementSymbol(int(z))


def _read_sdf_formal_charges(sdf: Path, atomic_numbers: list[int]) -> list[int]:
    mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
    if mol is None:
        raise ValueError(f"cannot read SDF {sdf}")
    if mol.GetNumAtoms() != len(atomic_numbers):
        raise ValueError(
            f"SDF has {mol.GetNumAtoms()} atoms but prmtop has {len(atomic_numbers)}"
        )
    formal = []
    for atom, z in zip(mol.GetAtoms(), atomic_numbers):
        if atom.GetSymbol() != _element(z):
            raise ValueError(
                f"element mismatch at index {atom.GetIdx()}: SDF {atom.GetSymbol()!r} "
                f"vs prmtop {_element(z)!r} (atom order must correspond)"
            )
        formal.append(atom.GetFormalCharge())
    return formal


def _prmtop_excluded_pairs(prm: PrmtopFile) -> set[tuple[int, int]]:
    """prmtop EXCLUDED_ATOMS_LIST（1-2 ∪ 1-3 ∪ 1-4，文件里是 1-based）。

    只记前向引用（j > i）；没有排除项的原子以 0 哨兵占位，跳过。
    """
    counts = prm["NUMBER_EXCLUDED_ATOMS"]
    entries = prm["EXCLUDED_ATOMS_LIST"]
    pairs: set[tuple[int, int]] = set()
    pos = 0
    for i, count in enumerate(counts):
        for _ in range(int(count)):
            j = int(entries[pos]) - 1
            pos += 1
            if j < 0:
                continue  # 0 哨兵：该原子没有（前向）排除项
            pairs.add((min(i, j), max(i, j)))
    return pairs


def _lj_by_type(prm: PrmtopFile) -> tuple[list[float], list[float]]:
    """每个 LJ 类型的 (Rmin/2 Å, ε kcal/mol)，取 A/B 系数矩阵对角元；并核对非对角元
    符合 Lorentz–Berthelot（IR 只支持 LB）。"""
    ntypes = prm.pointer("NTYPES")
    ico = prm["NONBONDED_PARM_INDEX"]
    acoef, bcoef = prm["LENNARD_JONES_ACOEF"], prm["LENNARD_JONES_BCOEF"]
    half_rmin, eps = [], []
    for t in range(ntypes):
        idx = ico[ntypes * t + t] - 1
        a, b = (acoef[idx], bcoef[idx]) if idx >= 0 else (0.0, 0.0)
        if a < 1e-10 or b < 1e-10:
            half_rmin.append(0.0)
            eps.append(0.0)
        else:
            factor = 2.0 * a / b
            half_rmin.append(factor ** (1.0 / 6.0) * 0.5)
            eps.append(b / 2.0 / factor)
    for t1 in range(ntypes):
        for t2 in range(t1 + 1, ntypes):
            idx = ico[ntypes * t1 + t2] - 1
            if idx < 0:
                raise ValueError(f"LJ 类型对 ({t1 + 1},{t2 + 1}) 使用 10-12 势，IR 不支持")
            rmin = half_rmin[t1] + half_rmin[t2]
            e = math.sqrt(eps[t1] * eps[t2])
            for got, want in ((acoef[idx], e * rmin**12), (bcoef[idx], 2.0 * e * rmin**6)):
                if abs(got - want) > _LB_REL_TOL * max(abs(want), 1e-30):
                    raise ValueError(
                        f"LJ 类型对 ({t1 + 1},{t2 + 1}) 的系数不符合 Lorentz–Berthelot "
                        f"组合（{got:.8e} vs {want:.8e}）：IR 只支持 LB"
                    )
    return half_rmin, eps


def _triples(values: list[int], width: int) -> list[tuple[int, ...]]:
    return [tuple(values[k:k + width]) for k in range(0, len(values), width)]


def from_amber(
    prmtop: Path,
    inpcrd: Path,
    sdf: Path,
    name: str | None = None,
) -> tuple[ForceFieldIR, np.ndarray]:
    """按 §2.2 换算导入，返回 (ForceFieldIR, 坐标 nm)。结尾调用 validate。"""
    prmtop = Path(prmtop)
    inpcrd = Path(inpcrd)
    sdf = Path(sdf)
    prm = read_prmtop(prmtop)
    natom = prm.pointer("NATOM")
    if "ATOMIC_NUMBER" not in prm.flags:
        raise ValueError("prmtop 缺少 ATOMIC_NUMBER（不按质量猜元素）")
    atomic_numbers = [int(z) for z in prm["ATOMIC_NUMBER"]]
    formal_charges = _read_sdf_formal_charges(sdf, atomic_numbers)
    half_rmin, eps_type = _lj_by_type(prm)
    types = prm["AMBER_ATOM_TYPE"]
    res_ptr = list(prm["RESIDUE_POINTER"]) + [natom + 1]
    res_of = [""] * natom
    for r, label in enumerate(prm["RESIDUE_LABEL"]):
        for i in range(res_ptr[r] - 1, res_ptr[r + 1] - 1):
            res_of[i] = label

    atoms = []
    for i in range(natom):
        t = prm["ATOM_TYPE_INDEX"][i] - 1
        sigma_nm = half_rmin[t] * 2.0 * 2.0 ** (-1.0 / 6.0) * 0.1   # Rmin/2 Å → σ nm
        atoms.append(
            Atom(
                index=i,
                name=prm["ATOM_NAME"][i],
                element=_element(atomic_numbers[i]),
                atomic_number=atomic_numbers[i],
                mass=float(prm["MASS"][i]),
                atom_type=types[i],
                formal_charge=formal_charges[i],
                partial_charge=float(prm["CHARGE"][i]) / _CHARGE_SCALE,
                sigma=float(sigma_nm),
                epsilon=float(eps_type[t]) * _LJ_EPS,
                residue_name=res_of[i],
            )
        )

    bond_k, bond_req = prm["BOND_FORCE_CONSTANT"], prm["BOND_EQUIL_VALUE"]
    bond_entries = (_triples(prm["BONDS_WITHOUT_HYDROGEN"], 3)
                    + _triples(prm["BONDS_INC_HYDROGEN"], 3))
    bonds = tuple(
        Bond(
            i=i // 3, j=j // 3,
            r0=float(bond_req[t - 1]) * 0.1,
            k=float(bond_k[t - 1]) * _K_BOND,
            source="gaff2-library",
            meta={"amber_types": f"{types[i // 3]}-{types[j // 3]}"},
        )
        for i, j, t in bond_entries
    )

    ang_k, ang_eq = prm["ANGLE_FORCE_CONSTANT"], prm["ANGLE_EQUIL_VALUE"]
    angle_entries = (_triples(prm["ANGLES_WITHOUT_HYDROGEN"], 4)
                     + _triples(prm["ANGLES_INC_HYDROGEN"], 4))
    angles = tuple(
        Angle(
            i=i // 3, j=j // 3, k_atom=k // 3,
            theta0=float(ang_eq[t - 1]),
            force_constant=float(ang_k[t - 1]) * _K_ANGLE,
            source="gaff2-library",
            meta={"amber_types": f"{types[i // 3]}-{types[j // 3]}-{types[k // 3]}"},
        )
        for i, j, k, t in angle_entries
    )

    nptra = prm.pointer("NPTRA")
    scee_all = prm["SCEE_SCALE_FACTOR"] if "SCEE_SCALE_FACTOR" in prm.flags else [1.2] * nptra
    scnb_all = prm["SCNB_SCALE_FACTOR"] if "SCNB_SCALE_FACTOR" in prm.flags else [2.0] * nptra
    dih_entries = (_triples(prm["DIHEDRALS_WITHOUT_HYDROGEN"], 5)
                   + _triples(prm["DIHEDRALS_INC_HYDROGEN"], 5))
    # (i, j, k, l, type_index, ignore_end, improper)
    dihedrals = [
        (i // 3, j // 3, abs(k) // 3, abs(l) // 3, t - 1, k < 0, l < 0)
        for i, j, k, l, t in dih_entries
    ]

    propers: list[ProperTorsion] = []
    impropers: list[ImproperTorsion] = []
    scee_scnb_seen: set[tuple[float, float]] = set()
    for i, j, k, l, t, ignore_end, improper in dihedrals:
        types4 = "-".join(types[a] for a in (i, j, k, l))
        params = dict(
            periodicity=int(prm["DIHEDRAL_PERIODICITY"][t]),
            phase=float(prm["DIHEDRAL_PHASE"][t]),
            force_constant=float(prm["DIHEDRAL_FORCE_CONSTANT"][t]) * _K_TORSION,
            source="gaff2-library",
            meta={"amber_types": types4},
        )
        if improper:
            impropers.append(ImproperTorsion(i=i, j=j, k=k, l=l, **params))
            continue
        if not ignore_end:
            scee_scnb_seen.add((float(scee_all[t]), float(scnb_all[t])))
        propers.append(ProperTorsion(i=i, j=j, k=k, l=l, **params))

    if len(scee_scnb_seen) > 1:
        raise ValueError(
            f"non-uniform 1-4 scaling across prmtop dihedral entries: "
            f"SCEE/SCNB pairs seen = {sorted(scee_scnb_seen)}"
        )
    if scee_scnb_seen:
        scee, scnb = next(iter(scee_scnb_seen))
    else:
        scee, scnb = 1.2, 2.0  # AMBER 家族缺省
    if not all(math.isfinite(x) and x >= 1 for x in (scee, scnb)):
        raise ValueError(f"unsupported SCEE/SCNB: {scee}, {scnb}")
    coulomb_14_scale = 1.0 / scee
    vdw_14_scale = 1.0 / scnb

    ir = ForceFieldIR(
        name=name or prmtop.stem,
        atoms=tuple(atoms),
        bonds=bonds,
        angles=angles,
        proper_torsions=tuple(propers),
        impropers=tuple(impropers),
        exclusions=(),
        pairs14=(),
        total_charge=int(sum(formal_charges)),
        forcefield_family="amber",
        combination_rule="lorentz-berthelot",
        coulomb_14_scale=coulomb_14_scale,
        vdw_14_scale=vdw_14_scale,
        provenance={
            "inputs": {
                "prmtop": {"path": str(prmtop), "sha256": sha256_file(prmtop)},
                "inpcrd": {"path": str(inpcrd), "sha256": sha256_file(inpcrd)},
                "sdf": {"path": str(sdf), "sha256": sha256_file(sdf)},
            },
            "prmtop_reader": "ffcompiler.formats.amber",
        },
    )

    exclusions, pairs14 = expected_exclusions_and_pairs(ir)
    ir = ForceFieldIR(
        name=ir.name, atoms=ir.atoms, bonds=ir.bonds, angles=ir.angles,
        proper_torsions=ir.proper_torsions, impropers=ir.impropers,
        exclusions=exclusions, pairs14=pairs14,
        total_charge=ir.total_charge,
        forcefield_family=ir.forcefield_family,
        combination_rule=ir.combination_rule,
        coulomb_14_scale=ir.coulomb_14_scale, vdw_14_scale=ir.vdw_14_scale,
        provenance=ir.provenance,
    )

    # prmtop 自身排除表核对：1-2/1-3（EXCLUDED_ATOMS_LIST）∪ 1-4（计算 1-4 的 proper）
    prmtop_pairs = _prmtop_excluded_pairs(prm)
    from_dihedrals: Counter[tuple[int, int]] = Counter()
    for i, _, _, l, _, ignore_end, improper in dihedrals:
        if improper or ignore_end:
            continue
        from_dihedrals[(min(i, l), max(i, l))] += 1
    ir_pairs = set(exclusions) | set(pairs14)
    if prmtop_pairs != ir_pairs:
        only_prmtop = sorted(prmtop_pairs - ir_pairs)[:10]
        only_ir = sorted(ir_pairs - prmtop_pairs)[:10]
        raise ValueError(
            "prmtop exclusion tables disagree with topology-derived 1-2/1-3/1-4 sets: "
            f"only in prmtop: {only_prmtop}; only expected: {only_ir}"
        )

    # Exclusion membership alone does not prove an active 1-4 interaction.
    # AMBER computes each such interaction through one non-ignored proper.
    if set(from_dihedrals) != set(pairs14) or any(n != 1 for n in from_dihedrals.values()):
        raise ValueError(
            "unsupported active 1-4 carriers: "
            f"missing={sorted(set(pairs14) - set(from_dihedrals))[:10]}, "
            f"extra={sorted(set(from_dihedrals) - set(pairs14))[:10]}, "
            f"duplicates={[(p, n) for p, n in from_dihedrals.items() if n != 1][:10]}"
        )

    validate(ir)

    coords = read_rst7(inpcrd)
    if coords.shape[0] != natom:
        raise ValueError(f"inpcrd has {coords.shape[0]} atoms but prmtop has {natom}")
    positions_nm = coords * 0.1
    return ir, positions_nm
