"""GROMACS 导出器（MVP-0 计划 T8）：自包含的 topol.top + molecule.itp + conf.gro。

improper 写周期形式 funct 4、参数布局 `phase_deg k n`（见
docs/decisions/0003-gromacs-improper-funct4.md）。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ffcompiler.ir import ForceFieldIR, validate


def _g(num: float) -> str:
    """至少 10 位有效数字。"""
    return f"{num:.10g}"


def _pair(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _molecule_name(ir: ForceFieldIR) -> str:
    return ir.name.replace(" ", "_") or "MOL"


def _check_impropers_single_term(ir: ForceFieldIR) -> None:
    seen: dict[tuple[int, int, int, int], int] = {}
    for term in ir.impropers:
        key = min((term.i, term.j, term.k, term.l), (term.l, term.k, term.j, term.i))
        seen[key] = seen.get(key, 0) + 1
    for key, count in seen.items():
        if count > 1:
            raise ValueError(
                f"GROMACS funct 4 allows only one improper per quartet; "
                f"quartet {key} has {count} terms"
            )


def write_gromacs(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    outdir: Path,
) -> dict[str, Path]:
    validate(ir)
    _check_impropers_single_term(ir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    molname = _molecule_name(ir)
    resname = ir.atoms[0].residue_name

    # ---- molecule.itp ------------------------------------------------
    itp_lines: list[str] = []

    unique_types: dict[str, ForceFieldIR] = {}
    for atom in ir.atoms:
        if atom.atom_type not in unique_types:
            unique_types[atom.atom_type] = atom
    itp_lines.append("[ atomtypes ]")
    itp_lines.append("; name  at.num  mass  charge  ptype  sigma  epsilon")
    for atom_type, atom in sorted(unique_types.items()):
        itp_lines.append(
            f"{atom_type:<6} {atom.atomic_number:>6} "
            f"{_g(atom.mass):>14} 0.0000000000 A {_g(atom.sigma):>14} {_g(atom.epsilon):>14}"
        )

    itp_lines.append("")
    itp_lines.append("[ moleculetype ]")
    itp_lines.append(f"{molname} 3")

    itp_lines.append("[ atoms ]")
    itp_lines.append("; nr  type  resnr  residue  atom  cgnr  charge  mass")
    for n, atom in enumerate(ir.atoms, start=1):
        itp_lines.append(
            f"{n:>5} {atom.atom_type:<6} 1 {resname:<6} {atom.name:<6} "
            f"{n:>5} {atom.partial_charge:>18.10f} {_g(atom.mass):>14}"
        )

    itp_lines.append("[ bonds ]")
    for bond in ir.bonds:
        itp_lines.append(
            f"{bond.i + 1:>5} {bond.j + 1:>5} 1 {_g(bond.r0):>14} {_g(bond.k):>14}"
        )

    itp_lines.append("[ pairs ]")
    for i, j in ir.pairs14:
        itp_lines.append(f"{i + 1:>5} {j + 1:>5} 1")

    itp_lines.append("[ angles ]")
    for angle in ir.angles:
        itp_lines.append(
            f"{angle.i + 1:>5} {angle.j + 1:>5} {angle.k_atom + 1:>5} 1 "
            f"{_g(np.degrees(angle.theta0)):>14} {_g(angle.force_constant):>14}"
        )

    itp_lines.append("[ dihedrals ]")
    for term in ir.proper_torsions:
        itp_lines.append(
            f"{term.i + 1:>5} {term.j + 1:>5} {term.k + 1:>5} {term.l + 1:>5} 9 "
            f"{_g(np.degrees(term.phase)):>14} {_g(term.force_constant):>14} {term.periodicity}"
        )
    for term in ir.impropers:
        itp_lines.append(
            f"{term.i + 1:>5} {term.j + 1:>5} {term.k + 1:>5} {term.l + 1:>5} 4 "
            f"{_g(np.degrees(term.phase)):>14} {_g(term.force_constant):>14} {term.periodicity}"
        )

    itp_path = outdir / "molecule.itp"
    itp_path.write_text("\n".join(itp_lines) + "\n", encoding="utf-8")

    # ---- topol.top ----------------------------------------------------
    top_lines = [
        "[ defaults ]",
        "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
        f"1 2 yes {_g(ir.vdw_14_scale)} {_g(ir.coulomb_14_scale)}",
        "",
        f'#include "molecule.itp"',
        "",
        "[ system ]",
        ir.name,
        "",
        "[ molecules ]",
        f"{molname} 1",
        "",
    ]
    top_path = outdir / "topol.top"
    top_path.write_text("\n".join(top_lines), encoding="utf-8")

    # ---- conf.gro ------------------------------------------------------
    pos = np.asarray(positions_nm, dtype=np.float64)
    gro_lines = [f"{ir.name} written by ffcompiler", f"{len(ir.atoms):5d}"]
    for n, atom in enumerate(ir.atoms, start=1):
        x, y, z = pos[atom.index]
        gro_lines.append(
            f"{1:>5d}{resname:<5s}{atom.name[:5]:>5s}{n:5d}"
            f"{x:15.10f}{y:15.10f}{z:15.10f}"
        )
    if ir.box_nm is not None:
        b = np.asarray(ir.box_nm, dtype=np.float64)
        off = (b[0, 1], b[0, 2], b[1, 0], b[1, 2], b[2, 0], b[2, 1])
        if max(abs(v) for v in off) < 1e-12:
            gro_lines.append(f"{b[0, 0]:10.5f}{b[1, 1]:10.5f}{b[2, 2]:10.5f}")
        else:
            # 三斜：对角在前，随后 v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)（决策 0023）
            vals = (b[0, 0], b[1, 1], b[2, 2], *off)
            gro_lines.append("".join(f"{v:10.5f}" for v in vals))
    else:
        gro_lines.append("   10.00000   10.00000   10.00000")
    gro_path = outdir / "conf.gro"
    gro_path.write_text("\n".join(gro_lines) + "\n", encoding="utf-8")

    return {"itp": itp_path, "top": top_path, "gro": gro_path}
