"""IR 校验（MVP-0 计划 T2）。发现任何问题就一次性抛出全部 issues。"""

from __future__ import annotations

import math
from typing import Iterable

from .schema import ForceFieldIR
from .topology import expected_exclusions_and_pairs

PERIODICITIES = set(range(1, 7))


class IRValidationError(ValueError):
    def __init__(self, issues: Iterable[str]) -> None:
        self.issues: list[str] = list(issues)
        super().__init__(" | ".join(self.issues) or "unknown IR validation error")


def _canon_quartet(i: int, j: int, k: int, l: int) -> tuple[int, int, int, int]:
    """(i,j,k,l) 与 (l,k,j,i) 视为同一个四元组，取字典序小者。"""
    forward = (i, j, k, l)
    backward = (l, k, j, i)
    return min(forward, backward)


def _check_indices(ir: ForceFieldIR, issues: list[str]) -> None:
    n = len(ir.atoms)
    for pos, atom in enumerate(ir.atoms):
        if atom.index != pos:
            issues.append(
                f"atom index broken: atoms[{pos}].index == {atom.index}, expected {pos} "
                "(indices must be exactly 0..N-1 in order)"
            )
            break
    valid = set(range(n))
    for bond in ir.bonds:
        for idx in (bond.i, bond.j):
            if idx not in valid:
                issues.append(f"bond references out-of-range atom index {idx} (N={n})")
    for angle in ir.angles:
        for idx in (angle.i, angle.j, angle.k_atom):
            if idx not in valid:
                issues.append(f"angle references out-of-range atom index {idx} (N={n})")
    for label, terms in (("proper", ir.proper_torsions), ("improper", ir.impropers)):
        for term in terms:
            for idx in (term.i, term.j, term.k, term.l):
                if idx not in valid:
                    issues.append(f"{label} torsion references out-of-range atom index {idx} (N={n})")


def _check_connectivity(ir: ForceFieldIR, issues: list[str]) -> None:
    bonds = {frozenset((b.i, b.j)) for b in ir.bonds}

    def has_bond(a: int, b: int) -> bool:
        return frozenset((a, b)) in bonds

    for angle in ir.angles:
        for a, b in ((angle.i, angle.j), (angle.j, angle.k_atom)):
            if not has_bond(a, b):
                issues.append(
                    f"angle ({angle.i},{angle.j},{angle.k_atom}) lacks bond ({a},{b})"
                )
    for term in ir.proper_torsions:
        for a, b in ((term.i, term.j), (term.j, term.k), (term.k, term.l)):
            if not has_bond(a, b):
                issues.append(
                    f"proper torsion ({term.i},{term.j},{term.k},{term.l}) lacks bond ({a},{b})"
                )
    for term in ir.impropers:
        # 中心原子可以在四元组的任何位置（决策 0002：tleap 写出的 improper
        # 中心原子不总在第三位，能量按存储顺序求值）。要求：恰好有一个原子
        # 与其余三个都成键。
        members = (term.i, term.j, term.k, term.l)
        centers = [
            m for n, m in enumerate(members)
            if all(has_bond(m, o) for o in members if o != m)
        ]
        if len(centers) != 1:
            issues.append(
                f"improper ({term.i},{term.j},{term.k},{term.l}) does not have exactly "
                f"one atom bonded to the other three (candidates: {centers})"
            )


def _check_duplicates(ir: ForceFieldIR, issues: list[str]) -> None:
    seen_bonds: set[frozenset[int]] = set()
    for bond in ir.bonds:
        key = frozenset((bond.i, bond.j))
        if key in seen_bonds:
            issues.append(f"duplicate bond ({bond.i},{bond.j})")
        seen_bonds.add(key)

    seen_angles: set[tuple[int, int, int]] = set()
    for angle in ir.angles:
        key = (min(angle.i, angle.k_atom), angle.j, max(angle.i, angle.k_atom))
        if key in seen_angles:
            issues.append(
                f"duplicate angle ({angle.i},{angle.j},{angle.k_atom})"
            )
        seen_angles.add(key)

    for label, terms in (("proper", ir.proper_torsions), ("improper", ir.impropers)):
        seen_quartets: dict[tuple[int, int, int, int], set[int]] = {}
        for term in terms:
            key = _canon_quartet(term.i, term.j, term.k, term.l)
            periodicities = seen_quartets.setdefault(key, set())
            if term.periodicity in periodicities:
                issues.append(
                    f"duplicate {label} torsion periodicity {term.periodicity} "
                    f"on quartet ({term.i},{term.j},{term.k},{term.l})"
                )
            periodicities.add(term.periodicity)


def _check_parameter_ranges(ir: ForceFieldIR, issues: list[str]) -> None:
    for atom in ir.atoms:
        for name, positive in (("mass", True), ("sigma", False), ("epsilon", False)):
            value = getattr(atom, name)
            if not math.isfinite(value) or (value <= 0 if positive else value < 0):
                issues.append(f"atom {atom.index} has invalid {name}={value}; must be finite and "
                              + ("> 0" if positive else ">= 0"))
        if not math.isfinite(atom.partial_charge):
            issues.append(f"atom {atom.index} has non-finite partial_charge={atom.partial_charge}")
    for bond in ir.bonds:
        if not math.isfinite(bond.k) or not bond.k > 0:
            issues.append(f"bond ({bond.i},{bond.j}) has k={bond.k}, must be > 0")
        if not math.isfinite(bond.r0) or not bond.r0 > 0:
            issues.append(f"bond ({bond.i},{bond.j}) has r0={bond.r0}, must be > 0")
    for angle in ir.angles:
        if not math.isfinite(angle.force_constant) or not angle.force_constant > 0:
            issues.append(
                f"angle ({angle.i},{angle.j},{angle.k_atom}) has "
                f"force_constant={angle.force_constant}, must be > 0"
            )
        if not 0 < angle.theta0 <= math.pi:
            issues.append(
                f"angle ({angle.i},{angle.j},{angle.k_atom}) has theta0={angle.theta0}, "
                f"must be in (0, pi]"
            )
    for label, terms in (("proper", ir.proper_torsions), ("improper", ir.impropers)):
        for term in terms:
            for name in ("phase", "force_constant"):
                if not math.isfinite(getattr(term, name)):
                    issues.append(f"{label} torsion ({term.i},{term.j},{term.k},{term.l}) "
                                  f"has non-finite {name}")
            if term.periodicity not in PERIODICITIES:
                issues.append(
                    f"{label} torsion ({term.i},{term.j},{term.k},{term.l}) has "
                    f"periodicity={term.periodicity}, must be in 1..6"
                )
    for label, terms in (("bond", ir.bonds), ("angle", ir.angles),
                         ("proper", ir.proper_torsions), ("improper", ir.impropers)):
        for term in terms:
            if term.uncertainty is not None and (
                not math.isfinite(term.uncertainty) or term.uncertainty < 0
            ):
                issues.append(f"{label} has invalid uncertainty={term.uncertainty}")


def _check_nonbonded_bookkeeping(ir: ForceFieldIR, issues: list[str]) -> None:
    exp_exc, exp_pairs = expected_exclusions_and_pairs(ir)
    if set(ir.exclusions) != set(exp_exc):
        issues.append(
            f"exclusions != expected exclusions (1-2 union 1-3): "
            f"got {len(ir.exclusions)} pairs, expected {len(exp_exc)}"
        )
    if set(ir.pairs14) != set(exp_pairs):
        issues.append(
            f"pairs14 != expected pairs14 (shortest path exactly 3 bonds): "
            f"got {len(ir.pairs14)} pairs, expected {len(exp_pairs)}"
        )
    overlap = set(ir.exclusions) & set(ir.pairs14)
    if overlap:
        issues.append(f"pairs appear in both exclusions and pairs14: {sorted(overlap)[:5]}")
    for name, pairs in (("exclusions", ir.exclusions), ("pairs14", ir.pairs14)):
        bad = [p for p in pairs if p[0] >= p[1]]
        if bad:
            issues.append(f"{name} must be stored with i<j; offenders: {bad[:5]}")


def _check_charge(ir: ForceFieldIR, issues: list[str]) -> None:
    total = sum(a.partial_charge for a in ir.atoms)
    if not math.isfinite(total) or not math.isfinite(ir.total_charge) or abs(total - ir.total_charge) > 1e-4:
        issues.append(
            f"charge not conserved: sum(partial_charge) = {total!r}, "
            f"total_charge = {ir.total_charge}"
        )


def _check_atom_type_consistency(ir: ForceFieldIR, issues: list[str]) -> None:
    by_type: dict[str, list] = {}
    for atom in ir.atoms:
        by_type.setdefault(atom.atom_type, []).append(atom)
    for atom_type, group in sorted(by_type.items()):
        for field_name in ("sigma", "epsilon", "mass", "atomic_number"):
            values = {getattr(a, field_name) for a in group}
            if len(values) > 1:
                issues.append(
                    f"atom_type {atom_type!r} has inconsistent {field_name}: {sorted(values)}"
                )


def _check_family(ir: ForceFieldIR, issues: list[str]) -> None:
    if ir.forcefield_family != "amber":
        issues.append(
            f"forcefield_family = {ir.forcefield_family!r}, MVP-0 only supports 'amber'"
        )
    if ir.combination_rule != "lorentz-berthelot":
        issues.append(
            f"combination_rule = {ir.combination_rule!r}, only 'lorentz-berthelot' is supported"
        )
    for name in ("coulomb_14_scale", "vdw_14_scale"):
        value = getattr(ir, name)
        if not 0 < value <= 1:
            issues.append(f"{name} = {value!r}, must be in (0, 1]")


def validate(ir: ForceFieldIR) -> None:
    issues: list[str] = []
    _check_indices(ir, issues)
    _check_connectivity(ir, issues)
    _check_duplicates(ir, issues)
    _check_parameter_ranges(ir, issues)
    _check_nonbonded_bookkeeping(ir, issues)
    _check_charge(ir, issues)
    _check_atom_type_consistency(ir, issues)
    _check_family(ir, issues)
    if issues:
        raise IRValidationError(issues)
