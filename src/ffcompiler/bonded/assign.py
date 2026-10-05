"""驱动函数 derive_bonded 与报告（MVP-1 计划 T21）。

原子顺序、拓扑一致性检查（Review Focus 1）、线性键角策略（Review Focus 2）
和写回 IR 都在这里；参数推导只调 T17/T18/T19 的方法函数。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from rdkit import Chem

from ffcompiler.bonded.drih import drih
from ffcompiler.bonded.internal import (
    LinearAngleError,
    angle_value,
    bond_value,
    linear_angles,
    proper_quartets,
)
from ffcompiler.bonded.mseminario import mseminario
from ffcompiler.bonded.symmetry import symmetrize as _symmetrize_terms
from ffcompiler.bonded.symmetry import term_classes
from ffcompiler.bonded.terms import DerivedTerm
from ffcompiler.bonded.m2seminario import m2seminario
from ffcompiler.hessian.base import (
    CartesianHessian,
    HessianDiagnostics,
    check_hessian,
)
from ffcompiler.ir import Angle, Bond, ForceFieldIR
from ffcompiler.ir.topology import bond_graph
from ffcompiler.ir.validation import validate

BondedMethod = Literal["mseminario", "m2seminario", "drih"]
# 命名按计划 §2.6：C2 与 Sobtop m2Seminario 全部一致（决策 0011）

_BOND_RATIO_MIN = 0.8
_BOND_RATIO_MAX = 1.25
_COVALENT_FACTOR = 1.1


class TopologyMismatchError(ValueError):
    """Hessian 的原子顺序或几何与 IR 拓扑不符（计划 Review Focus 1）。"""


class NonPositiveForceConstantError(ValueError):
    """方法在这个 Hessian 上给出非正 k（方法局限，决策 0014）。不自动回退到库参数。"""

    def __init__(self, method: str, terms: list[tuple[tuple[int, ...], float]]):
        self.method = method
        self.terms = terms
        listed = ", ".join(f"{atoms}: k={k:.4g}" for atoms, k in terms)
        super().__init__(
            f"{method} 在该 Hessian 上给出 {len(terms)} 个非正力常数（{listed}）。"
            "这是方法局限而非输入错误；可改用 mseminario 或 drih。"
        )


@dataclass(frozen=True)
class BondedReport:
    method: str
    hessian_provider: str
    hessian_method: str
    hessian_sha256: str | None
    diagnostics: HessianDiagnostics
    linear_angles_kept: tuple[tuple[int, int, int], ...]
    class_spread: dict[str, float]
    max_bond_ratio: float


def _check_topology(
    ir: ForceFieldIR,
    h: CartesianHessian,
    input_positions_nm: np.ndarray | None,
) -> float:
    pos = h.positions_nm
    max_dev = 0.0
    if input_positions_nm is not None:
        for bond in ir.bonds:
            r_input = bond_value(input_positions_nm, bond.i, bond.j)
            r_hess = bond_value(pos, bond.i, bond.j)
            ratio = r_hess / r_input
            max_dev = max(max_dev, abs(ratio - 1.0))
            if not _BOND_RATIO_MIN <= ratio <= _BOND_RATIO_MAX:
                raise TopologyMismatchError(
                    f"键 ({bond.i}, {bond.j}) 长度比超界：Hessian 几何 "
                    f"{r_hess:.4f} nm / 输入 {r_input:.4f} nm = {ratio:.3f}，"
                    f"允许 [{_BOND_RATIO_MIN}, {_BOND_RATIO_MAX}]"
                )
    table = Chem.GetPeriodicTable()
    radii = [table.GetRcovalent(int(z)) * 0.1 for z in h.atomic_numbers]
    bonded = bond_graph(ir)
    n = len(ir.atoms)
    for a in range(n):
        for b in range(a + 1, n):
            if b in bonded.get(a, frozenset()):
                continue
            distance = float(np.linalg.norm(pos[a] - pos[b]))
            threshold = _COVALENT_FACTOR * (radii[a] + radii[b])
            if distance < threshold:
                raise TopologyMismatchError(
                    f"非键原子对 ({a}, {b}) 距离 {distance:.4f} nm < "
                    f"{threshold:.4f} nm（1.1×共价半径和）：Hessian 几何与拓扑不符"
                )
    return max_dev


def check_hessian_topology(
    ir: ForceFieldIR, h: CartesianHessian, input_positions_nm: np.ndarray | None
) -> float:
    """Review Focus 1：原子顺序与拓扑一致性（T21 规格 1–2）。返回键长比最大偏离。

    derive_bonded 与 global_fit 共用。
    """
    if list(h.atomic_numbers) != [a.atomic_number for a in ir.atoms]:
        raise TopologyMismatchError(
            "Hessian 原子顺序与 IR 不同（不做自动重排）："
            f"Hessian={list(h.atomic_numbers)} IR={[a.atomic_number for a in ir.atoms]}"
        )
    return _check_topology(ir, h, input_positions_nm)


def _derive_terms(
    h: CartesianHessian,
    ir: ForceFieldIR,
    method: BondedMethod,
    bonds: list[tuple[int, int]],
    angles: list[tuple[int, int, int]],
    allow_nonstationary: bool,
) -> list[DerivedTerm]:
    if method == "mseminario":
        return mseminario(h, bonds, angles)
    if method == "m2seminario":
        return m2seminario(h, bonds, angles)
    if method == "drih":
        return drih(
            h, bonds, angles, dihedrals=proper_quartets(bond_graph(ir)),
            allow_nonstationary=allow_nonstationary,
        )
    raise ValueError(f"未知方法 {method!r}")


def derive_bonded(
    ir: ForceFieldIR,
    h: CartesianHessian,
    *,
    method: BondedMethod,
    symmetrize: bool = True,
    linear_angle_policy: Literal["error", "library"] = "error",
    allow_nonstationary: bool = False,
    input_positions_nm: np.ndarray | None = None,
) -> tuple[ForceFieldIR, np.ndarray, BondedReport]:
    """从 CartesianHessian 推导 bond/angle 并写回 IR。

    返回 (新 IR, Hessian 几何 nm, BondedReport)。torsion、improper、非键、
    exclusions、pairs14 不变；结尾调用 validate。
    """
    max_bond_ratio = check_hessian_topology(ir, h, input_positions_nm)
    diagnostics = check_hessian(
        h, np.array([a.mass for a in ir.atoms]), allow_nonstationary=allow_nonstationary
    )

    bonds = [(b.i, b.j) for b in ir.bonds]
    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]

    def _canonical_angle(triple: tuple[int, int, int]) -> tuple[int, int, int]:
        i, j, k = triple
        return (min(i, k), j, max(i, k))

    linear = linear_angles(h.positions_nm, angles)
    if linear and linear_angle_policy == "error":
        raise LinearAngleError(linear)
    kept_linear = {_canonical_angle(a) for a in linear}
    angles_for_method = [a for a in angles if _canonical_angle(a) not in kept_linear]

    terms = _derive_terms(h, ir, method, bonds, angles_for_method, allow_nonstationary)
    nonpositive = [(t.atoms, t.k) for t in terms if not t.k > 0.0]
    if nonpositive:
        raise NonPositiveForceConstantError(method, nonpositive)
    classes = term_classes(ir, [t.atoms for t in terms])
    averaged_terms, spread_by_class = _symmetrize_terms(terms, classes)
    labels: dict[int, str] = {}
    for term, cls in zip(terms, classes):
        labels.setdefault(cls, _term_label(term.atoms))
    class_spread = {labels[cls]: value for cls, value in spread_by_class.items()}
    # 每项所在等价类的离散度（规格 3：meta 含 class_spread）
    spread_of_term = {
        _term_key(term.atoms): spread_by_class[cls] for term, cls in zip(terms, classes)
    }
    if symmetrize:
        terms = averaged_terms

    meta_common = {
        "hessian_provider": h.provider,
        "hessian_method": h.method,
        "hessian_sha256": h.source_sha256 or "",
        "symmetrized": symmetrize,
    }
    new_bonds: list[Bond] = []
    seen_bond: dict[tuple[int, int], Bond] = {}
    for term in terms:
        if len(term.atoms) != 2:
            continue
        i, j = term.atoms
        bond = Bond(
            i=i, j=j, r0=term.value, k=term.k, source=method,
            meta=dict(meta_common, class_spread=spread_of_term[_term_key(term.atoms)]),
        )
        seen_bond[(min(i, j), max(i, j))] = bond
    for bond in ir.bonds:
        new_bonds.append(seen_bond[(min(bond.i, bond.j), max(bond.i, bond.j))])

    derived_angles: dict[tuple[int, int, int], DerivedTerm] = {}
    for term in terms:
        if len(term.atoms) == 3:
            i, j, k = term.atoms
            derived_angles[_canonical_angle((i, j, k))] = term
    new_angles: list[Angle] = []
    for angle in ir.angles:
        triple = _canonical_angle((angle.i, angle.j, angle.k_atom))
        if triple in kept_linear:
            new_angles.append(
                Angle(
                    i=angle.i, j=angle.j, k_atom=angle.k_atom,
                    theta0=angle.theta0, force_constant=angle.force_constant,
                    source=angle.source,
                    meta=dict(angle.meta, linear_angle=True),
                )
            )
            continue
        term = derived_angles[triple]
        new_angles.append(
            Angle(
                i=triple[0], j=triple[1], k_atom=triple[2],
                theta0=term.value, force_constant=term.k, source=method,
                meta=dict(meta_common, class_spread=spread_of_term[triple]),
            )
        )

    new_ir = ir.replace_terms(
        bonds=tuple(new_bonds),
        angles=tuple(new_angles),
        provenance_update={
            "bonded": {
                "method": method,
                "hessian_provider": h.provider,
                "hessian_method": h.method,
                "hessian_sha256": h.source_sha256,
                "symmetrized": symmetrize,
                "linear_angle_policy": linear_angle_policy,
                "linear_angles_kept": [list(a) for a in sorted(kept_linear)],
                "max_class_spread": max(class_spread.values(), default=0.0),
                "max_bond_ratio_dev": max_bond_ratio,
            }
        },
    )
    validate(new_ir)
    report = BondedReport(
        method=method,
        hessian_provider=h.provider,
        hessian_method=h.method,
        hessian_sha256=h.source_sha256,
        diagnostics=diagnostics,
        linear_angles_kept=tuple(sorted(kept_linear)),
        class_spread=class_spread,
        max_bond_ratio=max_bond_ratio,
    )
    return new_ir, h.positions_nm, report


def _term_key(atoms: tuple[int, ...]) -> tuple[int, ...]:
    """bond 取无序对、angle 取端点无序的规范形式。"""
    if len(atoms) == 2:
        return (min(atoms), max(atoms))
    i, j, k = atoms
    return (min(i, k), j, max(i, k))


def _term_label(atoms: tuple[int, ...]) -> str:
    if len(atoms) == 2:
        return f"bond({atoms[0]},{atoms[1]})"
    return f"angle({atoms[0]},{atoms[1]},{atoms[2]})"
