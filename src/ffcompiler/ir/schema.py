"""ForceFieldIR：引擎无关的完整分子 Hamiltonian（MVP-0 计划 T1）。

单位约定（计划 §1.2）：长度 nm、能量 kJ/mol、角度 rad、电荷 e、质量 amu；
文件中以度为单位的字段名以 ``_deg`` 结尾。IR 不含坐标。
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from typing import Any, Mapping

SCHEMA_VERSION = "ffir/0.1"


@dataclass(frozen=True)
class Atom:
    index: int
    name: str
    element: str
    atomic_number: int
    mass: float
    atom_type: str
    formal_charge: int
    partial_charge: float
    sigma: float
    epsilon: float
    residue_name: str = "MOL"


@dataclass(frozen=True)
class Bond:
    i: int
    j: int
    r0: float
    k: float
    source: str
    meta: Mapping[str, str | int | float | bool] = field(default_factory=dict)
    uncertainty: float | None = None


@dataclass(frozen=True)
class Angle:
    i: int
    j: int
    k_atom: int
    theta0: float
    force_constant: float
    source: str
    meta: Mapping[str, str | int | float | bool] = field(default_factory=dict)
    uncertainty: float | None = None


@dataclass(frozen=True)
class ProperTorsion:
    i: int
    j: int
    k: int
    l: int
    periodicity: int
    phase: float
    force_constant: float
    source: str
    meta: Mapping[str, str | int | float | bool] = field(default_factory=dict)
    uncertainty: float | None = None


@dataclass(frozen=True)
class ImproperTorsion:
    """与 ProperTorsion 字段完全相同；原子顺序就是求值顺序，不重排。"""

    i: int
    j: int
    k: int
    l: int
    periodicity: int
    phase: float
    force_constant: float
    source: str
    meta: Mapping[str, str | int | float | bool] = field(default_factory=dict)
    uncertainty: float | None = None


@dataclass(frozen=True)
class ForceFieldIR:
    name: str
    atoms: tuple[Atom, ...]
    bonds: tuple[Bond, ...]
    angles: tuple[Angle, ...]
    proper_torsions: tuple[ProperTorsion, ...]
    impropers: tuple[ImproperTorsion, ...]
    exclusions: tuple[tuple[int, int], ...]   # i<j，已排序
    pairs14: tuple[tuple[int, int], ...]      # i<j，已排序
    total_charge: int
    forcefield_family: str                    # MVP-0 只有 "amber"
    combination_rule: str                     # MVP-0 只有 "lorentz-berthelot"
    coulomb_14_scale: float
    vdw_14_scale: float
    # 周期盒（决策 0023）：3×3，行 = 晶胞向量 a/b/c [nm]；None = 孤立体系。
    # 导出器按引擎能力消费（GROMACS/OpenMM 全三斜；AMBER prmtop 仅正交）。
    box_nm: tuple[tuple[float, float, float], ...] | None = None
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def replace_terms(
        self,
        *,
        atoms: tuple[Atom, ...] | None = None,
        bonds: tuple[Bond, ...] | None = None,
        angles: tuple[Angle, ...] | None = None,
        proper_torsions: tuple[ProperTorsion, ...] | None = None,
        impropers: tuple[ImproperTorsion, ...] | None = None,
        provenance_update: Mapping[str, Any] | None = None,
    ) -> "ForceFieldIR":
        """返回新 IR；原对象不变。provenance_update 合并进新对象的 provenance。"""
        provenance = dict(self.provenance)
        if provenance_update:
            provenance.update(provenance_update)
        return ForceFieldIR(
            name=self.name,
            atoms=self.atoms if atoms is None else tuple(atoms),
            bonds=self.bonds if bonds is None else tuple(bonds),
            angles=self.angles if angles is None else tuple(angles),
            proper_torsions=(
                self.proper_torsions if proper_torsions is None else tuple(proper_torsions)
            ),
            impropers=self.impropers if impropers is None else tuple(impropers),
            exclusions=self.exclusions,
            pairs14=self.pairs14,
            total_charge=self.total_charge,
            forcefield_family=self.forcefield_family,
            combination_rule=self.combination_rule,
            coulomb_14_scale=self.coulomb_14_scale,
            vdw_14_scale=self.vdw_14_scale,
            provenance=provenance,
        )


__all__ = [
    "SCHEMA_VERSION",
    "Atom",
    "Bond",
    "Angle",
    "ProperTorsion",
    "ImproperTorsion",
    "ForceFieldIR",
]
