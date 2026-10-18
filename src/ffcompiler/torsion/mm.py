"""MM 工具（MVP-2 计划 §2.3，T26）：IR → ASE 计算器；去掉某根键二面角后的 MM 能量。

全部走 OpenMM Reference（双精度）。ASE 单位 eV、Å；对外 kJ/mol、nm。
"""

from __future__ import annotations

import numpy as np
import openmm
from ase.calculators.calculator import Calculator, all_changes
from openmm import unit

from ..exporters.openmm import to_openmm_system
from ..ir.schema import ForceFieldIR
from ..units import ANGSTROM_NM, EV_KJ_MOL
from ..validation.engines import reference_context

_KJ = unit.kilojoule_per_mole
_KJ_NM = unit.kilojoule_per_mole / unit.nanometer


class OpenMMCalculator(Calculator):
    """把 IR 的 MM 势包成 ASE 计算器（测试参考势、自洽检查用）。"""

    implemented_properties = ["energy", "forces"]

    def __init__(self, ir: ForceFieldIR, *, validate_ir: bool = True, **kwargs) -> None:
        super().__init__(**kwargs)
        self._context = reference_context(to_openmm_system(ir, validate_ir=validate_ir))

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self._context.setPositions(self.atoms.positions * ANGSTROM_NM)
        state = self._context.getState(getEnergy=True, getForces=True)
        e = state.getPotentialEnergy().value_in_unit(_KJ)
        f = state.getForces(asNumpy=True).value_in_unit(_KJ_NM)
        self.results = {
            "energy": e / EV_KJ_MOL,
            "forces": np.asarray(f) * ANGSTROM_NM / EV_KJ_MOL,  # kJ/mol/nm → eV/Å
        }


def strip_bond_torsions(ir: ForceFieldIR, bond: tuple[int, int]) -> ForceFieldIR:
    """去掉所有中心键为 {j,k} 的 proper torsion，其余不变。"""
    j, k = bond
    kept = tuple(t for t in ir.proper_torsions if {t.j, t.k} != {j, k})
    return ir.replace_terms(proper_torsions=kept)


def mm_energy(ir: ForceFieldIR, positions_nm):
    """完整 MM 能量（kJ/mol）；也用于不完整的扣除势，不在这里校验 IR。"""
    ctx = reference_context(to_openmm_system(ir, validate_ir=False))
    pos = np.asarray(positions_nm, dtype=np.float64)
    stack = pos[None] if pos.ndim == 2 else pos
    out = np.empty(len(stack))
    for m, x in enumerate(stack):
        ctx.setPositions(x)
        out[m] = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(_KJ)
    return float(out[0]) if pos.ndim == 2 else out


def mm_energy_without_bond(ir: ForceFieldIR, bond: tuple[int, int], positions_nm):
    """E_MM[IR \\ b]（kJ/mol）。positions (N,3) 返回 float，(M,N,3) 返回 (M,) 数组。"""
    return mm_energy(strip_bond_torsions(ir, bond), positions_nm)
