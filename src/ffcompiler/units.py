"""单位换算常数（MVP-1 计划 §2.1，锁定值）。

内部单位约定（与 ForceFieldIR 一致）：长度 nm、能量 kJ/mol、角度 rad、
电荷 e、质量 amu。Hessian 子系统的输入输出全部使用这些单位。
"""

from __future__ import annotations

HARTREE_KJ_MOL = 2625.499639479
BOHR_NM = 0.0529177210544
ANGSTROM_NM = 0.1
KCAL_KJ = 4.184
C_CM_S = 2.99792458e10

# Hartree/Bohr² → kJ/mol/nm²（fchk、xtb 的 Hessian 单位换算）
HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2 = HARTREE_KJ_MOL / BOHR_NM**2

__all__ = [
    "HARTREE_KJ_MOL",
    "BOHR_NM",
    "ANGSTROM_NM",
    "KCAL_KJ",
    "C_CM_S",
    "HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2",
]
