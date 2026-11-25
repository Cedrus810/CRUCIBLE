"""xtb GFN2 作二面角扫描参考势（决策 0028：离子/两性离子默认用 GFN2，不用 MACE）。

openmm_dev 里没有 tblite / xtb-python（生产环境不装包），所以计算器每次求能量与力都调用
外部 xtb 单点 + 梯度（``hessian.xtb.run_xtb_singlepoint``，``--gfn 2 --grad --chrg Q --uhf 0``）。
xtb 的工作文件放在系统临时目录（``TMPDIR`` 可覆盖），计算器回收时删除——它们是逐步的
瞬时文件，不是产物；扫描结果照常进 scan_cache。

身份（进缓存键与 provenance，refit 只读缓存时据此复原）：模型名 ``xtb-gfn2``；"模型文件"
记 xtb 可执行文件的路径/大小/修改时间（同 MACE 模型文件的记法，不用 hash）。
"""

from __future__ import annotations

import shutil
import tempfile
import weakref
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from ..ir.schema import ForceFieldIR
from ..units import BOHR_NM, HARTREE_KJ_MOL

MODEL_NAME = "xtb-gfn2"
SUPPORTED = frozenset(range(1, 87))   # GFN2-xTB 参数化覆盖 H–Rn

_EV_PER_HARTREE = 27.211386245988     # CODATA 2018
_ANGSTROM_PER_BOHR = BOHR_NM * 10.0

__all__ = ["MODEL_NAME", "SUPPORTED", "XtbGFN2Calculator", "is_ion_or_zwitterion",
           "xtb_binary"]


def is_ion_or_zwitterion(ir: ForceFieldIR) -> bool:
    """决策 0028 的判定：净电荷 ≠ 0，或任一原子带形式电荷（两性离子净电荷为 0 但有 ±）。"""
    return int(ir.total_charge) != 0 or any(int(a.formal_charge) != 0 for a in ir.atoms)


def xtb_binary() -> Path:
    """按 toolconfig 解析 xtb（FFC_XTB → 配置 → 旧默认，决策 0026）。"""
    from ..toolconfig import resolve

    return Path(resolve("xtb").path)


class XtbGFN2Calculator(Calculator):
    """ASE 计算器：energy（eV）与 forces（eV/Å），GFN2-xTB 单点 + 梯度。

    ``charge`` 是体系总电荷（闭壳层，uhf 0；开壳层由 run_xtb_singlepoint 运行前拒绝）。
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, charge: int, *, threads: int = 4, binary: Path | None = None):
        super().__init__()
        self.charge = int(charge)
        self.threads = int(threads)
        self.binary = Path(binary) if binary is not None else xtb_binary()
        self.workdir = Path(tempfile.mkdtemp(prefix="ffc_xtb_gfn2_"))
        self._cleanup = weakref.finalize(self, shutil.rmtree, self.workdir, True)
        self.n_calls = 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        from ..hessian.xtb import run_xtb_singlepoint

        out = run_xtb_singlepoint(
            np.asarray(self.atoms.numbers), np.asarray(self.atoms.positions) * 0.1,  # Å → nm
            method="gfn2", charge=self.charge, workdir=self.workdir,
            threads=self.threads, binary=self.binary)
        self.n_calls += 1
        self.results["energy"] = float(out["energy_hartree"]) * _EV_PER_HARTREE
        self.results["forces"] = (-np.asarray(out["gradient_hartree_per_bohr"])
                                  * (_EV_PER_HARTREE / _ANGSTROM_PER_BOHR))


# 单位一致性自检（导入时，纯算术）：Eh 换 kJ/mol 与 eV 的比值应为 96.485 kJ/mol/eV
assert abs(HARTREE_KJ_MOL / _EV_PER_HARTREE - 96.4853321) < 1e-6
