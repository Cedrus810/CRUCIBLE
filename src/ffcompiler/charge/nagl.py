"""NAGL 电荷（openff-gnn-am1bcc，MVP-3 计划 §2.3）。唯一 import openff 的模块，且只在函数内 import。"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ffcompiler.charge.base import ChargeError, ChargeProvider
from ffcompiler.torsion.rotors import check_mol_topology

DEFAULT_MODEL = "openff-gnn-am1bcc-1.0.0.pt"


def nagl_model_path(model: str) -> Path:
    import openff.nagl_models

    path = Path(openff.nagl_models.__file__).parent / "models" / "am1bcc" / model
    if not path.is_file():
        raise ChargeError(f"NAGL model not found: {path}")
    return path


def nagl_raw_charges(mol, model: str) -> np.ndarray:
    """RDKit 分子（显式氢）→ NAGL 原始电荷；适用域外或模型出错 → ChargeError。"""
    from openff.nagl import GNNModel
    from openff.toolkit import Molecule

    try:
        off = Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
        gnn = GNNModel.load(str(nagl_model_path(model)), eval_mode=True)
        q = gnn.compute_property(off, check_domains=True)
    except ChargeError:
        raise
    except Exception as exc:  # noqa: BLE001 - openff 的各种报错统一转成 ChargeError，保留原文
        raise ChargeError(f"NAGL {model} failed: {exc!r}") from exc
    return np.asarray(q, dtype=np.float64)


class NAGLCharges(ChargeProvider):
    method = "nagl"

    def __init__(self, model: str = DEFAULT_MODEL):
        self.model = model

    def compute(self, ir, positions_nm, mol) -> np.ndarray:
        check_mol_topology(ir, mol)  # 原子序数与键集合（只比元素会放过 H 对调，审阅 Important 1）
        return nagl_raw_charges(mol, self.model)

    def source_info(self) -> dict:
        import openff.nagl

        path = nagl_model_path(self.model)
        st = path.stat()
        return {"model": self.model, "path": str(path), "size": st.st_size, "mtime": st.st_mtime,
                "openff_nagl": openff.nagl.__version__}
