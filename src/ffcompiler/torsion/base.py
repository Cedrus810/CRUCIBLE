"""TorsionProvider 接口（设计 §10）。"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from rdkit import Chem

from ..ir.schema import ForceFieldIR


class TorsionProvider(ABC):
    @abstractmethod
    def parameterize(self, ir: ForceFieldIR, positions_nm: np.ndarray, mol: Chem.Mol):
        """返回 (新 IR, 报告)。新 IR 只改 proper torsion，并更新 provenance。"""
