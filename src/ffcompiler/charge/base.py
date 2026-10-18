"""电荷 provider 接口与残差均摊（MVP-3 计划 §2.1、§2.3）。"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass

import numpy as np


class ChargeError(RuntimeError):
    """电荷计算失败：外部程序出错、模型适用域外等。"""


class ChargeProvider(ABC):
    """provider 只算原始电荷（单位 e，未均摊）；写回、均摊、provenance 由 apply.assign_charges 负责。"""

    method: str = ""

    @abstractmethod
    def compute(self, ir, positions_nm: np.ndarray, mol) -> np.ndarray: ...

    @abstractmethod
    def source_info(self) -> dict: ...


def equalize_residual(q: np.ndarray, total_charge: int) -> tuple[np.ndarray, float]:
    """把 Q − Σq 均摊到每个原子，末原子吸收剩余（与 scripts/make_fixtures.py 的规则相同）。"""
    q = np.asarray(q, dtype=np.float64).copy()
    residual = float(total_charge - q.sum())
    if abs(residual) <= 1e-12:
        return q, residual
    q += residual / len(q)
    q[-1] += total_charge - q.sum()
    return q, residual


@dataclass
class ChargeReport:
    method: str
    raw_sum: float
    residual: float
    source: dict
    torsion_refit: dict | None = None

    def to_dict(self) -> dict:
        return asdict(self)
