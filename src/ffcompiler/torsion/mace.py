"""MACE-POLAR 参考势（MVP-2 计划 §0、§1.1，T27）。

直接用 mace 自带的 ASE 计算器（不经 openmm-ml：openmm-ml 1.6 不传 PolarMACE 需要的
fermi_level / external_field，决策 0016 修订）。torch / mace 只在本模块的函数里 import。
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

from ..ir.schema import ForceFieldIR

DEFAULT_MODEL = Path("/home/user/MLP/mace/MACE-POLAR-1-M.model")
MODEL_NAME = "MACE-POLAR-1-M"


def default_model_path() -> Path:
    path = Path(os.environ.get("FFC_MACE_MODEL", str(DEFAULT_MODEL)))
    if not path.is_file():
        raise FileNotFoundError(f"MACE 模型不存在：{path}（可用 FFC_MACE_MODEL 指定）")
    return path


def molecule_info(ir: ForceFieldIR) -> dict:
    """ASE atoms.info：总电荷原样传；只做闭壳层，spin（多重度）= 1。"""
    return {"charge": int(ir.total_charge), "spin": 1}


def polar_calculator(*, model_path: Path | None = None, device: str = "cpu", threads: int = 4):
    import torch
    from mace.calculators import MACECalculator

    torch.set_num_threads(int(threads))
    path = Path(model_path) if model_path is not None else default_model_path()
    return MACECalculator(model_paths=str(path), device=device, default_dtype="float64",
                          model_type="PolarMACE")


@lru_cache(maxsize=None)
def _atomic_numbers(path: str) -> frozenset[int]:
    import torch

    model = torch.load(path, map_location="cpu", weights_only=False)
    return frozenset(int(z) for z in model.atomic_numbers)


def supported_atomic_numbers(model_path: Path | None = None) -> frozenset[int]:
    path = Path(model_path) if model_path is not None else default_model_path()
    return _atomic_numbers(str(path))
