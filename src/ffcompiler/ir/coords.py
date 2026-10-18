"""坐标 sidecar：与 IR JSON 并行传递的坐标合同（2026-10-15 计划 C1，2026-10-17 施工）。

IR schema 按 MVP-0 设计**不含坐标**；跨命令传递"上一步的几何"（特别是
Hessian 优化几何）用本模块的独立 JSON：

- 显式配对：读取方必须传文件路径（CLI ``--coords-in``），不按同名文件猜测；
- 显式核对：坐标的原子序数/元素/原子名序列必须与 IR 逐位一致，否则报错拒绝
  （乱序或错拓扑不允许自动重排；带原子名才能查出同元素原子换序，
  2026-10-17 验收后由 ffcoords/0.1 升到 0.2）；
- 来源可追溯：``source`` 记录坐标从哪来（如 ``inpcrd:...``、``hessian:gfn2``），
  ``units`` 固定 ``nm``（schema 版本锁定）。

浮点数与 ir/io.py 同约定：``json`` 的 ``repr`` 最短往返表示，不做格式化截断。
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .schema import ForceFieldIR

COORDS_SCHEMA_VERSION = "ffcoords/0.2"


class CoordinatesMismatchError(ValueError):
    """坐标与 IR 的原子序数/元素/原子名序列不一致（C1 规格：拒绝，不自动重排）。"""


@dataclass(frozen=True)
class Coordinates:
    """一份与 IR 配对的 Cartesian 坐标（nm）。

    ``atom_order`` 是逐原子的 ``(atomic_number, element, name)``，写文件时一并
    落盘；读取后用 :func:`check_against_ir` 与目标 IR 核对。
    """

    positions_nm: np.ndarray                 # (N,3)
    source: str                              # 如 "inpcrd:/path"、"hessian:gfn2"
    atom_order: tuple[tuple[int, str, str], ...]
    meta: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        pos = np.asarray(self.positions_nm, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"positions_nm 形状应为 (N,3)，得到 {pos.shape}")
        if len(self.atom_order) != pos.shape[0]:
            raise ValueError(
                f"atom_order 长度 {len(self.atom_order)} != 坐标原子数 {pos.shape[0]}"
            )
        object.__setattr__(self, "positions_nm", pos)


def _ir_order(ir: ForceFieldIR) -> tuple[tuple[int, str, str], ...]:
    return tuple((a.atomic_number, a.element, a.name) for a in ir.atoms)


def from_positions(
    positions_nm: np.ndarray, ir: ForceFieldIR, source: str,
    **meta: str,
) -> Coordinates:
    """从数组和 IR 构造（原子序数/元素/原子名直接取自 IR）。"""
    return Coordinates(
        positions_nm=np.asarray(positions_nm, dtype=np.float64),
        source=source,
        atom_order=_ir_order(ir),
        meta=dict(meta),
    )


def coords_to_dict(c: Coordinates) -> dict[str, Any]:
    return {
        "units": "nm",
        "source": c.source,
        "atom_order": [[z, el, name] for z, el, name in c.atom_order],
        "positions_nm": c.positions_nm.tolist(),
        "meta": dict(c.meta),
    }


def coords_from_dict(data: dict[str, Any]) -> Coordinates:
    schema = data.get("schema")
    if schema != COORDS_SCHEMA_VERSION:
        raise ValueError(
            f"坐标文件 schema 版本不符：文件为 {schema!r}，本代码期望 "
            f"{COORDS_SCHEMA_VERSION!r}"
        )
    if data.get("units") != "nm":
        raise ValueError(f"坐标单位必须为 nm，文件标注 {data.get('units')!r}")
    return Coordinates(
        positions_nm=np.array(data["positions_nm"], dtype=np.float64),
        source=str(data["source"]),
        atom_order=tuple((int(z), str(el), str(name)) for z, el, name in data["atom_order"]),
        meta=dict(data.get("meta") or {}),
    )


def to_json(c: Coordinates, path: Path) -> None:
    payload = {"schema": COORDS_SCHEMA_VERSION, **coords_to_dict(c)}
    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def from_json(path: Path) -> Coordinates:
    return coords_from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def check_against_ir(c: Coordinates, ir: ForceFieldIR) -> None:
    """C1 规格：原子序数/元素/原子名序列必须逐位一致；不一致抛错，不做任何重排。"""
    ir_order = _ir_order(ir)
    if c.atom_order != ir_order:
        raise CoordinatesMismatchError(
            f"坐标与 IR 的原子顺序不一致（原子序数/元素/原子名序列不同）："
            f"coords={list(c.atom_order)[:8]}… ir={list(ir_order)[:8]}…"
            "（不做自动重排；请确认 coords 与 ir 来自同一次运行）"
        )
