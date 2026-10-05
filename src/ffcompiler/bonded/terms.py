"""bonded 参数推导的共享数据结构（计划 T17 定义，T20/T21 复用）。"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DerivedTerm:
    """一条推导出的 bond/angle 项（计划 T17 的接口）。

    atoms 与输入顺序相同：bond (i, j)、angle (i, j, k)（j 为顶点）。
    value 为参考几何上的平衡值（r0 nm / theta0 rad），k 为 ½k 约定的
    力常数（kJ/mol/nm² 或 kJ/mol/rad²）。
    """

    atoms: tuple[int, ...]
    value: float
    k: float


__all__ = ["DerivedTerm"]
