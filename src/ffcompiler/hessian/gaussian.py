"""Gaussian fchk 读取器（MVP-1 计划 T12）。

输出统一单位的 CartesianHessian：坐标 nm、Hessian kJ/mol/nm²、
梯度 kJ/mol/nm、能量 kJ/mol（计划 §2.1）。
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.units import (
    BOHR_NM,
    HARTREE_KJ_MOL,
    HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2,
)

# fchk 里 D 幂写法（Fortran 双精度）转 E
def _to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


class _MissingFieldError(ValueError):
    """fchk 中根本没有该字段（区别于字段存在但内容残缺/格式错误）。"""


def _find_array(lines: list[str], name: str) -> tuple[int, np.ndarray]:
    """返回（头行行号, 数组）。头行形如 ``Name   R   N=  count``。

    字段不存在 → ``_MissingFieldError``；字段存在但数值不足（被截断、
    提前遇到下一个字段头或文件尾）或数值无法解析 → 普通 ``ValueError``。
    """
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped.startswith(name):
            continue
        rest = stripped[len(name):].split()
        if len(rest) >= 3 and rest[-2] in ("N=", "N"):
            count = int(rest[-1])
            values: list[float] = []
            i = idx + 1
            while len(values) < count and i < len(lines):
                tokens = lines[i].split()
                if tokens and tokens[0][:1].isalpha():
                    break  # 遇到下一个字段头：本字段被截断
                for token in tokens:
                    if len(values) == count:
                        break
                    try:
                        values.append(_to_float(token))
                    except ValueError as exc:
                        raise ValueError(
                            f"fchk 字段 {name!r} 第 {i + 1} 行含非数值 {token!r}"
                        ) from exc
                i += 1
            if len(values) < count:
                raise ValueError(
                    f"fchk 字段 {name!r} 数值不足：需要 {count} 个，读到 {len(values)} 个"
                )
            return idx, np.array(values, dtype=np.float64)
    raise _MissingFieldError(f"fchk 缺少字段 {name!r}")


def _find_scalar(lines: list[str], name: str) -> float | None:
    """标量字段值写在头行上：``Name   R   value``。"""
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith(name):
            continue
        rest = stripped[len(name):].split()
        if len(rest) >= 2 and rest[-2] in ("R", "I"):
            return _to_float(rest[-1])
    return None


def _lower_triangle_to_full(values: np.ndarray, n3: int) -> np.ndarray:
    """fchk 下三角按行展开：(0,0),(1,0),(1,1),(2,0),…，镜像成满阵。"""
    if values.size != n3 * (n3 + 1) // 2:
        raise ValueError(
            f"Hessian 数值个数 {values.size} 不等于 3N(3N+1)/2 = {n3 * (n3 + 1) // 2}"
        )
    full = np.zeros((n3, n3))
    rows, cols = np.tril_indices(n3)
    full[rows, cols] = values
    return full + full.T - np.diag(np.diag(full))


def read_fchk(path: Path) -> CartesianHessian:
    """读取 Gaussian fchk（freq 任务）为 CartesianHessian。

    读取字段：Atomic numbers (I)、Current cartesian coordinates（Bohr）、
    Cartesian Force Constants（下三角，Ha/Bohr²）、Cartesian Gradient（可选）、
    Total Energy。method 取第二行（route 行）的方法名 token。
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: 不是 Gaussian fchk 文件（行数不足）")
    route_tokens = lines[1].split()
    method = route_tokens[1] if len(route_tokens) >= 2 else lines[1].strip()

    _, atomic_numbers = _find_array(lines, "Atomic numbers")
    n_atoms = atomic_numbers.size

    _, coords_bohr = _find_array(lines, "Current cartesian coordinates")
    if coords_bohr.size != 3 * n_atoms:
        raise ValueError(
            f"坐标个数 {coords_bohr.size} != 3N = {3 * n_atoms}"
        )

    hessian_values: np.ndarray | None = None
    try:
        _, hessian_values = _find_array(lines, "Cartesian Force Constants")
    except _MissingFieldError as exc:
        raise ValueError(
            f"{path}: 缺少 Cartesian Force Constants——"
            "这不是 freq 任务的 fchk，无法提供 Hessian"
        ) from exc

    try:
        _, gradient_ha = _find_array(lines, "Cartesian Gradient")
    except _MissingFieldError:
        gradient_ha = None  # 可选字段：只有真正缺失才视为没有梯度；残缺字段照常报错
    total_energy_hartree = _find_scalar(lines, "Total Energy")

    positions_nm = coords_bohr.reshape(n_atoms, 3) * BOHR_NM
    hessian = _lower_triangle_to_full(hessian_values, 3 * n_atoms) * (
        HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2
    )
    gradient = (
        None
        if gradient_ha is None
        else gradient_ha * (HARTREE_KJ_MOL / BOHR_NM)
    )
    energy_kj = (
        None if total_energy_hartree is None else total_energy_hartree * HARTREE_KJ_MOL
    )

    return CartesianHessian(
        atomic_numbers=atomic_numbers.astype(np.int64),
        positions_nm=positions_nm,
        hessian=hessian,
        gradient=gradient,
        energy_kj=energy_kj,
        provider="gaussian-fchk",
        method=method,
        source_file=str(path),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        program_version=None,
    )
