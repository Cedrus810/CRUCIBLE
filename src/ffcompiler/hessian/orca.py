"""ORCA .hess 读取器（路线图 D3；2026-10-17 施工，同日按真实样例重写）。

输出统一单位的 CartesianHessian：坐标 nm、Hessian kJ/mol/nm²（hessian/base
的内部约定）。梯度与能量不在 .hess 里（None）。

**格式依据**：ORCA 6.1.1 实际写出的文件
（``experiments/20261017_orca_hess_format/water.hess``，HF/STO-3G 水；固件
``tests/fixtures/hessian/orca/water.hess``）。布局：

- 文件头 ``$orca_hessian_file``（无内容），之后若干 ``$节``，``$end`` 结束；
- ``$hessian``：第一行是维度 3N；随后按列分块的**满阵**（Hartree/Bohr²），
  每块先是一行列号表头，再是 3N 行"行号 + 该块各列的值"；
- ``$atoms``：第一行原子数 N；之后每行"元素 质量 x y z"（坐标 Bohr，
  ORCA 已平移到质心）。

未知节跳过；块结构、维度、行号/列号不符一律显式报错，不猜。质量记进
``meta["orca_masses"]``（ORCA 计算频率用的质量，逗号分隔），不改变
CartesianHessian 的接口。
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from rdkit import Chem

from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.units import (
    BOHR_NM,
    HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2,
)

FILE_HEADER = "orca_hessian_file"


def _to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _sections(path: Path) -> tuple[str | None, dict[str, list[str]]]:
    """按 ``$name`` 分节；返回 (第一个节名, {节名: 该节的非空行（已 strip）})。"""
    lines_by_section: dict[str, list[str]] = {}
    first: str | None = None
    current: str | None = None
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("$"):
            name = stripped[1:].split()[0] if len(stripped) > 1 else ""
            if name == "end":
                current = None
                continue
            first = first or name
            current = name
            lines_by_section.setdefault(current, [])
            continue
        if current is not None:
            lines_by_section[current].append(stripped)
    return first, lines_by_section


def _element_to_z(token: str) -> int:
    if token.isdigit():
        return int(token)
    try:
        return Chem.GetPeriodicTable().GetAtomicNumber(token.capitalize())
    except Exception as exc:  # noqa: BLE001 - rdkit 对未知元素抛 RuntimeError/KeyError
        raise ValueError(f"无法识别的元素符号 {token!r}") from exc


def _is_int(token: str) -> bool:
    return token.lstrip("+-").isdigit()


def _parse_blocked_matrix(path: Path, lines: list[str]) -> np.ndarray:
    """``$hessian`` 节：维度行 + 按列分块的满阵（每块：列号表头 + n 行"行号 值…"）。"""
    if not lines or not _is_int(lines[0]):
        raise ValueError(f"{path}: $hessian 节的第一行应是维度 3N")
    n = int(lines[0])
    if n <= 0 or n % 3:
        raise ValueError(f"{path}: $hessian 维度 {n} 不是 3 的正整数倍")
    matrix = np.full((n, n), np.nan)
    pos = 1
    while pos < len(lines):
        header = lines[pos].split()
        if not all(_is_int(t) for t in header):
            raise ValueError(f"{path}: $hessian 第 {pos} 行应是列号表头，得到 {lines[pos]!r}")
        cols = [int(t) for t in header]
        if any(not 0 <= c < n for c in cols):
            raise ValueError(f"{path}: $hessian 列号 {cols} 超出 0..{n - 1}")
        if pos + n >= len(lines):
            raise ValueError(f"{path}: $hessian 列块 {cols} 不足 {n} 行")
        for r in range(n):
            tokens = lines[pos + 1 + r].split()
            if len(tokens) != len(cols) + 1 or not _is_int(tokens[0]) or int(tokens[0]) != r:
                raise ValueError(
                    f"{path}: $hessian 列块 {cols} 第 {r} 行格式不符：{lines[pos + 1 + r]!r}")
            matrix[r, cols] = [_to_float(t) for t in tokens[1:]]
        pos += 1 + n
    if np.isnan(matrix).any():
        raise ValueError(f"{path}: $hessian 列块不完整（有列未出现）")
    return matrix


def read_orca_hess(path: Path) -> CartesianHessian:
    """读取 ORCA .hess 为 CartesianHessian（格式见模块说明）。

    文件头不是 ``$orca_hessian_file``、缺 ``$hessian`` / ``$atoms``、维度与
    原子数不符、行列结构不符 → ValueError。
    """
    path = Path(path)
    first, sections = _sections(path)
    if first != FILE_HEADER:
        raise ValueError(
            f"{path}: 第一个节是 {first!r}，不是 ${FILE_HEADER}——不是 ORCA .hess 文件")
    if "hessian" not in sections:
        raise ValueError(f"{path}: 缺少 $hessian 节")
    hessian = _parse_blocked_matrix(path, sections["hessian"])
    n3 = hessian.shape[0]

    atom_lines = sections.get("atoms")
    if atom_lines is None:
        raise ValueError(f"{path}: 缺少 $atoms 节，无法得到元素与参考几何")
    if not atom_lines or not _is_int(atom_lines[0]):
        raise ValueError(f"{path}: $atoms 节的第一行应是原子数")
    n_atoms = int(atom_lines[0])
    if 3 * n_atoms != n3:
        raise ValueError(
            f"{path}: $atoms 的原子数 {n_atoms} 与 $hessian 维度 {n3} 不一致（应为 3N）")
    rows = [line.split() for line in atom_lines[1:]]
    if len(rows) != n_atoms or any(len(r) != 5 for r in rows):
        raise ValueError(
            f"{path}: $atoms 应有 {n_atoms} 行「元素 质量 x y z」，实际 {len(rows)} 行或列数不符")
    atomic_numbers = np.array([_element_to_z(r[0]) for r in rows], dtype=np.int64)
    masses = [_to_float(r[1]) for r in rows]
    coords_bohr = np.array([[_to_float(t) for t in r[2:5]] for r in rows], dtype=np.float64)

    return CartesianHessian(
        atomic_numbers=atomic_numbers,
        positions_nm=coords_bohr * BOHR_NM,
        hessian=hessian * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2,
        gradient=None,
        energy_kj=None,
        provider="orca-hess",
        method="unknown",  # .hess 不记录计算方法
        source_file=str(path),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        program_version=None,
        meta={"orca_masses": ",".join(repr(m) for m in masses)},
    )
