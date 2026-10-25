"""CP2K .out 的 vibrational analysis Hessian 读取器（决策 0023，吸收自 Sobtop）。

CP2K（``RUN_TYPE VIBRATIONAL_ANALYSIS``）输出中的分页矩阵段：

::

    VIB| Hessian in cartesian coordinates

                        1            2            3            4            5
         1  Cl       1.136309     0.000643    -0.001125    -0.267734    -0.088337
         2  Cl       0.000642     1.136802    ...
        （每页 5 列：页头为 1-based 列号行，页体为 3N 行「行号 元素 值×5」，
         页与页重复全部 3N 行直至 3N 列写满）

- 单位 E_h/a₀²（与 fchk/orca/xtb 同体系），读入时换算为 kJ/mol/nm²；
- 元素符号取自矩阵行（VIB 段自带）；
- 参考几何取文件中**最后一次** ``MODULE QUICKSTEP:  ATOMIC COORDINATES IN
  angstrom`` 坐标表（freq 用的最终结构）；没有该段（例如纯数值输出被裁剪）
  时抛 ValueError——调用方也可显式传 ``positions_nm`` 覆盖。

与 Sobtop ``loadcp2khess_out`` 行为对齐（它同样只从该文件读 Hessian，几何
来自输入结构文件）。
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np

from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.units import (
    ANGSTROM_NM,
    HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2,
)

_HESSIAN_HEADER = "VIB| Hessian in cartesian coordinates"
_COORD_HEADER = re.compile(r"MODULE QUICKSTEP:\s+ATOMIC COORDINATES IN angstrom")


def _to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _element_to_z(token: str) -> int:
    from ffcompiler.hessian.orca import _element_to_z as _z

    return _z(token)


def _is_int(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def _read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def _parse_paged_matrix(path: Path, lines: list[str], start: int) -> tuple[np.ndarray, list[str]]:
    """从 Hessian 表头下一行起解析分页矩阵，返回 (3N×3N 矩阵, 元素符号列表)。"""
    n_total: int | None = None
    matrix: np.ndarray | None = None
    elements: list[str] = []
    pos = start
    while pos < len(lines) and not lines[pos].strip():
        pos += 1  # 表头后的空行
    pages: list[tuple[list[int], list[tuple[int, str, list[float]]]]] = []
    while pos < len(lines):
        header = lines[pos].split()
        if not header or not all(_is_int(t) for t in header):
            break  # 矩阵结束（后续是别的 VIB 段）
        cols = [int(t) - 1 for t in header]
        page_rows: list[tuple[int, str, list[float]]] = []
        pos += 1
        while pos < len(lines):
            tokens = lines[pos].split()
            if (
                len(tokens) < 3
                or not _is_int(tokens[0])
                or not tokens[1].isalpha()
                or len(tokens) != 2 + len(cols)
            ):
                break
            page_rows.append(
                (int(tokens[0]) - 1, tokens[1], [_to_float(t) for t in tokens[2:]])
            )
            pos += 1
        if not page_rows:
            raise ValueError(f"{path}: Hessian 第 {pos} 页没有数据行")
        pages.append((cols, page_rows))
        # 跳过页间空行
        while pos < len(lines) and not lines[pos].strip():
            pos += 1
    if not pages:
        raise ValueError(f"{path}: VIB| Hessian 段之后没有矩阵页")
    n_total = pages[0][1][-1][0] + 1
    if any(r[0] + 1 != n_total for _, rows in pages for r in rows):
        pass  # 行号连续性在下面的覆盖检查中体现
    matrix = np.full((n_total, n_total), np.nan)
    seen_rows: set[int] = set()
    for cols, rows in pages:
        if len(rows) != n_total:
            raise ValueError(
                f"{path}: Hessian 列页 {cols[:2]}… 应有 {n_total} 行，实际 {len(rows)}"
            )
        for r, elem, values in rows:
            matrix[r, cols] = values
            seen_rows.add(r)
            if not elements:
                pass
    # 元素符号：取第一页的行序（每页重复同样 3N 行）
    first_rows = pages[0][1]
    elements = [e for _, e, _ in first_rows]
    if np.isnan(matrix).any():
        raise ValueError(f"{path}: VIB| Hessian 分页矩阵不完整（有元素未出现）")
    return matrix, elements


def _parse_last_coordinates(path: Path, lines: list[str]) -> np.ndarray | None:
    """最后一次 QUICKSTEP 坐标表（Å）。

    表头 ``Atom Kind Element X Y Z Z(eff) Mass``，每原子一行 9 列：
    ``idx kind element Z x y z z_eff mass``（坐标取第 5–7 列）。
    """
    last: int | None = None
    for idx, line in enumerate(lines):
        if _COORD_HEADER.search(line):
            last = idx
    if last is None:
        return None
    pos = last + 1
    while pos < len(lines) and not lines[pos].strip():
        pos += 1  # 表头后的空行
    if pos < len(lines) and lines[pos].split()[:1] == ["Atom"]:
        pos += 1  # 列名行 "Atom Kind Element X Y Z Z(eff) Mass"
    coords: list[tuple[str, tuple[float, float, float]]] = []
    while pos < len(lines):
        tokens = lines[pos].split()
        if (
            len(tokens) != 9
            or not _is_int(tokens[0])
            or not _is_int(tokens[1])
            or not tokens[2].isalpha()
            or not _is_int(tokens[3])
        ):
            break
        try:
            xyz = tuple(_to_float(t) for t in tokens[4:7])
        except ValueError:
            break
        coords.append((tokens[2], xyz))
        pos += 1
    if not coords:
        return None
    return np.array([c for _, c in coords], dtype=np.float64)


def read_cp2k_hess(
    path: Path,
    *,
    positions_nm: np.ndarray | None = None,
) -> CartesianHessian:
    """读取 CP2K .out 的 VIB Hessian 为 CartesianHessian。

    缺 ``VIB| Hessian in cartesian coordinates`` 段、矩阵不完整、坐标原子数
    与矩阵维度不一致 → ValueError。``positions_nm`` 可显式覆盖文件内坐标
    （调用方的结构文件才是权威几何时使用）。
    """
    path = Path(path)
    lines = _read_lines(path)
    header_idx = [i for i, ln in enumerate(lines) if ln.strip() == _HESSIAN_HEADER]
    if not header_idx:
        raise ValueError(
            f"{path}: 找不到 '{_HESSIAN_HEADER}'——确认 RUN_TYPE 为 "
            "VIBRATIONAL_ANALYSIS 且打印级别包含 VIB 段"
        )
    matrix, elements = _parse_paged_matrix(path, lines, header_idx[-1] + 1)
    # CP2K 的 VIB Hessian 是前向差分（VIB 段只有 +DX 复制品），raw 矩阵天然
    # 非对称；对角化/投影前对称化是标准处理，原始不对称度记入 meta。
    raw_asym = float(np.abs(matrix - matrix.T).max())
    matrix = 0.5 * (matrix + matrix.T)
    n_atoms, rem = divmod(matrix.shape[0], 3)
    if rem:
        raise ValueError(f"{path}: Hessian 维度 {matrix.shape[0]} 不是 3 的倍数")
    if len(elements) != 3 * n_atoms:
        raise ValueError(
            f"{path}: 矩阵行元素数 {len(elements)} 与维度 {3 * n_atoms} 不一致"
        )
    # 行 = 坐标分量（3 行/原子），元素符号按 r//3 折叠并校验同一原子三行一致
    per_atom = [elements[3 * k] for k in range(n_atoms)]
    for k, sym in enumerate(per_atom):
        if elements[3 * k + 1] != sym or elements[3 * k + 2] != sym:
            raise ValueError(
                f"{path}: 第 {k + 1} 个原子的三个坐标行元素符号不一致："
                f"{elements[3 * k:3 * k + 3]}"
            )
    atomic_numbers = np.array([_element_to_z(e) for e in per_atom], dtype=np.int64)
    if positions_nm is not None:
        pos = np.asarray(positions_nm, dtype=np.float64)
    else:
        coords_a = _parse_last_coordinates(path, lines)
        if coords_a is None:
            raise ValueError(
                f"{path}: 文件内没有可用的 QUICKSTEP 坐标表，且未显式传 "
                "positions_nm——请传入结构文件的几何"
            )
        if coords_a.shape[0] != n_atoms:
            raise ValueError(
                f"{path}: 坐标原子数 {coords_a.shape[0]} 与 Hessian 维度 "
                f"{n_atoms} 不一致；可用 positions_nm 显式覆盖"
            )
        pos = coords_a * ANGSTROM_NM

    return CartesianHessian(
        atomic_numbers=atomic_numbers,
        positions_nm=pos,
        hessian=matrix * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2,
        gradient=None,
        energy_kj=None,
        provider="cp2k-out",
        method="unknown",  # .out 的 VIB 段不记录方法名
        source_file=str(path),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        meta={
            "note": "CP2K forward-difference VIB Hessian, symmetrized on read",
            "raw_max_asymmetry_eh_bohr2": repr(raw_asym),
        },
    )
