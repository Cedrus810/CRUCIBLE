"""Sobtop param.txt 参考数据读取（MVP-1 计划 T16）。

param.txt 行格式（Sobtop 2026.1.16 实测）：
``     1        1S      2F         1.578599          2928.942          700.034``
原子 token = 编号+元素（如 ``1S``），编号 1 起算，读取时转 0 起算。
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_ATOM_TOKEN = re.compile(r"^(\d+)([A-Za-z]+)$")


@dataclass(frozen=True)
class SobtopParams:
    """Sobtop 输出的 bond/angle 刚性参数（原始单位，原子号 0 起算）。"""

    bonds: dict[tuple[int, int], tuple[float, float]]
    """(a, b)（升序）→ (键长 Å, k kJ/mol/Å²)"""
    angles: dict[tuple[int, int, int], tuple[float, float]]
    """(a, b, c)（b 为顶点，端点升序）→ (角度 deg, k kJ/mol/rad²)"""


def _parse_atom_token(token: str) -> int:
    match = _ATOM_TOKEN.match(token)
    if match is None:
        raise ValueError(f"无法解析原子 token {token!r}（应为 编号+元素，如 1S）")
    return int(match.group(1)) - 1


def read_sobtop_params(path: Path) -> SobtopParams:
    """读取 Sobtop param.txt。键和角的个数不足表头声明时报错。"""
    lines = Path(path).read_text().splitlines()

    bonds: dict[tuple[int, int], tuple[float, float]] = {}
    angles: dict[tuple[int, int, int], tuple[float, float]] = {}
    section: str | None = None
    declared: dict[str, int] = {}

    for line in lines:
        stripped = line.strip()
        match = re.match(r"Number of (bonds|angles|dihedrals):\s*(\d+)", stripped)
        if match:
            section = match.group(1)
            declared[section] = int(match.group(2))
            continue
        if not stripped or section is None or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if section == "bonds" and len(tokens) == 6:
            a, b = _parse_atom_token(tokens[1]), _parse_atom_token(tokens[2])
            key = (min(a, b), max(a, b))
            bonds[key] = (float(tokens[3]), float(tokens[4]))
        elif section == "angles" and len(tokens) == 7:
            a, b, c = (
                _parse_atom_token(tokens[1]),
                _parse_atom_token(tokens[2]),
                _parse_atom_token(tokens[3]),
            )
            if a > c:
                a, c = c, a
            angles[(a, b, c)] = (float(tokens[4]), float(tokens[5]))

    for kind, data in (("bonds", bonds), ("angles", angles)):
        want = declared.get(kind)
        if want is not None and len(data) != want:
            raise ValueError(
                f"{path}: {kind} 声明 {want} 个，实际解析出 {len(data)} 个"
            )
    return SobtopParams(bonds=bonds, angles=angles)
