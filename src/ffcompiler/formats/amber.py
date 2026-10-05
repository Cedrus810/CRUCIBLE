"""原生 AMBER prmtop / rst7 读写（决策 0015：不依赖 ParmEd）。

prmtop 由 ``%FLAG NAME`` / ``%FORMAT(fmt)`` / 定宽数据行组成。本模块只做格式层：
按 Fortran 格式（``20a4``、``5E16.8``、``10I8``、``1a80``…）解析、写出各 FLAG 的
原始数组，并按 POINTERS 核对数组长度。物理含义（单位、½k、1-4 语义）由
``importers.amber_prmtop`` 与 ``exporters.amber`` 处理。
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

POINTER_NAMES = (
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH", "MPHIA",
    "NHPARM", "NPARM", "NNB", "NRES", "NBONA", "NTHETA", "NPHIA", "NUMBND",
    "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT", "NBPER", "NGPER", "NDPER",
    "MBPER", "MGPER", "MDPER", "IFBOX", "NMXRS", "IFCAP", "NUMEXTRA", "NCOPY",
)

_FORMAT_RE = re.compile(r"^\s*(\d*)\s*([aAiIeEfF])\s*(\d+)(?:\.(\d+))?\s*$")

REQUIRED_FLAGS = (
    "POINTERS", "ATOM_NAME", "CHARGE", "MASS", "ATOM_TYPE_INDEX",
    "NUMBER_EXCLUDED_ATOMS", "NONBONDED_PARM_INDEX", "RESIDUE_LABEL", "RESIDUE_POINTER",
    "BOND_FORCE_CONSTANT", "BOND_EQUIL_VALUE", "ANGLE_FORCE_CONSTANT", "ANGLE_EQUIL_VALUE",
    "DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PERIODICITY", "DIHEDRAL_PHASE",
    "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF",
    "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN",
    "ANGLES_INC_HYDROGEN", "ANGLES_WITHOUT_HYDROGEN",
    "DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN",
    "EXCLUDED_ATOMS_LIST", "AMBER_ATOM_TYPE",
)


class PrmtopError(ValueError):
    """prmtop / rst7 格式错误（缺段、长度不符、无法解析的字段）。"""


@dataclass(frozen=True)
class FortranFormat:
    per_line: int
    kind: str      # "a" | "i" | "e" | "f"
    width: int
    decimals: int | None

    @classmethod
    def parse(cls, text: str) -> "FortranFormat":
        m = _FORMAT_RE.match(text)
        if m is None:
            raise PrmtopError(f"无法解析的 FORMAT: {text!r}")
        per_line = int(m[1]) if m[1] else 1
        return cls(per_line, m[2].lower(), int(m[3]), int(m[4]) if m[4] else None)

    def __str__(self) -> str:
        kind = {"a": "a", "i": "I", "e": "E", "f": "F"}[self.kind]
        tail = f".{self.decimals}" if self.decimals is not None else ""
        return f"{self.per_line}{kind}{self.width}{tail}"

    def format_value(self, value) -> str:
        if self.kind == "a":
            return f"{str(value):<{self.width}.{self.width}}"
        if self.kind == "i":
            return f"{int(value):>{self.width}d}"
        if self.kind == "e":
            return f"{float(value):>{self.width}.{self.decimals}E}"
        return f"{float(value):>{self.width}.{self.decimals}f}"

    def parse_value(self, text: str, flag: str):
        if self.kind == "a":
            return text.strip()
        try:
            if self.kind == "i":
                return int(text)
            return float(text.replace("D", "E").replace("d", "e"))
        except ValueError as exc:
            raise PrmtopError(f"{flag}: 无法解析字段 {text!r}") from exc


@dataclass
class Section:
    fmt: str
    values: list


@dataclass
class PrmtopFile:
    version: str
    flags: dict[str, Section] = field(default_factory=dict)

    def __getitem__(self, flag: str) -> list:
        return self.flags[flag].values

    def pointer(self, name: str) -> int:
        return int(self.flags["POINTERS"].values[POINTER_NAMES.index(name)])

    def set(self, flag: str, fmt: str, values) -> None:
        self.flags[flag] = Section(fmt, list(values))


def _parse_data(lines: list[str], fmt: FortranFormat, flag: str) -> list:
    values = []
    for line in lines:
        text = line.rstrip("\n")
        if fmt.kind == "a":
            text = text.rstrip()
        n_fields = math.ceil(len(text) / fmt.width) if text else 0
        for k in range(min(n_fields, fmt.per_line)):
            chunk = text[k * fmt.width:(k + 1) * fmt.width]
            if fmt.kind != "a" and not chunk.strip():
                raise PrmtopError(f"{flag}: 空的数值字段")
            values.append(fmt.parse_value(chunk, flag))
    return values


def _expected_lengths(prm: PrmtopFile) -> dict[str, int]:
    p = prm.pointer
    natom, ntypes = p("NATOM"), p("NTYPES")
    lengths = {
        "ATOM_NAME": natom, "CHARGE": natom, "MASS": natom, "ATOM_TYPE_INDEX": natom,
        "NUMBER_EXCLUDED_ATOMS": natom, "AMBER_ATOM_TYPE": natom,
        "ATOMIC_NUMBER": natom, "TREE_CHAIN_CLASSIFICATION": natom,
        "JOIN_ARRAY": natom, "IROTAT": natom, "RADII": natom, "SCREEN": natom,
        "NONBONDED_PARM_INDEX": ntypes * ntypes,
        "LENNARD_JONES_ACOEF": ntypes * (ntypes + 1) // 2,
        "LENNARD_JONES_BCOEF": ntypes * (ntypes + 1) // 2,
        "RESIDUE_LABEL": p("NRES"), "RESIDUE_POINTER": p("NRES"),
        "BOND_FORCE_CONSTANT": p("NUMBND"), "BOND_EQUIL_VALUE": p("NUMBND"),
        "ANGLE_FORCE_CONSTANT": p("NUMANG"), "ANGLE_EQUIL_VALUE": p("NUMANG"),
        "DIHEDRAL_FORCE_CONSTANT": p("NPTRA"), "DIHEDRAL_PERIODICITY": p("NPTRA"),
        "DIHEDRAL_PHASE": p("NPTRA"), "SCEE_SCALE_FACTOR": p("NPTRA"),
        "SCNB_SCALE_FACTOR": p("NPTRA"), "SOLTY": p("NATYP"),
        "BONDS_INC_HYDROGEN": 3 * p("NBONH"), "BONDS_WITHOUT_HYDROGEN": 3 * p("MBONA"),
        "ANGLES_INC_HYDROGEN": 4 * p("NTHETH"), "ANGLES_WITHOUT_HYDROGEN": 4 * p("MTHETA"),
        "DIHEDRALS_INC_HYDROGEN": 5 * p("NPHIH"), "DIHEDRALS_WITHOUT_HYDROGEN": 5 * p("MPHIA"),
        "EXCLUDED_ATOMS_LIST": p("NNB"),
        "HBOND_ACOEF": p("NPHB"), "HBOND_BCOEF": p("NPHB"), "HBCUT": p("NPHB"),
    }
    return lengths


def check_lengths(prm: PrmtopFile) -> None:
    for flag in REQUIRED_FLAGS:
        if flag not in prm.flags:
            raise PrmtopError(f"prmtop 缺少 %FLAG {flag}")
    if len(prm["POINTERS"]) < 31:
        raise PrmtopError(f"POINTERS 只有 {len(prm['POINTERS'])} 个值（至少 31 个）")
    for flag, n in _expected_lengths(prm).items():
        if flag in prm.flags and len(prm[flag]) != n:
            raise PrmtopError(f"%FLAG {flag} 有 {len(prm[flag])} 个值，POINTERS 要求 {n}")


def read_prmtop(path: Path) -> PrmtopFile:
    lines = Path(path).read_text().splitlines()
    if not lines:
        raise PrmtopError(f"{path}: 空文件")
    version = lines[0].rstrip() if lines[0].startswith("%VERSION") else ""
    prm = PrmtopFile(version=version)
    i = 1 if version else 0
    while i < len(lines):
        line = lines[i]
        if not line.startswith("%FLAG"):
            if line.strip() and not line.startswith("%COMMENT"):
                raise PrmtopError(f"{path}:{i + 1}: 期望 %FLAG，得到 {line[:40]!r}")
            i += 1
            continue
        flag = line[5:].strip()
        i += 1
        while i < len(lines) and lines[i].startswith("%COMMENT"):
            i += 1
        if i >= len(lines) or not lines[i].startswith("%FORMAT"):
            raise PrmtopError(f"%FLAG {flag} 后缺少 %FORMAT")
        fmt_text = lines[i].strip()[len("%FORMAT("):-1]
        fmt = FortranFormat.parse(fmt_text)
        i += 1
        start = i
        while i < len(lines) and not lines[i].startswith("%"):
            i += 1
        if flag in prm.flags:
            raise PrmtopError(f"%FLAG {flag} 重复")
        prm.flags[flag] = Section(str(fmt), _parse_data(lines[start:i], fmt, flag))
    check_lengths(prm)
    return prm


def write_prmtop(prm: PrmtopFile, path: Path) -> None:
    check_lengths(prm)
    out = [(prm.version or "%VERSION  VERSION_STAMP = V0001.000").ljust(80)]
    for flag, section in prm.flags.items():
        fmt = FortranFormat.parse(section.fmt)
        out.append(f"%FLAG {flag}".ljust(80))
        out.append(f"%FORMAT({section.fmt})".ljust(80))
        values = section.values
        if not values:
            out.append("")
            continue
        for k in range(0, len(values), fmt.per_line):
            out.append("".join(fmt.format_value(v) for v in values[k:k + fmt.per_line]))
    Path(path).write_text("\n".join(out) + "\n")


def read_rst7(path: Path) -> np.ndarray:
    """ASCII rst7 / inpcrd → 坐标 (N, 3)，单位 Å。速度、盒子行忽略。"""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise PrmtopError(f"{path}: rst7 太短")
    try:
        natom = int(lines[1].split()[0])
    except (IndexError, ValueError) as exc:
        raise PrmtopError(f"{path}: 第二行不是原子数") from exc
    values: list[float] = []
    for line in lines[2:]:
        for k in range(0, len(line.rstrip()), 12):
            chunk = line[k:k + 12]
            if chunk.strip():
                values.append(float(chunk))
        if len(values) >= 3 * natom:
            break
    if len(values) < 3 * natom:
        raise PrmtopError(f"{path}: coordinates 只有 {len(values)} 个值，需要 {3 * natom}")
    return np.array(values[:3 * natom], dtype=np.float64).reshape(natom, 3)


def write_rst7(coords_angstrom: np.ndarray, path: Path, *, title: str = "ffcompiler") -> None:
    coords = np.asarray(coords_angstrom, dtype=np.float64).reshape(-1, 3)
    if not np.all(np.isfinite(coords)) or np.abs(coords).max(initial=0.0) >= 1e4:
        raise PrmtopError("rst7 坐标必须有限且 |x| < 10000 Å（12.7f 定宽）")
    flat = coords.ravel()
    out = [title[:80], f"{len(coords):6d}"]
    for k in range(0, len(flat), 6):
        out.append("".join(f"{v:12.7f}" for v in flat[k:k + 6]))
    Path(path).write_text("\n".join(out) + "\n")
