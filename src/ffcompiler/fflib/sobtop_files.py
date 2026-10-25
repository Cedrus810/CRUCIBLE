"""Sobtop 格式力场数据文件的读取与匹配（决策 0022，吸收自 Sobtop 逆向）。

支持三种随 Sobtop 发行、也可由用户自编的数据文件：

- ``LJ_param.dat``：原子类型的 LJ 参数（vdW 半径 r0/2 [Å]、井深 ε [kcal/mol]），
  出厂内容为 UFF 全元素（``UF_*``，Rappe 1992）+ GAFF 1.81 + parm19 摘录；
- ``bonded_param.dat``：成键参数库，字段 ``$BOND/$ANGLE/$DIHEDRAL/$HARM_DIH/$IMPROPER``，
  AMBER 约定（键/角 E = K·Δq²，K 为 GROMACS ½k 约定的一半；扭转带 IDIVF）；
- ``bondcrit.dat``：按元素对的成键距离窗口覆盖（下限/上限 [Å]），
  注意是绝对距离区间，不是对半径判据的缩放。

匹配语义按逆向还原实现（见 docs/theory/sobtop_reverse_engineering.md §7）：

- 先试全具体类型（正反两个方向），再试含 ``X`` 通配的条目，通配越少越优先；
- ``$DIHEDRAL`` 的势垒先除以 IDIVF，同类型多条目全部保留（多项 Fourier）；
- 负周期性按 cos 偶函数等价改写：pn → |pn|，phase → −phase（mod 360）；
  phase 为 0 时数值不变（Sobtop 输出即此情形）；
- improper 的第三列（中心）类型集合决定哪些原子可生成 improper（数据驱动，无硬编码）。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..units import KCAL_KJ

__all__ = [
    "LJEntry",
    "BondEntry",
    "AngleEntry",
    "DihedralEntry",
    "ImproperEntry",
    "ForceFieldLibrary",
    "BondCriterion",
    "PeriodicTerm",
    "parse_lj_param",
    "parse_bonded_param",
    "parse_bondcrit",
    "lj_to_sigma_eps",
    "improper_center_types",
    "find_bond",
    "find_angle",
    "find_dihedrals",
    "find_improper",
    "find_lj",
    "GUESS_BOND_K_KJ_MOL_NM2",
    "GUESS_ANGLE_K_KJ_MOL_RAD2",
    "GUESS_IMPROPER_BARRIER_KCAL",
    "GUESS_IMPROPER_PHASE_DEG",
    "GUESS_IMPROPER_PN",
]

# Sobtop 缺参"乱猜"写死值（逆向实证，FeCl2 lazy / SiO2 等 itp 逐字核对）
GUESS_BOND_K_KJ_MOL_NM2 = 2.5e5      # r0 取当前距离
GUESS_ANGLE_K_KJ_MOL_RAD2 = 500.0    # a0 取当前角
GUESS_IMPROPER_BARRIER_KCAL = 1.1    # 缺 improper 用 X-X-ca-ha 等价值
GUESS_IMPROPER_PHASE_DEG = 180.0
GUESS_IMPROPER_PN = 2.0

_X = "X"


@dataclass(frozen=True)
class LJEntry:
    """LJ_param.dat 一行：类型、r0/2 [Å]、ε [kcal/mol]、描述。"""

    name: str
    r_half_a: float
    eps_kcal: float
    description: str = ""


@dataclass(frozen=True)
class BondEntry:
    """$BOND：t1-t2，K [kcal/mol/Å²]（AMBER K=k/2 约定），r0 [Å]。"""

    types: tuple[str, ...]
    k_kcal_a2: float
    r0_a: float

    @property
    def n_wildcards(self) -> int:
        return sum(t.strip() == _X for t in self.types)


@dataclass(frozen=True)
class AngleEntry:
    """$ANGLE：K [kcal/mol/rad²]（AMBER 约定），θ0 [°]。"""

    types: tuple[str, ...]
    k_kcal_rad2: float
    theta0_deg: float

    @property
    def n_wildcards(self) -> int:
        return sum(t.strip() == _X for t in self.types)


@dataclass(frozen=True)
class DihedralEntry:
    """$DIHEDRAL：IDIVF、势垒 [kcal/mol]、相位 [°]、周期性（可负）。"""

    types: tuple[str, ...]
    idivf: int
    barrier_kcal: float
    phase_deg: float
    pn: float

    @property
    def n_wildcards(self) -> int:
        return sum(t.strip() == _X for t in self.types)


@dataclass(frozen=True)
class ImproperEntry:
    """$IMPROPER：势垒 [kcal/mol]、相位 [°]、周期性。第三列是 sp2 中心。"""

    types: tuple[str, ...]
    barrier_kcal: float
    phase_deg: float
    pn: float


@dataclass(frozen=True)
class PeriodicTerm:
    """一条展开后的周期扭转/improper 项（GROMACS funct 9/4 语义）。"""

    types: tuple[str, ...]
    kd_kj: float
    phase_deg: float
    pn: float


@dataclass
class BondCriterion:
    """bondcrit.dat 的元素对距离窗口 [Å]（含边界成键）。"""

    element1: str
    element2: str
    lower_a: float
    upper_a: float


@dataclass
class ForceFieldLibrary:
    """bonded_param.dat 解析结果（按字段分组，保持文件顺序）。"""

    bonds: list[BondEntry] = field(default_factory=list)
    angles: list[AngleEntry] = field(default_factory=list)
    dihedrals: list[DihedralEntry] = field(default_factory=list)
    impropers: list[ImproperEntry] = field(default_factory=list)


def _parse_fields(line: str) -> tuple[tuple[str, ...], list[str]]:
    """行 → (类型元组, 其余数值字段字符串列表)。

    类型 token 从行首连续的非数值字段拼接（Sobtop 匹配忽略类型内空格）。
    """
    tokens = line.split()
    type_token = tokens[0]
    idx = 1
    while idx < len(tokens):
        try:
            float(tokens[idx])
        except ValueError:
            type_token += tokens[idx]
            idx += 1
        else:
            break
    parts = tuple(p.strip() for p in type_token.split("-") if p.strip())
    return parts, tokens[idx:]


def _data_lines(path: str | Path) -> list[str]:
    out: list[str] = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s and not s.startswith(";"):
            out.append(line)
    return out


def parse_lj_param(path: str | Path) -> list[LJEntry]:
    """解析 LJ_param.dat（'Ctest 1.9255 0.105 ;描述' 形式的行）。"""
    entries: list[LJEntry] = []
    for line in _data_lines(path):
        body = line.split(";", 1)
        fields = body[0].split()
        if len(fields) < 3:
            continue
        desc = body[1].strip() if len(body) > 1 else ""
        entries.append(LJEntry(fields[0], float(fields[1]), float(fields[2]), desc))
    return entries


def lj_to_sigma_eps(r_half_a: float, eps_kcal: float) -> tuple[float, float]:
    """库存 (r0/2 [Å], ε [kcal/mol]) → GROMACS (σ [nm], ε [kJ/mol])。

    σ = 2·R·2^(-1/6)·0.1 nm；ε ×4.184。实测 UF_C (1.9255, 0.105)
    → (0.3430851, 0.43932)，与 Sobtop itp 输出逐位一致。
    """
    sigma_nm = 2.0 * r_half_a * 2.0 ** (-1.0 / 6.0) * 0.1
    return sigma_nm, eps_kcal * KCAL_KJ


def parse_bonded_param(path: str | Path) -> ForceFieldLibrary:
    """解析 bonded_param.dat 的五个字段；未知字段跳过，重复字段后者追加。"""
    lib = ForceFieldLibrary()
    field_name: str | None = None
    for line in _data_lines(path):
        token = line.split()[0]
        upper = token.upper()
        if upper.startswith("$") and not upper.startswith("$END"):
            field_name = upper
            continue
        if upper.startswith("$END"):
            field_name = None
            continue
        types, fields = _parse_fields(line)
        if field_name == "$BOND" and len(fields) >= 2:
            lib.bonds.append(BondEntry(types, float(fields[0]), float(fields[1])))
        elif field_name == "$ANGLE" and len(fields) >= 2:
            lib.angles.append(
                AngleEntry(types, float(fields[0]), float(fields[1]))
            )
        elif field_name == "$DIHEDRAL" and len(fields) >= 4:
            lib.dihedrals.append(
                DihedralEntry(
                    types,
                    int(round(float(fields[0]))),
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                )
            )
        elif field_name == "$HARM_DIH":
            continue  # 出厂为空；谐振扭转库目前只读不匹配
        elif field_name == "$IMPROPER" and len(fields) >= 3:
            lib.impropers.append(
                ImproperEntry(types, float(fields[0]), float(fields[1]), float(fields[2]))
            )
    return lib


def parse_bondcrit(path: str | Path) -> list[BondCriterion]:
    """解析 bondcrit.dat：'元素1 元素2 下限 上限'（Å 距离窗口）。"""
    out: list[BondCriterion] = []
    for line in _data_lines(path):
        parts = line.split()
        if len(parts) == 4 and not parts[0].isdigit():
            try:
                out.append(
                    BondCriterion(
                        parts[0], parts[1], float(parts[2]), float(parts[3])
                    )
                )
            except ValueError:
                continue
    return out


def improper_center_types(impropers: list[ImproperEntry]) -> set[str]:
    """中心类型集合 = $IMPROPER 条目第三列去重。

    Sobtop 的 improper 感知完全数据驱动：配位数 3 且类型在该集合中的原子
    才生成 improper（GAFF ca/c/c2/n/na/n2 + AMBER 大写集来自出厂库）。
    """
    return {e.types[2].strip() for e in impropers if len(e.types) == 4}


def _matches(pattern: tuple[str, ...], types: tuple[str, ...]) -> bool:
    if len(pattern) != len(types):
        return False
    return all(p.strip() == _X or p.strip() == t.strip() for p, t in zip(pattern, types))


def _find_wildcard[T](
    entries: list[T],
    types: tuple[str, ...],
    getattr_types,
    getattr_nwild,
) -> list[T]:
    """先按特异性分层（通配数少者优先），层内正序先于反序。

    与 Sobtop 实证一致：oh-c3-c3-hc 命中反序的精确条目 hc-c3-c3-oh，
    而不是正序就能匹配的 X-c3-c3-X 通配。
    """
    for n_wild in range(len(types)):
        best: list[T] = []
        for direction in (types, tuple(reversed(types))):
            for entry in entries:
                if getattr_nwild(entry) == n_wild and _matches(
                    getattr_types(entry), direction
                ):
                    best.append(entry)
            if best:
                return best
    return []


def find_bond(lib: ForceFieldLibrary, t1: str, t2: str) -> BondEntry | None:
    """精确条目（正反序）优先于含 X 的条目；返回首条匹配。"""
    hit = _find_wildcard(lib.bonds, (t1, t2), lambda e: e.types, lambda e: e.n_wildcards)
    return hit[0] if hit else None


def find_angle(lib: ForceFieldLibrary, t1: str, t2: str, t3: str) -> AngleEntry | None:
    hit = _find_wildcard(
        lib.angles, (t1, t2, t3), lambda e: e.types, lambda e: e.n_wildcards
    )
    return hit[0] if hit else None


def _expand_periodic(
    types: tuple[str, ...], barrier_kcal: float, phase_deg: float, pn: float
) -> PeriodicTerm:
    """IDIVF 已除入；负 pn 按 cos 偶函数等价改写为正 pn。"""
    if pn < 0:
        pn = -pn
        phase_deg = (-phase_deg) % 360.0
    return PeriodicTerm(types, barrier_kcal * KCAL_KJ, phase_deg % 360.0, pn)


def find_dihedrals(
    lib: ForceFieldLibrary, t1: str, t2: str, t3: str, t4: str
) -> list[PeriodicTerm]:
    """扭转匹配：势垒/IDIVF 后转 kJ；同类型多条全部展开（多项 Fourier）。

    实证：hc-c3-c3-oh（0.16/−3 与 0.25/1 两条）→ kd 0.66944/1.046 kJ、
    phase 0、pn 3/1；X -c3-c3-X（IDIVF 9）→ 1.400/9×4.184 = 0.65084。
    """
    hits = _find_wildcard(
        lib.dihedrals, (t1, t2, t3, t4), lambda e: e.types, lambda e: e.n_wildcards
    )
    return [
        _expand_periodic(e.types, e.barrier_kcal / e.idivf, e.phase_deg, e.pn)
        for e in hits
    ]


def find_improper(
    lib: ForceFieldLibrary, t1: str, t2: str, t3: str, t4: str
) -> PeriodicTerm | None:
    """improper 无 IDIVF；缺失时调用方用 GUESS_IMPROPER_* 常量兜底。"""
    types = (t1, t2, t3, t4)
    for direction in (types, tuple(reversed(types))):
        for entry in lib.impropers:
            if _matches(entry.types, direction):
                return _expand_periodic(
                    entry.types, entry.barrier_kcal, entry.phase_deg, entry.pn
                )
    return None


def find_lj(entries: list[LJEntry], name: str) -> LJEntry | None:
    """按类型名精确查找 LJ 条目。"""
    for entry in entries:
        if entry.name.strip() == name.strip():
            return entry
    return None
