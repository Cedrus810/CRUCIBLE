"""覆盖集诊断的纯逻辑层（2026-11-20 计划 §2）：记录 schema、verdict 判定、汇总报告。

计划 §2.2 的七项诊断字段映射到 :class:`MoleculeRecord` 的分级 stage 记录；驱动脚本
（``scripts/coverage_diagnose.py``）负责跑外部配方并填表，本模块只做不依赖外部工具的
部分：frcmod 缺项解析、verdict 判定、逐类别汇总。

verdict 规则（计划 §2.2 第 7 项，施工前锁定）：
- ``blocked``  —— entry/typing 失败（SDF 合同或 antechamber/parmchk2/tleap 配方失败），
  或导出等价检查失败/出错（三后端不一致 = 产物不可信）；
  例外（用户决定 2026-11-24，同决策 0013）：超差行全部来自混合精度 GROMACS 腿时，导出
  阶段记 ok + ``gmx_mixed_precision`` 缺项，不判 blocked；
- ``ok_with_gaps`` —— 管线走通但有缺项：frcmod 有 ATTN（GAFF2 库缺参数）、参考势
  适用性检查不通过（元素不在模型里、开壳层）、或后续拟合步产生了兜底项；
- ``ok`` —— 其余。

不进 pipeline、不加门（计划 §2.4）：本模块的输出只用于诊断报告。
"""

from __future__ import annotations

import csv
import io
from dataclasses import asdict, dataclass, field
from pathlib import Path

VERDICT_OK = "ok"
VERDICT_GAPS = "ok_with_gaps"
VERDICT_BLOCKED = "blocked"

# 计划 §2.1 的八个化学类别（键名 ASCII，便于 CSV/JSON；说明文字进报告）
CATEGORIES: dict[str, str] = {
    "alcohol_ether": "醇/醚",
    "carbonyl_ester_amide": "羰基/酯/酰胺",
    "amine": "胺",
    "heteroaromatic": "杂芳环",
    "aromatic_conjugated_fused": "芳香共轭/稠环/桥环",
    "small_ring": "小环",
    "heteroatom_halogen": "含硫/磷/卤/氰",
    "ion_zwitterion": "离子/两性离子",
}

_FRCMOD_SECTIONS = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON", "LJEDIT")


def parse_frcmod_gaps(frcmod: Path | str) -> dict:
    """解析 parmchk2 的 frcmod：GAFF2 库缺参数以 ``ATTN, need revision`` 标注。

    返回 ``{"n_attn": int, "sections": {"DIHE": 3, ...}, "attn_lines": [...]}``；
    文件不存在或不是 frcmod 结构时报 :class:`ValueError`（entry 阶段失败时调用方
    不会走到这一步，所以这里按"有文件但内容异常"处理）。
    """
    text = Path(frcmod).read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"{frcmod}: 空文件，不是 frcmod 结构")
    sections: dict[str, int] = {}
    attn_lines: list[str] = []
    current = None
    for raw in text.splitlines():
        line = raw.strip()
        if line in _FRCMOD_SECTIONS:
            current = line
            sections.setdefault(current, 0)
            continue
        if "ATTN" in line.upper():
            attn_lines.append(line)
            if current is not None:
                sections[current] = sections.get(current, 0) + 1
    return {
        "n_attn": len(attn_lines),
        "sections": {k: v for k, v in sections.items() if v},
        "attn_lines": attn_lines,
    }


def verdict_from_stages(stages: dict[str, dict]) -> tuple[str, list[str]]:
    """按锁定的规则从各 stage 状态判 verdict，返回 ``(verdict, blockers)``。

    stage 状态约定：``{"status": "ok"|"fail"|"skip", ...}``（``fail`` 必带 ``error``）；
    gap 标记放在 ``gaps`` 键（列表，原因分类与计划 §2.2 的字段一一对应）。
    """
    blockers: list[str] = []
    gaps: list[str] = []
    entry = stages.get("entry", {})
    if entry.get("status") == "fail":
        blockers.append(f"entry: {entry.get('error', 'unspecified')}")
    export = stages.get("export", {})
    if export.get("status") == "fail":
        blockers.append(f"export_not_equivalent: {export.get('error', 'unspecified')}")
    if blockers:
        return VERDICT_BLOCKED, blockers
    for name in ("export", "rotors", "frcmod", "refcheck", "fallback", "mace_smoke", "xtb_sp",
                 "scan_probe"):
        stage = stages.get(name)
        if stage is None or stage.get("status") == "skip":
            continue
        if stage.get("status") == "fail":
            gaps.append(f"{name}: {stage.get('error', 'unspecified')}")
        gaps.extend(stage.get("gaps", []))
    if gaps:
        return VERDICT_GAPS, gaps
    return VERDICT_OK, []


@dataclass
class MoleculeRecord:
    """单分子的覆盖诊断记录（计划 §2.2 七项字段的结构化形式）。"""

    name: str
    category: str
    net_charge: int
    n_atoms: int
    sdf_sha256: str
    stages: dict = field(default_factory=dict)   # stage 名 -> {"status", ...}
    verdict: str = VERDICT_BLOCKED               # 由 verdict_from_stages 判定后写入
    blockers: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.category not in CATEGORIES:
            raise ValueError(
                f"{self.name}: 未知类别 {self.category!r}；可选：{', '.join(CATEGORIES)}")

    def finalize(self) -> None:
        self.verdict, self.blockers = verdict_from_stages(self.stages)

    def to_dict(self) -> dict:
        return asdict(self)


def load_table(path: Path | str) -> list[dict]:
    """读覆盖集清单 JSON：``[{"name", "category", "sdf", "net_charge"}]``。

    sdf 路径相对清单文件解析（绝对路径原样）；逐条检查 name/category/sdf 齐全、
    name 在清单内唯一、category 在 :data:`CATEGORIES` 内。
    """
    import json

    path = Path(path)
    items = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(items, list) or not items:
        raise ValueError(f"{path}: 清单必须是含至少一条记录的 JSON 数组")
    seen: set[str] = set()
    out: list[dict] = []
    for i, item in enumerate(items):
        missing = [k for k in ("name", "category", "sdf") if k not in item]
        if missing:
            raise ValueError(f"{path}[{i}]: 缺字段 {missing}")
        name = str(item["name"])
        if name in seen:
            raise ValueError(f"{path}[{i}]: 分子名重复 {name!r}")
        seen.add(name)
        sdf = Path(item["sdf"])
        if not sdf.is_absolute():
            sdf = (path.parent / sdf).resolve()
        charge = int(item.get("net_charge", 0))
        row = {"name": name, "category": item["category"], "sdf": sdf,
               "net_charge": charge}
        MoleculeRecord(name=name, category=row["category"], net_charge=charge,
                       n_atoms=0, sdf_sha256="")  # 借 __post_init__ 校验类别
        out.append(row)
    return out


def write_summary_csv(records: list[MoleculeRecord], path: Path) -> None:
    """逐分子一行：name/category/net_charge/n_atoms/verdict/各 stage 状态/gaps。"""
    stage_names = ("entry", "rotors", "frcmod", "export", "refcheck",
                   "fallback", "mace_smoke", "xtb_sp", "scan_probe")
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["name", "category", "net_charge", "n_atoms", "verdict",
                     "blockers", *stage_names])
    for rec in records:
        writer.writerow(
            [rec.name, rec.category, rec.net_charge, rec.n_atoms, rec.verdict,
             "; ".join(rec.blockers),
             *[rec.stages.get(s, {}).get("status", "-") for s in stage_names]])
    path.write_text(buf.getvalue(), encoding="utf-8")


def summarize(records: list[MoleculeRecord]) -> str:
    """汇总报告（markdown）：逐类别统计表 + 逐分子状态表 + blocker 分类。"""
    lines = ["# 覆盖集诊断汇总", ""]
    lines.append("| 类别 | 分子数 | ok | ok_with_gaps | blocked |")
    lines.append("|---|---|---|---|---|")
    for cat, label in CATEGORIES.items():
        rows = [r for r in records if r.category == cat]
        if not rows:
            continue
        counts = {v: sum(1 for r in rows if r.verdict == v)
                  for v in (VERDICT_OK, VERDICT_GAPS, VERDICT_BLOCKED)}
        lines.append(f"| {label} ({cat}) | {len(rows)} | {counts[VERDICT_OK]} "
                     f"| {counts[VERDICT_GAPS]} | {counts[VERDICT_BLOCKED]} |")
    unknown = [r for r in records if r.category not in CATEGORIES]
    if unknown:  # 理论不可达（__post_init__ 已挡）；防御式保留
        lines.append(f"| （未知类别） | {len(unknown)} | | | |")
    lines.append("")
    lines.append("| 分子 | 类别 | q | verdict | 说明 |")
    lines.append("|---|---|---|---|---|")
    for rec in sorted(records, key=lambda r: (CATEGORIES.get(r.category, ""), r.name)):
        note = "; ".join(rec.blockers) if rec.blockers else ""
        lines.append(f"| {rec.name} | {CATEGORIES.get(rec.category, rec.category)} "
                     f"| {rec.net_charge:+d} | {rec.verdict} | {note} |")
    lines.append("")
    lines.append("verdict 规则与字段定义见 `docs/superpowers/plans/"
                 "2026-11-20-coverage-validation-batch1-plan.md` §2.2；"
                 "本报告只用于诊断，不进 pipeline、不加门（计划 §2.4）。")
    return "\n".join(lines) + "\n"
