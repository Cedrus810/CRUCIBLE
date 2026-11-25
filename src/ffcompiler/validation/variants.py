"""对照矩阵的 IR 组合（2026-11-20 计划 §C1"分别替换"）。

GAFF2 基线（入口骨架 IR）与完整管线 IR 之间的替换矩阵：
``substitute_terms(base, final, parts=...)`` 从 ``final`` 取指定部分、其余保留
``base``，得到"只换 bonded / 只换 torsion / 只换 charge / 全换"等变体。
每个变体都是可导出、可求值的完整 IR（经 ``validate``），逐类别指标在报告层算。

parts 语义（可组合）：
- ``bonded``  —— bonds + angles（含 global-fit 的 bond/angle 改动）；
- ``torsion`` —— proper_torsions + impropers（MLP 拟合项与兜底项一起换）；
- ``charge``  —— 逐原子 partial_charge（等价类原子的等价性由各 provider 保证）；
- ``vdw``     —— 逐原子 sigma/epsilon（GAFF2 LJ 目前两边一致，为对称性保留）。

不做：键合扩展项（耦合/非谐）——IR 还没有这些项（计划 §4 待基准暴露系统误差后立项）。
"""

from __future__ import annotations

from dataclasses import replace
from typing import Iterable

from ..ir.schema import ForceFieldIR
from ..ir.validation import validate

PARTS = frozenset({"bonded", "torsion", "charge", "vdw"})


def substitute_terms(
    base: ForceFieldIR, final: ForceFieldIR, parts: Iterable[str]
) -> ForceFieldIR:
    """从 ``final`` 取 ``parts``、其余保留 ``base``。parts 为空 = base 原样（再加 provenance 标记）。"""
    parts = set(parts)
    unknown = parts - PARTS
    if unknown:
        raise ValueError(f"未知 parts：{sorted(unknown)}；可选 {sorted(PARTS)}")
    if len(base.atoms) != len(final.atoms):
        raise ValueError(
            f"原子数不一致：base {len(base.atoms)}，final {len(final.atoms)}（必须同一分子）")

    atoms = base.atoms
    if "charge" in parts or "vdw" in parts:
        atoms = tuple(
            replace(
                a,
                partial_charge=f.partial_charge if "charge" in parts else a.partial_charge,
                sigma=f.sigma if "vdw" in parts else a.sigma,
                epsilon=f.epsilon if "vdw" in parts else a.epsilon,
            )
            for a, f in zip(base.atoms, final.atoms)
        )
    if "charge" in parts and base.total_charge != final.total_charge:
        # replace_terms 保留 base.total_charge；两边总电荷不同说明不是同一分子的电荷
        raise ValueError(
            f"总电荷不一致：base {base.total_charge}，final {final.total_charge}")
    ir = base.replace_terms(
        atoms=atoms,
        bonds=final.bonds if "bonded" in parts else base.bonds,
        angles=final.angles if "bonded" in parts else base.angles,
        proper_torsions=final.proper_torsions if "torsion" in parts else base.proper_torsions,
        impropers=final.impropers if "torsion" in parts else base.impropers,
        provenance_update={"variant": {
            "parts": sorted(parts),
            "from": {"base": base.provenance.get("name", base.name),
                     "final": final.provenance.get("name", final.name)},
        }},
    )
    validate(ir)
    return ir


def comparison_matrix(base: ForceFieldIR, final: ForceFieldIR) -> dict[str, ForceFieldIR]:
    """计划 §C1 的标准对照矩阵：GAFF2 基线 + 单替换 + 完整组合。

    键即报告口径：gaff2 / bonded / torsion / charge / vdw / full。
    """
    out = {"gaff2": base}
    for part in sorted(PARTS):
        out[part] = substitute_terms(base, final, {part})
    out["full"] = substitute_terms(base, final, PARTS)
    return out
