"""拟合依赖的快照与过期检查（2026-10-15 计划 C3/C4，2026-10-17 施工）。

两类拟合各有必须记住的输入：

- **global_fit 的 H_fixed**：拟合时扣掉的固定部分由 IR 的电荷、二面角
  （proper/improper）、LJ、1-4 缩放、exclusions/pairs14 决定（§6.4.2）。
  拟合把这些写进 ``provenance["bonded"]["h_fixed_dependency"]``；
- **MLP 二面角的 MM 扣除**：逐键扫描时减掉的 MM 能量用当时的
  bond/angle/improper/电荷/LJ（``mm_energy_without_bond``），写进
  ``provenance["torsion"]["mm_dependency"]``。**其他键的 proper 二面角也在
  扣除能量里，但不进快照**：单遍拟合本来就用输入 IR 的 GAFF2 项扣除、最后
  再把各键拟合项组合（10-15 计划 B2 要诊断的近似），记进来会让每个刚拟合
  完的多可转键 IR 立即判为过期。这一近似由 B2/B3 处理，不在过期检查里。
  可选顺序更新另记不参与更新的固定 proper；互相更新的目标项由 fit_protocol
  和逐轮报告描述，不声称有限轮后达到了联合最优。

之后任何一步改了这些输入而没有重拟合，对应依赖即过期：
:func:`staleness` 把**记录的**快照与**当前** IR 的快照逐字段比对，报告过期
项。快照存逐项数值（JSON 可往返），**不用 hash 判断**（计划 C3：hash 无法
说明哪里变了，且对 key 顺序敏感）；字段相等用严格 ``==``（float 经
``repr`` 最短往返后逐位一致）。
"""

from __future__ import annotations

from typing import Any

from ffcompiler.ir import ForceFieldIR

# term 行只记影响物理的数值字段；meta（拟合诊断等）不属于依赖
_TORSION_FIELDS = ("i", "j", "k", "l", "periodicity", "phase", "force_constant", "source")


def _torsion_row(t: Any) -> list:
    return [getattr(t, name) for name in _TORSION_FIELDS]


def _nonbonded_rows(ir: ForceFieldIR) -> list[list[float]]:
    return [[float(a.sigma), float(a.epsilon)] for a in ir.atoms]


def h_fixed_snapshot(ir: ForceFieldIR) -> dict[str, Any]:
    """进入 H_fixed 的全部输入（global_fit 拟合时记录、过期检查时重算）。"""
    return {
        "charges": [float(a.partial_charge) for a in ir.atoms],
        "charge_source": _plain(ir.provenance.get("charge")),
        "torsions": [_torsion_row(t) for t in ir.proper_torsions],
        "impropers": [_torsion_row(t) for t in ir.impropers],
        "lj": _nonbonded_rows(ir),
        "combination_rule": ir.combination_rule,
        "coulomb_14_scale": ir.coulomb_14_scale,
        "vdw_14_scale": ir.vdw_14_scale,
        "exclusions": [list(p) for p in ir.exclusions],
        "pairs14": [list(p) for p in ir.pairs14],
    }


def mm_subtraction_snapshot(ir: ForceFieldIR) -> dict[str, Any]:
    """MLP 二面角 MM 扣除的输入（bond/angle/improper/电荷/LJ 与缩放）。

    其他键的 proper 二面角不记（单遍拟合的已知近似，见模块说明）。"""
    return {
        "charges": [float(a.partial_charge) for a in ir.atoms],
        "bonds": [[b.i, b.j, float(b.r0), float(b.k)] for b in ir.bonds],
        "angles": [[a.i, a.j, a.k_atom, float(a.theta0), float(a.force_constant)]
                   for a in ir.angles],
        "impropers": [_torsion_row(t) for t in ir.impropers],
        "lj": _nonbonded_rows(ir),
        "combination_rule": ir.combination_rule,
        "coulomb_14_scale": ir.coulomb_14_scale,
        "vdw_14_scale": ir.vdw_14_scale,
        "exclusions": [list(p) for p in ir.exclusions],
        "pairs14": [list(p) for p in ir.pairs14],
    }


def fixed_proper_snapshot(ir: ForceFieldIR, fit_bonds) -> list[list]:
    """顺序更新时保留不动的 proper（含原兜底键），不含互相更新的目标键。"""
    selected = {frozenset(b) for b in fit_bonds}
    return [_torsion_row(t) for t in ir.proper_torsions
            if frozenset((t.j, t.k)) not in selected]


def _plain(value: Any) -> Any:
    """provenance 子树转 JSON 可往返的纯数据（tuple→list 递归）。"""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


_LABELS = {
    "charges": "电荷",
    "charge_source": "电荷来源",
    "torsions": "proper 二面角",
    "impropers": "improper",
    "bonds": "bond",
    "angles": "angle",
    "lj": "LJ",
    "combination_rule": "混合规则",
    "coulomb_14_scale": "1-4 静电缩放",
    "vdw_14_scale": "1-4 LJ 缩放",
    "exclusions": "exclusions",
    "pairs14": "pairs14",
}


def _diff_fields(recorded: dict, current: dict) -> list[str]:
    changed = []
    for key in sorted(set(recorded) | set(current)):
        if recorded.get(key) != current.get(key):
            changed.append(_LABELS.get(key, key))
    return changed


def staleness(ir: ForceFieldIR) -> list[str]:
    """过期项清单（空 = 没有过期）。

    - bonded 为 global_fit 且记录了 h_fixed_dependency：与当前
      :func:`h_fixed_snapshot` 比对；
    - torsion 为 mlp 且记录了 mm_dependency：与当前
      :func:`mm_subtraction_snapshot` 比对；
    - 没记录依赖（旧 IR）→ 报"依赖未记录"，同样算过期（不允许静默沿用）。
    """
    items: list[str] = []

    bonded = ir.provenance.get("bonded") or {}
    if bonded.get("method") == "global_fit":
        dep = bonded.get("h_fixed_dependency")
        if dep is None:
            items.append("global_fit: H_fixed 依赖未记录（旧版本 IR），需重拟合")
        else:
            # dep = {"inputs": 快照, "hessian": 身份, "fd_step_nm": …}；只比快照
            for field in _diff_fields(dep["inputs"], h_fixed_snapshot(ir)):
                items.append(f"global_fit: H_fixed 的 {field} 在拟合后已改变")

    torsion = ir.provenance.get("torsion") or {}
    if torsion.get("provider") == "mlp" and any(
        t.source == "mlp" for t in ir.proper_torsions
    ):
        dep = torsion.get("mm_dependency")
        if dep is None:
            items.append("mlp torsion: MM 扣除依赖未记录（旧版本 IR），需重拟合")
        else:
            for field in _diff_fields(dep, mm_subtraction_snapshot(ir)):
                items.append(f"mlp torsion: MM 扣除的 {field} 在拟合后已改变")
        protocol = torsion.get("fit_protocol", {})
        if (not isinstance(protocol, dict)
                or (protocol and (type(protocol.get("refinement_rounds")) is not int
                    or not 0 <= protocol["refinement_rounds"] <= 3
                    or protocol.get("name") != (
                        "sequential" if protocol["refinement_rounds"] else "single_pass")))):
            items.append("mlp torsion: fit_protocol 无效，需重新参数化")
            return items
        if protocol.get("name") == "sequential":
            fixed = torsion.get("fixed_proper_dependency")
            expected_order = sorted([list(r["bond"]) for r in torsion.get("rotors", [])
                                     if r.get("status") == "mlp"])
            if protocol.get("order") != expected_order:
                items.append("mlp torsion: fit_protocol 顺序与原 mlp 键不一致")
            elif fixed is None:
                items.append("mlp torsion: 顺序更新的固定 proper 依赖未记录，需重拟合")
            elif fixed != fixed_proper_snapshot(ir, protocol["order"]):
                items.append("mlp torsion: 顺序更新的固定 proper 在拟合后已改变")

    return items
