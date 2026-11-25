"""多起点扫描诊断报告（设计 §10.4 / 2026-11-20 计划 §3.1，MVP-7）。

只读已有产物，不调参考势、不重扫、不重拟合：

    --report          多起点 derive-torsions 的 torsion_report.json（行内带 "multistart" 诊断）
    --scan-cache      （可选）同一次运行的扫描缓存目录：读 x0 腿 ``{j}-{k}.npz`` 与各起点腿
                      ``{j}-{k}.{label}.npz``，补逐点网格角、x0 能量与逐腿滞后
    --baseline-report （可选）同分子单起点（现行协议）的 torsion_report.json：对照 fit_rmse

报告内容（计划 §3.1 原文）：现有 hysteresis + 逐点分支数 + 分支间最大能量差 + 包络相对
现行单起点协议的能量变化。**只列数，不判定**：不设阈值、不做自动兜底（2026-10-17 裁定
"先看分布"；设计 §10.4、§26）。分支只是诊断，拟合目标始终是包络。

输出（--outdir）：``report.md``、``per_point.csv``（逐键逐格点）、``summary.json``。

用法：
    python scripts/torsion_multistart_report.py --report RUN/torsion_report.json \
        --scan-cache RUN/scans [--baseline-report BASE/torsion_report.json] --outdir OUT
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def _load_leg(path: Path, n_points: int, quartet: list[int] | None) -> dict | None:
    """读一条扫描腿的 npz（不核缓存键：本脚本只读同一次运行的缓存，核点数与四元组）。"""
    if not path.is_file():
        return None
    with np.load(path) as d:
        energy = np.asarray(d["energy"], dtype=np.float64)
        if len(energy) != n_points:
            raise ValueError(f"{path}: 格点数 {len(energy)} ≠ 报告中的 {n_points}（不是同一次运行？）")
        if quartet is not None and [int(a) for a in d["quartet"]] != [int(a) for a in quartet]:
            raise ValueError(f"{path}: 四元组 {list(d['quartet'])} ≠ 报告中的 {quartet}")
        return {
            "phi_target_deg": np.degrees(np.asarray(d["phi_target_rad"], dtype=np.float64)),
            "energy": energy,
            "hysteresis": np.asarray(d["hysteresis"], dtype=np.float64),
            "converged": np.asarray(d["converged"], dtype=bool),
        }


def _baseline_rows(path: Path | None) -> dict[tuple[int, int], dict]:
    if path is None:
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {tuple(r["bond"]): r for r in data.get("rotors", [])}


def build(report_path: Path, scan_cache: Path | None,
          baseline_path: Path | None) -> tuple[list[dict], list[dict]]:
    """返回 (逐键汇总, 逐键逐点行)。"""
    data = json.loads(report_path.read_text(encoding="utf-8"))
    if "multistart" not in data:
        raise ValueError(f"{report_path}: 没有 multistart 块——不是多起点运行的报告")
    baseline = _baseline_rows(baseline_path)
    summaries: list[dict] = []
    points: list[dict] = []
    for row in data.get("rotors", []):
        diag = row.get("multistart")
        if diag is None:
            continue
        bond = tuple(int(a) for a in row["bond"])
        j, k = bond
        n = int(diag["n_points"])
        labels = list(diag["start_labels"])
        quartet = row.get("scan_quartet")
        legs: dict[str, dict | None] = {}
        if scan_cache is not None:
            for label in labels:
                name = f"{j}-{k}.npz" if label == "x0" else f"{j}-{k}.{label}.npz"
                legs[label] = _load_leg(scan_cache / name, n, quartet)
        x0 = legs.get("x0")
        delta = np.asarray(diag["envelope_delta_kJ_mol"], dtype=np.float64)
        base = baseline.get(bond)
        summary = {
            "bond": list(bond),
            "status": row.get("status"),
            "reason": row.get("reason"),
            "rotor_class": row.get("rotor_class"),
            "start_labels": labels,
            "n_points": n,
            "threshold_nm": diag["threshold_nm"],
            "branch_count_max": diag["branch_count_max"],
            "n_points_multibranch": diag["n_points_multibranch"],
            "branch_gap_max_kJ_mol": diag["branch_gap_max"],
            "envelope_gain_max_kJ_mol": diag["envelope_gain_max"],
            "n_points_improved": diag["n_points_improved"],
            "unconverged_points": diag["unconverged_points"],
            "envelope_max_hysteresis_kJ_mol": row.get("max_hysteresis"),
            "leg_max_hysteresis_kJ_mol": {
                label: (float(leg["hysteresis"].max()) if leg is not None else None)
                for label, leg in legs.items()},
            "missing_legs": [label for label, leg in legs.items() if leg is None],
            "fit_rmse_envelope": row.get("fit_rmse"),
            "fit_rmse_baseline": base.get("fit_rmse") if base else None,
            "baseline_status": base.get("status") if base else None,
            "baseline_max_hysteresis_kJ_mol": base.get("max_hysteresis") if base else None,
        }
        summaries.append(summary)
        x0_rel = env_rel = None
        if x0 is not None:
            x0_e = x0["energy"]
            env_e = x0_e + delta
            x0_rel = x0_e - x0_e.min()
            env_rel = env_e - env_e.min()
        for s in range(n):
            points.append({
                "bond": f"{j}-{k}",
                "point": s,
                "phi_target_deg": float(x0["phi_target_deg"][s]) if x0 is not None else "",
                "branch_count": int(diag["branch_counts"][s]),
                "branch_gap_kJ_mol": float(diag["branch_gap_kJ_mol"][s]),
                "envelope_delta_kJ_mol": float(delta[s]),
                "x0_rel_kJ_mol": float(x0_rel[s]) if x0_rel is not None else "",
                "envelope_rel_kJ_mol": float(env_rel[s]) if env_rel is not None else "",
                "x0_hysteresis_kJ_mol": float(x0["hysteresis"][s]) if x0 is not None else "",
                "envelope_converged": s not in diag["unconverged_points"],
            })
    if not summaries:
        raise ValueError(f"{report_path}: 没有任何键带 multistart 诊断")
    return summaries, points


def _fmt(x, nd=3) -> str:
    if x is None or x == "":
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}"
    return str(x)


def render(summaries: list[dict], points: list[dict], *, report_path: Path,
           scan_cache: Path | None, baseline_path: Path | None) -> str:
    lines = ["# 多起点扫描诊断报告（MVP-7，设计 §10.4）", ""]
    lines.append(f"- 多起点报告：`{report_path}`")
    lines.append(f"- 扫描缓存：`{scan_cache}`" if scan_cache else "- 扫描缓存：未给（无逐点角度/x0 能量/逐腿滞后）")
    lines.append(f"- 单起点基线：`{baseline_path}`" if baseline_path else "- 单起点基线：未给（无 fit_rmse 对照）")
    lines.append("- 口径：包络 = 逐点全部 起点×方向 中收敛候选的最低能量；envelope_delta = 包络 − x0 腿"
                 "（≤ 0）；分支 = 同一格点候选按对齐重原子 RMSD 贪心聚类（阈值见表）。")
    lines.append("- 只列数、不判定：不设阈值、不自动兜底；分支只进诊断，拟合目标是包络。")
    lines.append("")
    lines.append("| 键 | 状态 | 起点 | 阈值 nm | 最大分支数 | 多分支点 | 分支间最大能差 | 包络最大增益 "
                 "| 改善点数 | 未收敛点 | x0 腿最大滞后 | 包络最大滞后 | fit_rmse 包络 | fit_rmse 单起点 |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in summaries:
        j, k = s["bond"]
        status = s["status"] + (f"（{s['reason']}）" if s.get("reason") else "")
        lines.append(
            f"| {j}-{k} | {status} | {len(s['start_labels'])} | {s['threshold_nm']} "
            f"| {s['branch_count_max']} | {s['n_points_multibranch']}/{s['n_points']} "
            f"| {_fmt(s['branch_gap_max_kJ_mol'])} | {_fmt(s['envelope_gain_max_kJ_mol'])} "
            f"| {s['n_points_improved']} | {len(s['unconverged_points'])} "
            f"| {_fmt(s['leg_max_hysteresis_kJ_mol'].get('x0'))} "
            f"| {_fmt(s['envelope_max_hysteresis_kJ_mol'])} "
            f"| {_fmt(s['fit_rmse_envelope'])} | {_fmt(s['fit_rmse_baseline'])} |")
    lines.append("")
    lines.append("能量单位 kJ/mol。\"包络最大滞后\" 是胜出候选所在起点自身正反两遍的能差（逐点），"
                 "与 x0 腿滞后口径相同。")
    for s in summaries:
        j, k = s["bond"]
        lines.append("")
        lines.append(f"## 键 {j}-{k}")
        lines.append("")
        lines.append(f"- 起点：{', '.join(s['start_labels'])}")
        legs = s["leg_max_hysteresis_kJ_mol"]
        if legs:
            lines.append("- 逐腿最大滞后（kJ/mol）：" + "，".join(
                f"{label} {_fmt(v)}" for label, v in legs.items()))
        if s["missing_legs"]:
            lines.append(f"- 缓存里缺的腿：{', '.join(s['missing_legs'])}")
        if s["baseline_status"] is not None:
            lines.append(f"- 单起点基线：status={s['baseline_status']}，最大滞后 "
                         f"{_fmt(s['baseline_max_hysteresis_kJ_mol'])}")
        lines.append("")
        lines.append("| 点 | φ 目标 ° | 分支数 | 分支间能差 | 包络−x0 | x0 相对能 | 包络相对能 | x0 滞后 | 收敛 |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for p in points:
            if p["bond"] != f"{j}-{k}":
                continue
            lines.append(
                f"| {p['point']} | {_fmt(p['phi_target_deg'], 1)} | {p['branch_count']} "
                f"| {_fmt(p['branch_gap_kJ_mol'])} | {_fmt(p['envelope_delta_kJ_mol'])} "
                f"| {_fmt(p['x0_rel_kJ_mol'])} | {_fmt(p['envelope_rel_kJ_mol'])} "
                f"| {_fmt(p['x0_hysteresis_kJ_mol'])} | {'是' if p['envelope_converged'] else '否'} |")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--scan-cache", type=Path, default=None)
    parser.add_argument("--baseline-report", type=Path, default=None)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.scan_cache is not None and not args.scan_cache.is_dir():
        parser.error(f"--scan-cache 不是目录：{args.scan_cache}")

    summaries, points = build(args.report, args.scan_cache, args.baseline_report)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "summary.json").write_text(
        json.dumps({"argv": sys.argv, "bonds": summaries}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")
    with open(args.outdir / "per_point.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(points[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)
    (args.outdir / "report.md").write_text(
        render(summaries, points, report_path=args.report, scan_cache=args.scan_cache,
               baseline_path=args.baseline_report), encoding="utf-8")
    print(f"[multistart-report] {len(summaries)} 键 → {args.outdir / 'report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
