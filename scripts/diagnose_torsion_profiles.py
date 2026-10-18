"""B1/B2：MLP 二面角剖面诊断（2026-10-15 计划 §4；2026-10-17 施工）。

只读一个 G2 运行目录（IR + scans/），不加载 MACE、不调用 xtb、不重扫。

B1：用缓存按 G2 的单遍流程重做参数化，核对与存档 IR 逐项相同；按 G2 原口径
    （各减自身最小值的 RMSE）重算 fit_rmse / profile_rmse_fit / profile_rmse_gaff2，
    与 metrics.csv（6 位有效数字）比较。
B2：对每根键 r，在它自己的扫描构型上比较三套 IR 的逐项能量：
    base = 原 GAFF2；only = 只把 r 换成拟合项；all = G2 最终 IR（全部 mlp 键都换）。
    分项：r 上的 proper、其他 proper、improper、bond、angle、LJ、静电、1-4 LJ、1-4 静电。
    检查三个恒等式：
      (1) 分项之和 = NumPy 整体求值 = OpenMM 整体求值；
      (2) E_all − E_only = Σ_{s≠r}[T_s^fit − T_s^base]（其他项逐帧相同）；
      (3) only 的最佳常数平移误差 RMS = Fourier 拟合残差 RMS（单遍拟合的定义）。

    python scripts/diagnose_torsion_profiles.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_diagnosis_20261017 [--molecules methyl_benzoate ...]

输出：manifest.json、summary.csv、frames/<分子>_<键>.csv、plots/<分子>_<键>.png、report.md。
不用 RunContext（它探测 git，NFS 工作区禁跑 git）；源码身份记 sha256 清单（只作身份，不作复用判断）。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ffcompiler.ir import io as ir_io  # noqa: E402
from ffcompiler.ir.energy import ENERGY_KEYS, evaluate  # noqa: E402

FIXTURE_DIR = ROOT / "tests" / "fixtures" / "gaff2"
COMPONENTS = ("target", "other_proper", "improper", "bond", "angle", "lj", "coulomb", "lj14", "coulomb14")
VARIANTS = ("base", "only", "all")


def _load_fixture(name: str):
    from rdkit import Chem

    from ffcompiler.importers import from_amber

    mdir = FIXTURE_DIR / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


def _on_bond(t, bond) -> bool:
    return {t.j, t.k} == set(bond)


def _components(ir, bond, x) -> dict[str, float]:
    """一帧的分项能量；r 上的 proper 与其他 proper 分开（各用一个只含该组 proper 的 IR 求值）。"""
    full = evaluate(ir, x)
    target_ir = ir.replace_terms(proper_torsions=tuple(t for t in ir.proper_torsions if _on_bond(t, bond)))
    target = evaluate(target_ir, x)["proper"]
    out = {"target": target, "other_proper": full["proper"] - target}
    for key in ENERGY_KEYS:
        if key != "proper":
            out[key] = full[key]
    out["total_numpy"] = full["total"]
    return out


def _proper_energy(ir, bond, x) -> float:
    sub = ir.replace_terms(proper_torsions=tuple(t for t in ir.proper_torsions if _on_bond(t, bond)))
    return evaluate(sub, x)["proper"]


def _openmm_profile(ir, positions_stack) -> np.ndarray:
    from openmm import unit

    from ffcompiler.exporters.openmm import to_openmm_system
    from ffcompiler.validation.engines import reference_context

    ctx = reference_context(to_openmm_system(ir))
    out = []
    for x in positions_stack:
        ctx.setPositions(x)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    return np.asarray(out)


def _minshift_rmse(a, b) -> float:  # G2 原口径（scripts/gate_mvp2.py::_shifted_rmse）
    return float(np.sqrt(np.mean(((a - a.min()) - (b - b.min())) ** 2)))


def _bestshift_rmse(a, b) -> float:
    d = a - b
    return float(np.sqrt(np.mean((d - d.mean()) ** 2)))


def _barrier(e) -> float:
    return float(np.max(e) - np.min(e))


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_metrics(g2_run: Path) -> dict[tuple[str, str], dict]:
    with open(g2_run / "metrics.csv", newline="") as fh:
        return {(r["molecule"], r["bond"]): r for r in csv.DictReader(fh)}


def diagnose_molecule(name: str, g2_run: Path, outdir: Path, metrics: dict) -> list[dict]:
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider

    base, pos, mol = _load_fixture(name)
    stored = ir_io.from_json(g2_run / "ir" / f"{name}.json")
    provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                  scan_cache=g2_run / "scans" / name, model_path=default_model_path(),
                                  cache_only=True)
    redo, report = provider.parameterize(base, pos, mol)
    # 与存档 IR 比较：项的原子/周期/来源须逐项相同，数值报最大差（跨机器 BLAS 末位不同，不要求逐位）
    same_layout = [(t.i, t.j, t.k, t.l, t.periodicity, t.source) for t in redo.proper_torsions] == \
        [(t.i, t.j, t.k, t.l, t.periodicity, t.source) for t in stored.proper_torsions]
    max_dk = max((abs(a.force_constant - b.force_constant)
                  for a, b in zip(redo.proper_torsions, stored.proper_torsions)), default=0.0)
    max_dphase = max((abs(np.angle(np.exp(1j * (a.phase - b.phase))))
                      for a, b in zip(redo.proper_torsions, stored.proper_torsions)), default=0.0)

    rows = []
    for row in report.rotors:
        bond = tuple(row["bond"])
        label = f"{bond[0]}-{bond[1]}"
        rec = {"molecule": name, "bond": label, "status": row["status"], "reason": row["reason"] or "",
               "g2_ir_same_terms": same_layout, "g2_ir_max_dk": max_dk, "g2_ir_max_dphase_rad": max_dphase}
        scan = report.scans.get(bond)
        if row["status"] != "mlp" or scan is None:
            rows.append(rec)
            continue
        fitted_r = [t for t in redo.proper_torsions if _on_bond(t, bond)]
        only = base.replace_terms(proper_torsions=tuple(
            [t for t in base.proper_torsions if not _on_bond(t, bond)] + fitted_r))
        irs = {"base": base, "only": only, "all": redo}
        x = scan.positions
        e_ref = np.asarray(scan.energy, dtype=np.float64)
        comps = {v: [_components(irs[v], bond, xi) for xi in x] for v in VARIANTS}
        tot = {v: np.array([sum(c[k] for k in COMPONENTS) for c in comps[v]]) for v in VARIANTS}
        tot_np = {v: np.array([c["total_numpy"] for c in comps[v]]) for v in VARIANTS}
        tot_mm = {v: _openmm_profile(irs[v], x) for v in VARIANTS}

        # (1) 分项和 = NumPy 整体 = OpenMM
        sum_err = max(float(np.max(np.abs(tot[v] - tot_np[v]))) for v in VARIANTS)
        omm_err = max(float(np.max(np.abs(tot_np[v] - tot_mm[v]))) for v in VARIANTS)
        # (2) E_all − E_only 只来自其他键的 proper
        other_delta = np.array([ca["other_proper"] - co["other_proper"]
                                for ca, co in zip(comps["all"], comps["only"])])
        non_other = max(float(np.max(np.abs(np.array([ca[k] - co[k] for ca, co in
                                                      zip(comps["all"], comps["only"])]))))
                        for k in COMPONENTS if k != "other_proper")
        ident_err = float(np.max(np.abs((tot["all"] - tot["only"]) - other_delta)))
        # (3) only 的最佳平移误差 = 拟合残差
        from ffcompiler.torsion.fourier import fit_fourier
        from ffcompiler.torsion.mm import mm_energy_without_bond

        y = e_ref - mm_energy_without_bond(base, bond, x)
        fit = fit_fourier(scan.phi_rad, y)

        # 其他各键的贡献：沿本键扫描，GAFF2 项、拟合项与二者之差各自的变化幅度
        others = {}
        for other in report.rotors:
            sb = tuple(other["bond"])
            if sb == bond or other["status"] != "mlp":
                continue
            t_b = np.array([_proper_energy(base, sb, xi) for xi in x])
            t_f = np.array([_proper_energy(redo, sb, xi) for xi in x])
            others[f"{sb[0]}-{sb[1]}"] = {"range_gaff2": _barrier(t_b), "range_fit": _barrier(t_f),
                                          "range_delta": _barrier(t_f - t_b),
                                          "bestshift_rms_delta": _bestshift_rmse(t_f - t_b, 0.0 * t_b)}
        m = metrics.get((name, label), {})
        rec.update({
            "n_frames": int(len(e_ref)), "scan_quartet": "-".join(map(str, scan.quartet)),
            "all_converged": bool(np.all(scan.converged)),
            "max_hysteresis": float(np.max(scan.hysteresis)), "max_bond_change_pm": scan.max_bond_change_pm,
            "fit_rmse_report": row["fit_rmse"], "fit_rmse_recomputed": fit.rmse,
            "minshift_base": _minshift_rmse(tot_mm["base"], e_ref),
            "minshift_only": _minshift_rmse(tot_mm["only"], e_ref),
            "minshift_all": _minshift_rmse(tot_mm["all"], e_ref),
            "bestshift_base": _bestshift_rmse(tot_mm["base"], e_ref),
            "bestshift_only": _bestshift_rmse(tot_mm["only"], e_ref),
            "bestshift_all": _bestshift_rmse(tot_mm["all"], e_ref),
            "barrier_ref": _barrier(e_ref), "barrier_base": _barrier(tot_mm["base"]),
            "barrier_only": _barrier(tot_mm["only"]), "barrier_all": _barrier(tot_mm["all"]),
            "phi_min_ref_deg": float(np.degrees(scan.phi_rad[np.argmin(e_ref)])),
            "phi_min_all_deg": float(np.degrees(scan.phi_rad[np.argmin(tot_mm["all"])])),
            "other_delta_range": _barrier(other_delta),
            "other_delta_bestshift_rms": _bestshift_rmse(other_delta, np.zeros_like(other_delta)),
            "other_bonds": json.dumps(others, sort_keys=True),
            "check_sum_vs_numpy": sum_err, "check_numpy_vs_openmm": omm_err,
            "check_identity_all_minus_only": ident_err, "check_nonproper_all_vs_only": non_other,
            "g2_fit_rmse": m.get("fit_rmse", ""), "g2_profile_rmse_fit": m.get("profile_rmse_fit", ""),
            "g2_profile_rmse_gaff2": m.get("profile_rmse_gaff2", ""),
        })
        rows.append(rec)

        frames = outdir / "frames" / f"{name}_{label}.csv"
        frames.parent.mkdir(parents=True, exist_ok=True)
        with open(frames, "w", newline="") as fh:
            w = csv.writer(fh)
            head = ["frame", "phi_target_deg", "phi_deg", "converged", "hysteresis", "E_ref",
                    "fit_resid"]
            head += [f"{v}_{k}" for v in VARIANTS for k in COMPONENTS]
            head += [f"{v}_total_openmm" for v in VARIANTS]
            w.writerow(head)
            from ffcompiler.torsion.fourier import fourier_energy

            resid = y - fourier_energy(fit, scan.phi_rad)
            for i in range(len(e_ref)):
                line = [i, f"{np.degrees(scan.phi_target_rad[i]):.6f}", f"{np.degrees(scan.phi_rad[i]):.6f}",
                        int(scan.converged[i]), repr(float(scan.hysteresis[i])), repr(float(e_ref[i])),
                        repr(float(resid[i]))]
                line += [repr(float(comps[v][i][k])) for v in VARIANTS for k in COMPONENTS]
                line += [repr(float(tot_mm[v][i])) for v in VARIANTS]
                w.writerow(line)
        _plot(outdir / "plots" / f"{name}_{label}.png", name, label, scan.phi_rad, e_ref, tot_mm)
    return rows


def _plot(path: Path, name: str, label: str, phi, e_ref, tot) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    order = np.argsort(phi)
    deg = np.degrees(np.asarray(phi))[order]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(deg, (e_ref - e_ref.min())[order], "k-o", ms=3, label="MACE (ref)")
    for v, style in (("base", "C0--"), ("only", "C1-"), ("all", "C3-")):
        e = tot[v]
        ax.plot(deg, (e - e.min())[order], style, label=f"MM {v}")
    ax.set_xlabel("scan dihedral (deg)")
    ax.set_ylabel("E − min (kJ/mol)")
    ax.set_title(f"{name} {label}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g2-run", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--molecules", nargs="+", default=None)
    args = ap.parse_args(argv)
    t0 = time.time()
    g2 = args.g2_run
    names = args.molecules or sorted(p.stem for p in (g2 / "ir").glob("*.json")
                                     if (g2 / "scans" / p.stem).is_dir())
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics = _read_metrics(g2)
    rows = []
    for name in names:
        rows.extend(diagnose_molecule(name, g2, args.outdir, metrics))
        print(f"[{time.time() - t0:7.1f}s] {name} done", flush=True)

    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(args.outdir / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    sources = [Path(__file__)] + sorted((ROOT / "src" / "ffcompiler").rglob("*.py"))
    manifest = {
        "task": "B1/B2 torsion profile diagnosis (2026-10-15 plan §4)",
        "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
        "started_unix": t0, "wall_s": time.time() - t0,
        "g2_run": str(g2), "molecules": names,
        "inputs_sha256": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): _sha256(p)
                          for p in [g2 / "metrics.csv", *sorted((g2 / "ir").glob("*.json")),
                                    *sorted((g2 / "scans").rglob("*.npz"))]},
        "source_sha256": {str(p.relative_to(ROOT)): _sha256(p) for p in sources},
        "note": "sha256 只作身份记录，不用于复用判断",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {len(rows)} rows, {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
