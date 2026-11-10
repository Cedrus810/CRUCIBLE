"""受约束松弛检查：拟合前预测重新松弛后的剖面误差（2026-11-01 计划阶段 1）。

只读 G2 扫描缓存，只算 MM，不加载 MACE、不调用 xtb、不重扫。对键 r 的每个扫描几何
x_s、参考能量 E_ref(s)，IR ∈ {only, all}（同 B2 定义：only 只把 r 换成拟合项，all 为 G2 最终 IR）：

    e_fix(s)      = E_IR(x_s) − E_ref(s)                         固定几何误差（G2/B2 现有口径）
    e_rel(s)      = min{E_IR(y) : φ_q(y) = φ_q(x_s)} − E_ref(s)   只约束扫描四元组 q，逐帧从 x_s
                                                                   出发松弛（OpenFF 式重新松弛）
    relax(s)      = E_bg(x_s) − min{E_bg(y) : Φ_r(y) = Φ_r(x_s)}  E_bg = 该 IR 去掉 r 上的 proper
                                                                   （only → GAFF2 背景，all → 最终 IR
                                                                   背景），与 θ_r 无关
    e_rel_pred(s) = e_fix(s) − [relax(s) − mean(relax)]           预测量 X

    P0 = RMSE_bs(e_fix)；Y = RMSE_bs(e_rel) − P0；X = RMSE_bs(e_rel_pred) − P0；P1 = 扫描最大滞后。
    RMSE_bs = 最佳常数平移 RMSE，kJ/mol。误差曲线（e_fix、e_rel、e_rel_pred）本身已是 E_IR − E_ref，
    其最佳常数平移口径就是去均值后的 RMS，不再减 E_ref。e_rel 的松弛用完整 IR（保留 r 上的拟合
    扭转项）；X 的松弛固定 Φ_r 上全部四元组（线性无关子集，被去掉的由恒等式随之固定，逐帧核对）。

    OMP_NUM_THREADS=1 python scripts/torsion_relaxed_check.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_relaxed_check_20261111 --molecules methyl_benzoate

每键一个 frames/<分子>_<键>.csv，重跑同一命令时跳过已有 CSV 的键（按文件存在判断）。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from ase.constraints import FixConstraint  # noqa: E402

from torsion_force_bound import _load_fixture, dihedral_gradient, quartets_on_bond  # noqa: E402
from torsion_relax_energy import independent_quartets  # noqa: E402

from ffcompiler.ir.energy import dihedral_angle, evaluate  # noqa: E402
from ffcompiler.units import ANGSTROM_NM, EV_KJ_MOL  # noqa: E402


class MultiDihedral(FixConstraint):
    """同时固定多个二面角（弧度目标，ffcompiler 口径）：位置做同时 Gauss-Newton 投影，
    力投影掉 range(Jᵀ)。

    ASE FixInternals 对每个约束做一次质量加权单步校正、外层 50 轮逐约束迭代、容差
    epsilon=1e-7°；dmso C–S 键的 bg 松弛要同时固定 4 个同轴二面角，在该方案下不收敛
    （2026-11-11 实测：默认参数、maxstep 0.02、epsilon 1e-5 均失败；epsilon 1e-4 能过
    但约束漂移 1.44e-6 rad，超出验收门 1e-6 rad）。本类改为对全部约束同时解
    J δ = Δφ（J 为解析二面角梯度堆叠，行满秩，Newton 二次收敛，容差 1e-12 rad）；
    可行集合、fmax 判据与 LBFGS 优化器都不变。J 在 Å 坐标上计算（梯度公式与长度单位无关），
    力投影为 F − Jᵀ(JJᵀ)⁻¹JF，与 FixInternals 的欧氏投影同口径。
    """

    def __init__(self, quartets, targets_rad, tol=1e-12, max_iter=30):
        self.quartets = [tuple(int(a) for a in q) for q in quartets]
        self.targets = [float(t) for t in targets_rad]
        self.tol = float(tol)
        self.max_iter = int(max_iter)

    def _phi(self, x, q):
        return dihedral_angle(*(x[a] for a in q))

    def _jac(self, x):
        return np.stack([dihedral_gradient(x, q).ravel() for q in self.quartets])

    def adjust_positions(self, atoms, newpos):
        dphi = None
        for _ in range(self.max_iter):
            dphi = np.array([float(_wrap(self._phi(newpos, q) - t))
                             for q, t in zip(self.quartets, self.targets)])
            if float(np.max(np.abs(dphi))) < self.tol:
                return
            jac = self._jac(newpos)
            delta, *_ = np.linalg.lstsq(jac, -dphi, rcond=None)
            newpos += delta.reshape(newpos.shape)
        raise ValueError("MultiDihedral.adjust_positions did not converge "
                         f"(max |dphi| = {float(np.max(np.abs(dphi))):.2e} rad)")

    def adjust_forces(self, atoms, forces):
        jac = self._jac(atoms.positions)
        flat = forces.ravel()
        flat -= jac.T @ np.linalg.solve(jac @ jac.T, jac @ flat)
        forces[:] = flat.reshape(forces.shape)

    def copy(self):
        return MultiDihedral(self.quartets, self.targets, self.tol, self.max_iter)

    def index_shuffle(self, atoms, ind):
        raise IndexError("MultiDihedral not part of slice")

    def get_indices(self):
        return sorted({a for q in self.quartets for a in q})

    def __repr__(self):
        return f"MultiDihedral(n={len(self.quartets)}, tol={self.tol:g})"


def _on_bond(t, bond) -> bool:
    return {t.j, t.k} == set(bond)


def _wrap(a) -> float:
    return float(np.angle(np.exp(1j * a)))


def _bestshift_rmse(a, b) -> float:
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean((d - d.mean()) ** 2)))


def _relax_q(ir, numbers, x, q, others, fmax_ev_a, max_steps):
    """只约束扫描四元组 q（φ 取 x 的值，弧度），从 x 出发松弛完整 IR；e0 在加约束前求值。

    others：Φ_r 上除 q 外的四元组，其漂移是数据（机制证据），不是错误。
    返回 (e0, e1, conv, dphi_q, dphi_others, disp)，能量 kJ/mol，disp 为 Kabsch 对齐后最大位移 nm。
    """
    from ase import Atoms
    from ase.optimize import LBFGS

    from ffcompiler.torsion.mm import OpenMMCalculator

    atoms = Atoms(numbers=numbers, positions=x / ANGSTROM_NM)
    atoms.calc = OpenMMCalculator(ir, validate_ir=False)
    e0 = atoms.get_potential_energy() * EV_KJ_MOL
    atoms.set_constraint(MultiDihedral([q], [dihedral_angle(*(x[t] for t in q))]))
    conv = bool(LBFGS(atoms, logfile=None).run(fmax=fmax_ev_a, steps=max_steps))
    e1 = atoms.get_potential_energy() * EV_KJ_MOL
    y = atoms.positions * ANGSTROM_NM
    atoms.set_constraint()
    dq = abs(_wrap(dihedral_angle(*(y[t] for t in q)) - dihedral_angle(*(x[t] for t in q))))
    doth = max((_wrap(dihedral_angle(*(y[t] for t in o)) - dihedral_angle(*(x[t] for t in o)))
                for o in others), default=0.0)
    a, b = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    disp = float(np.max(np.linalg.norm(b @ rot - a, axis=1)))
    return e0, e1, conv, float(dq), float(doth), disp


def _relax_bg(ir_fixed, numbers, x, quartets, all_quartets, fmax_ev_a, max_steps):
    """Φ_r（线性无关子集）固定在 x 的值上松弛 E_bg；能量差以 E_fixed(x) 为起点。

    与 torsion_relax_energy._relax 同一定义（该脚本 11-01 用它算出 relax 报告），但约束
    改用 MultiDihedral——ASE FixInternals 在 dmso 的 4-同轴约束上不收敛（见类 docstring）。
    返回 (drop, conv, dphi, disp)。
    """
    from ase import Atoms
    from ase.optimize import LBFGS

    from ffcompiler.torsion.mm import OpenMMCalculator

    atoms = Atoms(numbers=numbers, positions=x / ANGSTROM_NM)
    atoms.calc = OpenMMCalculator(ir_fixed, validate_ir=False)
    e0 = atoms.get_potential_energy() * EV_KJ_MOL
    atoms.set_constraint(MultiDihedral(quartets,
                                       [dihedral_angle(*(x[t] for t in q)) for q in quartets]))
    conv = bool(LBFGS(atoms, logfile=None).run(fmax=fmax_ev_a, steps=max_steps))
    e1 = atoms.get_potential_energy() * EV_KJ_MOL
    y = atoms.positions * ANGSTROM_NM
    atoms.set_constraint()
    dphi = max(abs(_wrap(dihedral_angle(*(y[t] for t in q)) - dihedral_angle(*(x[t] for t in q))))
               for q in all_quartets)
    a, b = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    disp = float(np.max(np.linalg.norm(b @ rot - a, axis=1)))
    return e0 - e1, conv, float(dphi), disp


def _write_frames(fpath: Path, frames: list[dict]) -> None:
    fpath.parent.mkdir(parents=True, exist_ok=True)
    with open(fpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(frames[0]))
        w.writeheader()
        w.writerows(frames)


def _read_frames(fpath: Path) -> list[dict]:
    with open(fpath, newline="") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        rec = {}
        for k, v in r.items():
            if k in ("molecule", "bond"):
                rec[k] = v
            else:
                rec[k] = float(v)
        out.append(rec)
    return out


def _summarise(name: str, label: str, frames: list[dict]) -> list[dict]:
    """frames 含 only/all 两套逐帧量；每变体出一行 summary。"""
    def demean_rms(c):
        # 误差曲线的最佳常数平移 RMSE：曲线本身已是 E_IR − E_ref，平移对齐到零 = 去均值
        return float(np.sqrt(np.mean((c - c.mean()) ** 2)))

    a = lambda k: np.array([f[k] for f in frames], dtype=float)
    p1 = float(a("hysteresis").max())
    rows = []
    for v in ("only", "all"):
        efix, erel = a(f"{v}_e_fix"), a(f"{v}_e_rel")
        bg = a(f"{v}_bg_relax")
        epred = efix - (bg - bg.mean())
        p0 = demean_rms(efix)
        y = demean_rms(erel) - p0
        x = demean_rms(epred) - p0
        rms_diff_pred = float(np.sqrt(np.mean((erel - epred) ** 2)))
        rms_diff_fix = float(np.sqrt(np.mean((erel - efix) ** 2)))
        rows.append({
            "molecule": name, "bond": label, "variant": v, "n_frames": len(frames),
            "scan_all_converged": float(min(f["converged"] for f in frames)),
            "p1_max_hysteresis": p1,
            "p0_rmse_bs": p0, "rmse_bs_e_rel": demean_rms(erel),
            "Y": y, "rmse_bs_e_rel_pred": demean_rms(epred), "X": x,
            "rms_e_rel_minus_pred": rms_diff_pred, "rms_e_rel_minus_fix": rms_diff_fix,
            "rms_ratio_pred_over_fix": (rms_diff_pred / rms_diff_fix) if rms_diff_fix > 0 else float("nan"),
            "bg_relax_median": float(np.median(bg)), "bg_relax_min": float(bg.min()),
            "bg_relax_max": float(bg.max()), "bg_relax_range": float(bg.max() - bg.min()),
            "relax_drop_min": float(a(f"{v}_relax_drop").min()),
            "relax_drop_median": float(np.median(a(f"{v}_relax_drop"))),
            "all_conv_q": float(min(f[f"{v}_relax_conv"] for f in frames)),
            "all_conv_bg": float(min(f[f"{v}_bg_conv"] for f in frames)),
            "check_dphi_q_max": float(a(f"{v}_dphi_q_rad").max()),
            "check_dphi_bg_max": float(a(f"{v}_bg_dphi_rad").max()),
            "check_efix_numpy_vs_omm_max": float(np.max(np.abs(a(f"{v}_e_fix_numpy") - efix))),
            "check_relax_drop_ge0_min": float(a(f"{v}_relax_drop").min()),
            "check_bg_relax_ge0_min": float(bg.min()),
            "dphi_others_max": float(a(f"{v}_dphi_others_rad").max()),
            "disp_q_max_nm": float(a(f"{v}_max_disp_nm").max()),
            "disp_bg_max_nm": float(a(f"{v}_bg_max_disp_nm").max()),
        })
    return rows


def analyse_molecule(name: str, g2_run: Path, outdir: Path, fmax_ev_a: float, max_steps: int) -> list[dict]:
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider
    from ffcompiler.torsion.mm import strip_bond_torsions

    base, pos, mol = _load_fixture(name)
    provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                  scan_cache=g2_run / "scans" / name, model_path=default_model_path(),
                                  cache_only=True)
    redo, report = provider.parameterize(base, pos, mol)
    numbers = [a.atomic_number for a in base.atoms]
    summary = []
    for row in report.rotors:
        bond = tuple(row["bond"])
        label = f"{bond[0]}-{bond[1]}"
        fpath = outdir / "frames" / f"{name}_{label}.csv"
        if fpath.exists():
            frames = _read_frames(fpath)
            print(f"  {label}: {fpath} 已存在，跳过", flush=True)
            summary.extend(_summarise(name, label, frames))
            continue
        if row["status"] != "mlp":
            continue
        scan = report.scans[bond]
        fitted_r = [t for t in redo.proper_torsions if _on_bond(t, bond)]
        irs = {"only": base.replace_terms(proper_torsions=tuple(
                   [t for t in base.proper_torsions if not _on_bond(t, bond)] + fitted_r)),
               "all": redo}
        bgs = {"only": strip_bond_torsions(base, bond), "all": strip_bond_torsions(redo, bond)}
        all_q = quartets_on_bond(base, bond)
        q = tuple(scan.quartet)
        assert q in all_q or tuple(reversed(q)) in all_q
        others = [t for t in all_q if t != q and t != tuple(reversed(q))]
        e_ref = np.asarray(scan.energy, dtype=float)
        efix_np = {v: np.array([evaluate(irs[v], xi)["total"] for xi in scan.positions])
                   for v in ("only", "all")}
        frames = []
        for s, x in enumerate(scan.positions):
            qs = independent_quartets(x, all_q)
            rec = {"frame": s, "phi_deg": float(np.degrees(scan.phi_rad[s])),
                   "converged": float(bool(scan.converged[s])), "hysteresis": float(scan.hysteresis[s]),
                   "E_ref": float(e_ref[s])}
            for v in ("only", "all"):
                e0, e1, conv, dq, doth, disp = _relax_q(irs[v], numbers, x, q, others,
                                                        fmax_ev_a, max_steps)
                drop_bg, conv_bg, dphi_bg, disp_bg = _relax_bg(bgs[v], numbers, x, qs, all_q,
                                                                fmax_ev_a, max_steps)
                rec.update({
                    f"{v}_e_fix": e0 - e_ref[s], f"{v}_e_rel": e1 - e_ref[s],
                    f"{v}_relax_drop": e0 - e1, f"{v}_relax_conv": float(conv),
                    f"{v}_dphi_q_rad": dq, f"{v}_dphi_others_rad": doth,
                    f"{v}_max_disp_nm": disp,
                    f"{v}_bg_relax": drop_bg, f"{v}_bg_conv": float(conv_bg),
                    f"{v}_bg_dphi_rad": dphi_bg, f"{v}_bg_max_disp_nm": disp_bg,
                    f"{v}_e_fix_numpy": float(efix_np[v][s]) - e_ref[s],
                })
            frames.append(rec)
            print(f"  {label} s={s:2d} only: fix={rec['only_e_fix']:8.4f} rel={rec['only_e_rel']:8.4f} "
                  f"bg={rec['only_bg_relax']:7.4f} | all: fix={rec['all_e_fix']:8.4f} "
                  f"rel={rec['all_e_rel']:8.4f} bg={rec['all_bg_relax']:7.4f}", flush=True)
        _write_frames(fpath, frames)
        summary.extend(_summarise(name, label, frames))
    return summary


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g2-run", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--molecules", nargs="+", required=True)
    ap.add_argument("--fmax", type=float, default=0.001, help="松弛收敛 fmax，eV/Å")
    ap.add_argument("--max-steps", type=int, default=5000)
    args = ap.parse_args(argv)
    t0 = time.time()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in args.molecules:
        rows.extend(analyse_molecule(name, args.g2_run, args.outdir, args.fmax, args.max_steps))
        print(f"[{time.time() - t0:7.1f}s] {name} done", flush=True)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(args.outdir / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    sources = [Path(__file__), ROOT / "scripts" / "torsion_force_bound.py",
               ROOT / "scripts" / "torsion_relax_energy.py"] + sorted((ROOT / "src" / "ffcompiler").rglob("*.py"))
    manifest = {
        "task": "constrained relaxation check: pre-fit prediction of re-relaxed profile error "
                "(2026-11-01 relaxed-check plan, phase 1)",
        "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
        "started_unix": t0, "wall_s": time.time() - t0,
        "g2_run": str(args.g2_run), "molecules": args.molecules,
        "fmax_eV_per_A": args.fmax, "max_steps": args.max_steps,
        "sigma_kJ_mol": 0.41,
        "inputs_sha256": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): _sha256(p)
                          for p in [args.g2_run / "metrics.csv",
                                    *sorted(p for p in (args.g2_run / "ir").glob("*.json")
                                            if p.stem in args.molecules),
                                    *sorted(p for p in (args.g2_run / "scans").rglob("*.npz")
                                            if p.relative_to(args.g2_run / "scans").parts[0] in args.molecules)]},
        "source_sha256": {str(p.relative_to(ROOT)): _sha256(p) for p in sources},
        "note": "sha256 只作身份记录，不用于复用判断；续跑按 frames CSV 文件存在判断",
        "units": "energy kJ/mol",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {len(rows)} rows, {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
