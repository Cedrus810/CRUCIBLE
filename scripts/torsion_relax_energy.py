"""能量度量 vs 受约束 MM 松弛（docs/theory/innovation_screen_20261029.md §8 用户边界第二条）。

只读 G2 扫描缓存，只算 MM。对键 r 的每个扫描几何 x，E_fixed = GAFF2 去掉 r 上的 proper，Φ = r 上的
全部四元组（取线性无关子集；被去掉的那个由恒等式随之固定，脚本里逐帧核对）：

    relax   ：E_fixed(x) − min{E_fixed(y) : Φ(y) = Φ(x)}。ASE FixInternals 固定 Φ，LBFGS。与扭转参数无关。
    second  ：½ g_Zᵀ (Zᵀ H_L Z)⁻¹ g_Z。H_L = H − Σ μ_i ∇²φ_i，μ 为 ∇E = Jᵀμ 的最小二乘乘子；
              Z 为 [J; 6 个刚体方向] 零空间的正交基。仅当 Zᵀ H_L Z 正定时有定义（局部稳定）。
    nocurv  ：同 second，但 H_L 换成普通 MM Hessian H（不含约束曲率项）。

H 由 OpenMM 力的中心差分得到，∇²φ 由解析梯度的中心差分得到。这里的松弛不同于 OpenFF 只约束扫描二面角
且保留拟合扭转项的重新松弛。

    python scripts/torsion_relax_energy.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_force_bound_20261101/relax --molecules methyl_benzoate
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from torsion_force_bound import _forces, _load_fixture, dihedral_gradient, quartets_on_bond  # noqa: E402

from ffcompiler.ir.energy import dihedral_angle  # noqa: E402
from ffcompiler.units import ANGSTROM_NM, EV_KJ_MOL  # noqa: E402

H_HESS_NM = 1e-5      # MM Hessian 差分步长
H_DIHED_NM = 1e-6     # 二面角 Hessian 差分步长
PERTURB_NM = 0.01     # 负曲率最低模的扰动幅度（单位向量 × 0.01 nm，单原子位移 ≤ 0.01 nm）


def independent_quartets(x, quartets):
    keep, rows = [], []
    for q in quartets:
        trial = rows + [dihedral_gradient(x, q).ravel()]
        if np.linalg.matrix_rank(np.stack(trial)) == len(trial):
            keep.append(q)
            rows = trial
    return keep


def rigid_modes(x):
    n = len(x)
    c = x - x.mean(axis=0)
    modes = []
    for a in range(3):
        t = np.zeros((n, 3))
        t[:, a] = 1.0
        modes.append(t.ravel())
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        modes.append(np.cross(e, c).ravel())
    return np.stack(modes)


def mm_hessian(ctx, x):
    n3 = x.size
    h = np.empty((n3, n3))
    for m in range(n3):
        xp, xm = x.ravel().copy(), x.ravel().copy()
        xp[m] += H_HESS_NM
        xm[m] -= H_HESS_NM
        h[m] = -(_forces(ctx, xp.reshape(x.shape)).ravel() - _forces(ctx, xm.reshape(x.shape)).ravel()) / (2 * H_HESS_NM)
    asym = float(np.max(np.abs(h - h.T)))
    return 0.5 * (h + h.T), asym


def dihedral_hessian(x, q):
    n3 = x.size
    h = np.empty((n3, n3))
    for m in range(n3):
        xp, xm = x.ravel().copy(), x.ravel().copy()
        xp[m] += H_DIHED_NM
        xm[m] -= H_DIHED_NM
        h[m] = (dihedral_gradient(xp.reshape(x.shape), q).ravel()
                - dihedral_gradient(xm.reshape(x.shape), q).ravel()) / (2 * H_DIHED_NM)
    return 0.5 * (h + h.T)


def second_order(g, hess, jac, x):
    c = np.vstack([jac, rigid_modes(x)])
    _, s, vt = np.linalg.svd(c)
    rank = int(np.sum(s > s.max() * max(c.shape) * np.finfo(float).eps))
    z = vt[rank:].T
    gz = z.T @ g
    hz = z.T @ hess @ z
    w = np.linalg.eigvalsh(hz)
    val = 0.5 * float(gz @ np.linalg.solve(hz, gz)) if w.min() > 0 else float("nan")
    _, vecs = np.linalg.eigh(hz)
    return val, float(w.min()), z.shape[1], z @ vecs[:, 0]


def out_of_plane_fraction(x, v, heavy):
    """v 在重原子最佳拟合平面法向上的平方范数占比。"""
    c = x[heavy] - x[heavy].mean(axis=0)
    normal = np.linalg.svd(c)[2][-1]
    vv = v.reshape(-1, 3)
    return float(np.sum((vv @ normal) ** 2) / np.sum(vv ** 2))


def _relax(ir_fixed, numbers, x, quartets, all_quartets, fmax_ev_a, max_steps, start=None):
    """从 start（默认 x）出发、Φ 固定在 x 的值上松弛；能量差总以 E_fixed(x) 为起点。"""
    from ase import Atoms
    from ase.constraints import FixInternals
    from ase.optimize import LBFGS

    from ffcompiler.torsion.mm import OpenMMCalculator

    atoms = Atoms(numbers=numbers, positions=x / ANGSTROM_NM)
    atoms.calc = OpenMMCalculator(ir_fixed, validate_ir=False)
    e0 = atoms.get_potential_energy() * EV_KJ_MOL
    targets = [[atoms.get_dihedral(*q), list(q)] for q in quartets]
    if start is not None:
        atoms.positions = start / ANGSTROM_NM
    atoms.set_constraint(FixInternals(dihedrals_deg=targets))
    conv = bool(LBFGS(atoms, logfile=None).run(fmax=fmax_ev_a, steps=max_steps))
    e1 = atoms.get_potential_energy() * EV_KJ_MOL
    y = atoms.positions * ANGSTROM_NM
    atoms.set_constraint()
    dphi = max(abs(np.angle(np.exp(1j * (dihedral_angle(*(y[t] for t in q)) - dihedral_angle(*(x[t] for t in q))))))
               for q in all_quartets)
    # Kabsch 对齐后的最大原子位移
    a, b = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    disp = float(np.max(np.linalg.norm(b @ rot - a, axis=1)))
    return e0 - e1, conv, float(dphi), disp


def analyse_molecule(name, g2_run: Path, outdir: Path, fmax_ev_a: float, max_steps: int) -> list[dict]:
    from ffcompiler.exporters.openmm import to_openmm_system
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider
    from ffcompiler.torsion.mm import strip_bond_torsions
    from ffcompiler.validation.engines import reference_context

    base, pos, mol = _load_fixture(name)
    provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                  scan_cache=g2_run / "scans" / name, model_path=default_model_path(),
                                  cache_only=True)
    _, report = provider.parameterize(base, pos, mol)
    numbers = [a.atomic_number for a in base.atoms]
    heavy = [i for i, z in enumerate(numbers) if z > 1]
    summary = []
    for row in report.rotors:
        if row["status"] != "mlp":
            continue
        bond = tuple(row["bond"])
        label = f"{bond[0]}-{bond[1]}"
        scan = report.scans[bond]
        ir_fixed = strip_bond_torsions(base, bond)
        ctx = reference_context(to_openmm_system(ir_fixed, validate_ir=False))
        all_q = quartets_on_bond(base, bond)
        frames = []
        for s, x in enumerate(scan.positions):
            qs = independent_quartets(x, all_q)
            jac = np.stack([dihedral_gradient(x, q).ravel() for q in qs])
            g = -_forces(ctx, x).ravel()
            mu = np.linalg.lstsq(jac.T, g, rcond=None)[0]
            hess, asym = mm_hessian(ctx, x)
            h_l = hess - sum(m * dihedral_hessian(x, q) for m, q in zip(mu, qs))
            e2, wmin2, nz, vmin = second_order(g, h_l, jac, x)
            e2n, wmin2n, _, _ = second_order(g, hess, jac, x)
            er, conv, dphi, disp = _relax(ir_fixed, numbers, x, qs, all_q, fmax_ev_a, max_steps)
            # 负曲率时，沿最低模 ±PERTURB_NM 扰动再松弛，检查是否停在受约束鞍点
            er_p, conv_p, dphi_p, disp_p = er, conv, dphi, disp
            if wmin2 < 0:
                for sign in (1.0, -1.0):
                    start = x + sign * PERTURB_NM * vmin.reshape(x.shape)
                    cand = _relax(ir_fixed, numbers, x, qs, all_q, fmax_ev_a, max_steps, start=start)
                    if cand[0] > er_p:
                        er_p, conv_p, dphi_p, disp_p = cand
            frames.append({"frame": s, "phi_deg": float(np.degrees(scan.phi_rad[s])),
                           "n_quartets": len(all_q), "n_independent": len(qs), "tangent_dim": nz,
                           "relax": er, "relax_converged": conv, "relax_max_dphi_rad": dphi,
                           "relax_max_disp_nm": disp,
                           "relax_pert": er_p, "relax_pert_converged": conv_p,
                           "relax_pert_max_dphi_rad": dphi_p, "relax_pert_max_disp_nm": disp_p,
                           "minmode_out_of_plane": out_of_plane_fraction(x, vmin, heavy),
                           "second": e2, "second_min_eig": wmin2,
                           "nocurv": e2n, "nocurv_min_eig": wmin2n,
                           "check_hess_asym": asym})
            print(f"  {label} s={s:2d} relax={er:8.4f} pert={er_p:8.4f} second={e2:8.4f} nocurv={e2n:8.4f} "
                  f"wmin={wmin2:9.1f} conv={conv}/{conv_p}", flush=True)
        fpath = outdir / f"{name}_{label}.csv"
        fpath.parent.mkdir(parents=True, exist_ok=True)
        with open(fpath, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(frames[0]))
            w.writeheader()
            w.writerows(frames)
        a = lambda k: np.array([f[k] for f in frames], dtype=float)
        rel = a("relax")
        summary.append({"molecule": name, "bond": label,
                        "relax_min": float(rel.min()), "relax_median": float(np.median(rel)),
                        "relax_max": float(rel.max()), "relax_range": float(rel.max() - rel.min()),
                        "relax_all_converged": bool(all(f["relax_converged"] for f in frames)),
                        "relax_max_dphi_rad": float(a("relax_max_dphi_rad").max()),
                        "relax_max_disp_nm": float(a("relax_max_disp_nm").max()),
                        "relax_pert_min": float(a("relax_pert").min()),
                        "relax_pert_median": float(np.median(a("relax_pert"))),
                        "relax_pert_max": float(a("relax_pert").max()),
                        "relax_pert_range": float(a("relax_pert").max() - a("relax_pert").min()),
                        "relax_pert_gain_max": float((a("relax_pert") - rel).max()),
                        "relax_pert_all_converged": bool(all(f["relax_pert_converged"] for f in frames)),
                        "relax_pert_max_dphi_rad": float(a("relax_pert_max_dphi_rad").max()),
                        "relax_pert_max_disp_nm": float(a("relax_pert_max_disp_nm").max()),
                        "minmode_out_of_plane_median": float(np.median(a("minmode_out_of_plane"))),
                        "second_nan_frames": int(np.isnan(a("second")).sum()),
                        "second_over_relax_median": float(np.nanmedian(a("second") / rel)),
                        "second_max_abs_err": float(np.nanmax(np.abs(a("second") - rel))),
                        "nocurv_nan_frames": int(np.isnan(a("nocurv")).sum()),
                        "nocurv_over_relax_median": float(np.nanmedian(a("nocurv") / rel)),
                        "nocurv_max_abs_err": float(np.nanmax(np.abs(a("nocurv") - rel))),
                        "check_hess_asym_max": float(a("check_hess_asym").max())})
    return summary


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
    with open(args.outdir / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    manifest = {"task": "energy metric vs constrained MM relaxation (innovation_screen_20261029 §8)",
                "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
                "wall_s": time.time() - t0, "fmax_eV_per_A": args.fmax, "max_steps": args.max_steps,
                "h_hess_nm": H_HESS_NM, "h_dihed_nm": H_DIHED_NM, "units": "energy kJ/mol"}
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {len(rows)} rows, {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
