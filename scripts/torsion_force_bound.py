"""候选 A 冒烟：扭转参数无法修复的力下界（docs/theory/innovation_screen_20261029.md §3.1）。

只读一个 G2 运行目录的扫描缓存，只算 MM，不加载 MACE、不调用 xtb、不重扫。

量：在缓存的扫描几何 x 上，令 J 为可调二面角集合 Φ 的 Cartesian Jacobian，P⊥ 为 range(Jᵀ)
的欧氏正交补投影。扫描用 ASE FixInternals 只约束扫描四元组 q；它的 adjust_forces 是欧氏投影，
投掉 ∇φ_q 与 6 个刚体方向。参考势平移/转动不变，所以未投影参考力 F_ref = λ∇φ_q + δ，
δ 是优化器看到的受约束力，收敛时每个原子 |δ_i| < fmax，‖δ‖ ≤ √N·fmax。q ∈ Φ 时 P⊥∇φ_q = 0，
于是对 r = F_ref − F_fixed：

    ‖P⊥F_fixed‖ − √N·fmax ≤ ‖P⊥r‖ ≤ ‖P⊥F_fixed‖ + √N·fmax

两种冻结问题：
    bond：Φ = 键 r 上的全部四元组；E_fixed = GAFF2 去掉 r 上的 proper。
    all ：Φ = 全部 mlp 键上的全部四元组；E_fixed = GAFF2 去掉这些键上的 proper。
秩用 numpy 默认数值容差（只反映双精度舍入，不是自定的判据）。

    python scripts/torsion_force_bound.py --g2-run jobs/gate_mvp2_20261013 \
        --outdir jobs/torsion_force_bound_20261101 --molecules methyl_benzoate
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

from ffcompiler.ir.energy import dihedral_angle  # noqa: E402
from ffcompiler.units import ANGSTROM_NM, EV_KJ_MOL  # noqa: E402

FIXTURE_DIR = ROOT / "tests" / "fixtures" / "gaff2"


def _load_fixture(name: str):
    from rdkit import Chem

    from ffcompiler.importers import from_amber

    mdir = FIXTURE_DIR / name
    ir, pos = from_amber(*(mdir / f"{name}.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / f"{name}.sdf"), removeHs=False)
    return ir, pos, mol


def quartets_on_bond(ir, bond) -> list[tuple[int, int, int, int]]:
    """分子图上中心键为 bond 的全部四元组（不论 IR 里有没有对应的项）。"""
    nbr: dict[int, set[int]] = {}
    for b in ir.bonds:
        nbr.setdefault(b.i, set()).add(b.j)
        nbr.setdefault(b.j, set()).add(b.i)
    j, k = bond
    return [(i, j, k, l) for i in sorted(nbr[j] - {k}) for l in sorted(nbr[k] - {j}) if i != l]


def dihedral_gradient(x: np.ndarray, q) -> np.ndarray:
    """∂φ/∂x，(N,3)，nm⁻¹（Bekker/Blondel–Karplus 解析式；符号由有限差分自检确定）。"""
    i, j, k, l = q
    f, g, h = x[i] - x[j], x[j] - x[k], x[l] - x[k]
    a, b = np.cross(f, g), np.cross(h, g)
    gn = np.linalg.norm(g)
    aa, bb = a @ a, b @ b
    out = np.zeros_like(x)
    out[i] = -gn / aa * a
    out[l] = gn / bb * b
    fg, hg = f @ g, h @ g
    out[j] = gn / aa * a + fg / (aa * gn) * a - hg / (bb * gn) * b
    out[k] = -gn / bb * b - fg / (aa * gn) * a + hg / (bb * gn) * b
    return out


def _fd_gradient(x: np.ndarray, q, h: float = 1e-6) -> np.ndarray:
    out = np.zeros_like(x)
    for a in q:
        for c in range(3):
            xp, xm = x.copy(), x.copy()
            xp[a, c] += h
            xm[a, c] -= h
            d = dihedral_angle(*(xp[t] for t in q)) - dihedral_angle(*(xm[t] for t in q))
            out[a, c] = np.angle(np.exp(1j * d)) / (2 * h)
    return out


def perp_projector_apply(jac: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, int, np.ndarray]:
    """jac：(m, 3N) 行为二面角梯度。返回 P⊥v、数值秩、奇异值。"""
    u, s, vt = np.linalg.svd(jac, full_matrices=False)
    rank = int(np.linalg.matrix_rank(jac))
    basis = vt[:rank]
    return v - basis.T @ (basis @ v), rank, s


def _forces(ctx, x) -> np.ndarray:
    from openmm import unit

    ctx.setPositions(x)
    return ctx.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)


def decompose_bond(name: str, bond, scan, ir_fixed, quartets, outdir: Path) -> list[dict]:
    """bond 变体：P⊥F_fixed 按 MM 分项拆开（OpenMM force group），附各项能量。

    share_c = ⟨P⊥F_c, P⊥F⟩ / ‖P⊥F‖²，各项相加为 1。
    """
    from openmm import unit

    from ffcompiler.exporters.openmm import FORCE_GROUPS, to_openmm_system
    from ffcompiler.validation.engines import reference_context

    ctx = reference_context(to_openmm_system(ir_fixed, decompose=True, validate_ir=False))
    label = f"{bond[0]}-{bond[1]}"
    rows = []
    for s, x in enumerate(scan.positions):
        jac = np.stack([dihedral_gradient(x, q).ravel() for q in quartets])
        ctx.setPositions(x)
        total = _forces(ctx, x).ravel()
        perp_tot = perp_projector_apply(jac, total)[0]
        n2 = float(perp_tot @ perp_tot)
        row = {"molecule": name, "bond": label, "frame": s, "phi_deg": float(np.degrees(scan.phi_rad[s])),
               "norm_perp": float(np.sqrt(n2))}
        fsum = np.zeros_like(total)
        for comp, grp in FORCE_GROUPS.items():
            st = ctx.getState(getForces=True, getEnergy=True, groups={grp})
            fc = st.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer).ravel()
            fsum += fc
            pc = perp_projector_apply(jac, fc)[0]
            row[f"{comp}_perp"] = float(np.linalg.norm(pc))
            row[f"{comp}_share"] = float(pc @ perp_tot / n2)
            row[f"{comp}_E"] = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        row["check_groupsum_vs_total"] = float(np.max(np.abs(fsum - total)))
        rows.append(row)
    fpath = outdir / "decompose" / f"{name}_{label}_bond.csv"
    fpath.parent.mkdir(parents=True, exist_ok=True)
    with open(fpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def analyse_molecule(name: str, g2_run: Path, outdir: Path, decompose: bool = False) -> list[dict]:
    from ffcompiler.exporters.openmm import to_openmm_system
    from ffcompiler.torsion.mace import default_model_path
    from ffcompiler.torsion.mlp import MLPTorsionProvider
    from ffcompiler.validation.engines import reference_context

    base, pos, mol = _load_fixture(name)
    provider = MLPTorsionProvider(lambda ir_: None, supported=frozenset(range(1, 119)),
                                  scan_cache=g2_run / "scans" / name, model_path=default_model_path(),
                                  cache_only=True)
    _, report = provider.parameterize(base, pos, mol)
    rotors = [tuple(r["bond"]) for r in report.rotors if r["status"] == "mlp"]

    def strip(bonds):
        keep = tuple(t for t in base.proper_torsions if not any({t.j, t.k} == set(b) for b in bonds))
        return base.replace_terms(proper_torsions=keep)

    ctx_all = reference_context(to_openmm_system(strip(rotors), validate_ir=False))
    q_all = [q for b in rotors for q in quartets_on_bond(base, b)]

    rows = []
    for bond in rotors:
        scan = report.scans[bond]
        label = f"{bond[0]}-{bond[1]}"
        fmax_kj_nm = float(scan.params["fmax_eV_per_A"]) * EV_KJ_MOL / ANGSTROM_NM
        n_atoms = scan.positions.shape[1]
        delta_max = np.sqrt(n_atoms) * fmax_kj_nm
        variants = {
            "bond": (reference_context(to_openmm_system(strip([bond]), validate_ir=False)),
                     quartets_on_bond(base, bond)),
            "all": (ctx_all, q_all),
        }
        if decompose:
            decompose_bond(name, bond, scan, strip([bond]), quartets_on_bond(base, bond), outdir)
            continue
        for vname, (ctx, quartets) in variants.items():
            assert tuple(scan.quartet) in quartets or tuple(reversed(scan.quartet)) in quartets
            frames = []
            for s, x in enumerate(scan.positions):
                jac = np.stack([dihedral_gradient(x, q).ravel() for q in quartets])
                if s == 0:  # 解析梯度对有限差分自检
                    fd_err = max(float(np.max(np.abs(dihedral_gradient(x, q) - _fd_gradient(x, q))))
                                 for q in quartets)
                    jmax = float(np.max(np.abs(jac)))
                f = _forces(ctx, x)
                fv = f.ravel()
                perp, rank, sv = perp_projector_apply(jac, fv)
                p_scan = perp_projector_apply(jac, dihedral_gradient(x, scan.quartet).ravel())[0]
                per_atom = np.linalg.norm(perp.reshape(-1, 3), axis=1)
                top = np.argsort(per_atom)[::-1][:3]
                net_f = np.linalg.norm(f.sum(axis=0))
                net_t = np.linalg.norm(np.cross(x - x.mean(axis=0), f).sum(axis=0))
                n_perp = float(np.linalg.norm(perp))
                frames.append({
                    "frame": s, "phi_deg": float(np.degrees(scan.phi_rad[s])),
                    "converged": bool(scan.converged[s]), "hysteresis": float(scan.hysteresis[s]),
                    "norm_F_fixed": float(np.linalg.norm(fv)), "norm_perp": n_perp,
                    "lower": max(0.0, n_perp - delta_max), "upper": n_perp + delta_max,
                    "max_atom_perp": float(per_atom.max()),
                    "top_atoms": ";".join(f"{int(a)}:{per_atom[a]:.1f}" for a in top),
                    "rank": rank, "n_quartets": len(quartets), "sv_min_over_max": float(sv.min() / sv.max()),
                    "perp_of_scan_grad": float(np.linalg.norm(p_scan)),
                    "net_force": float(net_f), "net_torque": float(net_t),
                })
            fpath = outdir / "frames" / f"{name}_{label}_{vname}.csv"
            fpath.parent.mkdir(parents=True, exist_ok=True)
            with open(fpath, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(frames[0]))
                w.writeheader()
                w.writerows(frames)
            conv = [fr for fr in frames if fr["converged"]]
            np_ = np.array([fr["norm_perp"] for fr in conv])
            rows.append({
                "molecule": name, "bond": label, "variant": vname, "scan_quartet": "-".join(map(str, scan.quartet)),
                "n_quartets": len(quartets), "rank_min": min(fr["rank"] for fr in frames),
                "n_frames": len(frames), "n_converged": len(conv),
                "fmax_kJ_mol_nm": fmax_kj_nm, "delta_max": delta_max,
                "perp_min": float(np_.min()), "perp_median": float(np.median(np_)), "perp_max": float(np_.max()),
                "lower_max": float(max(fr["lower"] for fr in conv)),
                "frames_lower_gt0": int(sum(fr["lower"] > 0 for fr in conv)),
                "norm_F_fixed_median": float(np.median([fr["norm_F_fixed"] for fr in conv])),
                "phi_at_perp_max_deg": conv[int(np.argmax(np_))]["phi_deg"],
                "top_atoms_at_perp_max": conv[int(np.argmax(np_))]["top_atoms"],
                "check_fd_grad_max_abs": fd_err, "check_jac_max_abs": jmax,
                "check_perp_scan_grad_max": float(max(fr["perp_of_scan_grad"] for fr in frames)),
                "check_net_force_max": float(max(fr["net_force"] for fr in frames)),
                "check_net_torque_max": float(max(fr["net_torque"] for fr in frames)),
                "sv_ratio_min": float(min(fr["sv_min_over_max"] for fr in frames)),
            })
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--g2-run", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--molecules", nargs="+", required=True)
    ap.add_argument("--decompose", action="store_true", help="只做 bond 变体的分项拆解，写 decompose/")
    args = ap.parse_args(argv)
    t0 = time.time()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in args.molecules:
        rows.extend(analyse_molecule(name, args.g2_run, args.outdir, args.decompose))
        print(f"[{time.time() - t0:7.1f}s] {name} done", flush=True)
    if args.decompose:
        print(f"decompose done, {time.time() - t0:.1f} s → {args.outdir / 'decompose'}")
        return 0
    with open(args.outdir / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    manifest = {"task": "candidate A smoke: torsion-unfixable force bound (innovation_screen_20261029 §3.1)",
                "argv": sys.argv, "host": platform.node(), "python": sys.version.split()[0],
                "wall_s": time.time() - t0, "g2_run": str(args.g2_run), "molecules": args.molecules,
                "units": "force kJ/mol/nm, angle deg"}
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"done: {len(rows)} rows, {time.time() - t0:.1f} s → {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
