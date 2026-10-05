"""Sobtop oracle 比对（MVP-1 计划 §2.8；T17–T19 测试与决策记录共用）。

`load_case` 读取 `tests/fixtures/hessian/sobtop_refs/<case>/` 的 mol2 + Hessian，
键序与 mol2 BOND 段相同（Sobtop mSeminario 的键方向依赖这个顺序）。
`oracle_mismatches` 按 §2.8 容差逐项比对 k、键长、键角，返回全部超差项。
直接运行时打印三种方法 × 全部 case 的比对表（写进决策记录 0011/0012）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from ffcompiler.bonded.drih import drih  # noqa: E402
from ffcompiler.bonded.internal import proper_quartets  # noqa: E402
from ffcompiler.bonded.m2seminario import m2seminario  # noqa: E402
from ffcompiler.bonded.mseminario import mseminario  # noqa: E402
from ffcompiler.bonded.sobtop_ref import SobtopParams, read_sobtop_params  # noqa: E402
from ffcompiler.bonded.terms import DerivedTerm  # noqa: E402
from ffcompiler.hessian.base import CartesianHessian  # noqa: E402
from ffcompiler.hessian.gaussian import read_fchk  # noqa: E402
from ffcompiler.units import HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2  # noqa: E402

REFS = REPO / "tests" / "fixtures" / "hessian" / "sobtop_refs"
KNOWN_LIMITATIONS = REPO / "tests" / "fixtures" / "hessian" / "sobtop_known_limitations"
GATE_CASES = tuple(
    sorted(p.name for p in REFS.iterdir() if p.is_dir() and p.name != "sf6_fchk")
)

# §2.8 容差（param.txt 单位：k kJ/mol/Å² 或 kJ/mol/rad²，键长 Å，键角 deg）
K_REL_TOL = 1e-4
K_ABS_TOL = 5e-3
LENGTH_TOL_A = 2e-6
ANGLE_TOL_DEG = 2e-4


def load_case(
    case: str, root: Path = REFS
) -> tuple[CartesianHessian, list[tuple[int, int]], list[tuple[int, int, int]]]:
    cdir = root / case
    bonds: list[tuple[int, int]] = []
    coords: list[list[float]] = []
    section = ""
    for line in (cdir / "input.mol2").read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue
        if not line.strip():
            continue
        p = line.split()
        if section == "@<TRIPOS>ATOM":
            coords.append([float(p[2]), float(p[3]), float(p[4])])
        elif section == "@<TRIPOS>BOND":
            bonds.append((int(p[1]) - 1, int(p[2]) - 1))
    neighbors: dict[int, set[int]] = {}
    for a, b in bonds:
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    angles = []
    for v, nbrs in neighbors.items():
        ordered = sorted(nbrs)
        for x in range(len(ordered)):
            for y in range(x + 1, len(ordered)):
                angles.append((ordered[x], v, ordered[y]))
    positions_nm = np.array(coords) * 0.1
    if (cdir / "input.fchk").is_file():
        h = read_fchk(cdir / "input.fchk")
    else:
        n3 = 3 * positions_nm.shape[0]
        values: list[float] = []
        in_block = False
        for line in (cdir / "hessian").read_text().splitlines():
            s = line.strip()
            if s.startswith("$hessian"):
                in_block = True
                continue
            if in_block and s.startswith("$"):
                break
            if in_block:
                values.extend(float(t) for t in s.split())
        hess = np.array(values).reshape(n3, n3) * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2
        h = CartesianHessian(
            atomic_numbers=np.full(positions_nm.shape[0], 6), positions_nm=positions_nm,
            hessian=hess, gradient=None, energy_kj=None, provider="xtb", method="gfn2",
        )
    return h, bonds, angles


def quartets_from(bonds: list[tuple[int, int]]) -> list[tuple[int, int, int, int]]:
    nbrs: dict[int, set[int]] = {}
    for a, b in bonds:
        nbrs.setdefault(a, set()).add(b)
        nbrs.setdefault(b, set()).add(a)
    return proper_quartets({k: frozenset(v) for k, v in nbrs.items()})


def oracle_mismatches(terms: list[DerivedTerm], params: SobtopParams) -> list[str]:
    """§2.8 逐项比对；返回超差/缺失项的描述（空列表 = 全部一致）。"""
    bad: list[str] = []
    n_bonds = n_angles = 0
    for t in terms:
        if len(t.atoms) == 2:
            n_bonds += 1
            key = tuple(sorted(t.atoms))
            if key not in params.bonds:
                bad.append(f"bond {t.atoms} 不在 param.txt")
                continue
            length_ref, k_ref = params.bonds[key]
            ours_k, ours_v, tol_v = t.k / 100.0, t.value * 10.0, LENGTH_TOL_A
            v_ref = length_ref
        else:
            n_angles += 1
            i, j, k_at = t.atoms
            key = (min(i, k_at), j, max(i, k_at))
            if key not in params.angles:
                bad.append(f"angle {t.atoms} 不在 param.txt")
                continue
            v_ref, k_ref = params.angles[key]
            ours_k, ours_v, tol_v = t.k, float(np.degrees(t.value)), ANGLE_TOL_DEG
        if abs(ours_k - k_ref) > K_REL_TOL * abs(k_ref) + K_ABS_TOL:
            bad.append(f"{t.atoms} k: ours={ours_k:.6f} sobtop={k_ref:.6f}")
        if abs(ours_v - v_ref) > tol_v:
            bad.append(f"{t.atoms} 平衡值: ours={ours_v:.6f} sobtop={v_ref:.6f}")
    if n_bonds != len(params.bonds) or n_angles != len(params.angles):
        bad.append(
            f"项数不符: ours {n_bonds} bond/{n_angles} angle,"
            f" sobtop {len(params.bonds)}/{len(params.angles)}"
        )
    return bad


def max_rel_k_error(terms: list[DerivedTerm], params: SobtopParams) -> float:
    worst = 0.0
    for t in terms:
        if len(t.atoms) == 2:
            _, k_ref = params.bonds[tuple(sorted(t.atoms))]
            ours = t.k / 100.0
        else:
            i, j, k_at = t.atoms
            _, k_ref = params.angles[(min(i, k_at), j, max(i, k_at))]
            ours = t.k
        worst = max(worst, abs(ours - k_ref) / abs(k_ref))
    return worst


METHODS = {
    "mseminario": ("param_m2.txt", lambda h, b, a: mseminario(h, b, a)),
    "m2seminario": ("param_m3.txt", lambda h, b, a: m2seminario(h, b, a)),
    "drih": ("param_m4.txt", lambda h, b, a: drih(h, b, a, dihedrals=quartets_from(b))),
}


def main() -> None:
    print("| case | 方法 | 项数 | 超差项 | 最大相对 k 偏差 |")
    print("|---|---|---|---|---|")
    for case in GATE_CASES:
        h, bonds, angles = load_case(case)
        for method, (fname, fn) in METHODS.items():
            params = read_sobtop_params(REFS / case / fname)
            terms = fn(h, bonds, angles)
            bad = oracle_mismatches(terms, params)
            print(
                f"| {case} | {method} | {len(terms)} | {len(bad)} |"
                f" {max_rel_k_error(terms, params):.1e} |"
            )


if __name__ == "__main__":
    main()
