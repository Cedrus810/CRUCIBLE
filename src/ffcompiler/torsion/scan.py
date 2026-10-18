"""ASE 约束松弛二面角扫描（MVP-2 计划 §2.2，T26）。

对任意 ASE 计算器：先无约束松弛得 x0、φ0；网格 φ_s = φ0 + step·s；正反两遍顺序扫描，
每步在上一步几何上把 k 侧刚性转到 φ_s，再用 FixInternals 硬约束该二面角做 LBFGS 松弛。
每个网格点取两遍中能量较低的几何。能量单位 kJ/mol，坐标 nm，角度 rad（IR 定义）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import LBFGS

from ..ir.energy import dihedral_angle
from ..units import ANGSTROM_NM, EV_KJ_MOL


class ScanError(RuntimeError):
    """参考势求值抛错或给出非有限能量/力。"""


@dataclass
class ScanResult:
    quartet: tuple[int, int, int, int]
    phi_target_rad: np.ndarray      # (n,)
    phi_rad: np.ndarray             # (n,)
    energy: np.ndarray              # (n,) kJ/mol，参考势本身（硬约束无约束能）
    positions: np.ndarray           # (n, N, 3) nm
    hysteresis: np.ndarray          # (n,) |E_正 − E_反|，s = 0 为 0
    converged: np.ndarray           # (n,) bool
    max_bond_change_pm: float
    params: dict = field(default_factory=dict)


def _rotating_side(quartet, bonds, n_atoms: int) -> list[bool]:
    _, j, k, _ = quartet
    adj: dict[int, set[int]] = {a: set() for a in range(n_atoms)}
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    side, stack = {k}, [k]
    while stack:
        a = stack.pop()
        for b in adj[a]:
            if b in side or (a == k and b == j):
                continue
            side.add(b)
            stack.append(b)
    if j in side:
        raise ScanError(f"键 {j}-{k} 在环内，无法刚性转动")
    return [a in side for a in range(n_atoms)]


def _phi(pos_angstrom: np.ndarray, quartet) -> float:
    return dihedral_angle(*(pos_angstrom[a] for a in quartet))


def relaxed_scan(
    calc,
    numbers: Sequence[int],
    positions_nm: np.ndarray,
    quartet: tuple[int, int, int, int],
    bonds: Sequence[tuple[int, int]],
    *,
    info: Mapping | None = None,
    step_deg: float = 10.0,
    fmax: float = 0.005,
    max_steps: int = 2000,
) -> ScanResult:
    n_grid = round(360.0 / step_deg)
    if not math.isclose(n_grid * step_deg, 360.0):
        raise ValueError(f"step_deg={step_deg} 不能整除 360")
    quartet = tuple(int(a) for a in quartet)
    atoms = Atoms(numbers=list(numbers), positions=np.asarray(positions_nm) / ANGSTROM_NM)
    atoms.info.update(dict(info or {}))
    atoms.calc = calc
    mask = _rotating_side(quartet, bonds, len(atoms))
    bond_idx = np.asarray(bonds, dtype=int)

    def relax(where: str) -> tuple[bool, float]:
        try:
            opt = LBFGS(atoms, logfile=None)
            conv = bool(opt.run(fmax=fmax, steps=max_steps))
            e = atoms.get_potential_energy() * EV_KJ_MOL
            f = atoms.get_forces()
        except ScanError:
            raise
        except Exception as exc:  # 参考势求值失败：带上位置信息重抛
            raise ScanError(f"四元组 {quartet} {where}：参考势求值失败：{exc!r}") from exc
        if not (math.isfinite(e) and np.all(np.isfinite(f))):
            raise ScanError(f"四元组 {quartet} {where}：能量或力非有限（E={e}）")
        return conv, e

    conv0, e0 = relax("起点")
    x0 = atoms.positions.copy()
    phi0_ase = atoms.get_dihedral(*quartet)
    r0 = np.linalg.norm(x0[bond_idx[:, 0]] - x0[bond_idx[:, 1]], axis=1)

    passes: list[dict[int, tuple]] = []
    for direction in (+1, -1):
        atoms.set_constraint()
        atoms.positions = x0.copy()
        got: dict[int, tuple] = {}
        for m in range(1, n_grid):
            s = (direction * m) % n_grid
            target = (phi0_ase + step_deg * s) % 360.0
            atoms.set_constraint()
            atoms.set_dihedral(*quartet, target, mask=mask)
            phi_target = _phi(atoms.positions, quartet)
            atoms.set_constraint(FixInternals(dihedrals_deg=[[atoms.get_dihedral(*quartet), list(quartet)]]))
            conv, e = relax(f"网格点 s={s}")
            got[s] = (e, atoms.positions.copy(), phi_target, conv)
        atoms.set_constraint()
        passes.append(got)

    energy = np.empty(n_grid)
    phi = np.empty(n_grid)
    phi_t = np.empty(n_grid)
    pos = np.empty((n_grid, len(atoms), 3))
    hyst = np.zeros(n_grid)
    converged = np.empty(n_grid, dtype=bool)
    phi_x0 = _phi(x0, quartet)
    energy[0], pos[0], phi[0], phi_t[0], converged[0] = e0, x0, phi_x0, phi_x0, conv0
    for s in range(1, n_grid):
        fwd, bwd = passes[0][s], passes[1][s]
        cands = [c for c in (fwd, bwd) if c[3]] or [fwd, bwd]
        best = min(cands, key=lambda c: c[0])
        energy[s] = best[0]
        pos[s] = best[1]
        phi_t[s] = best[2]
        phi[s] = _phi(best[1], quartet)
        converged[s] = best[3]
        hyst[s] = abs(fwd[0] - bwd[0])
    r_all = np.linalg.norm(pos[:, bond_idx[:, 0]] - pos[:, bond_idx[:, 1]], axis=2)
    return ScanResult(
        quartet=quartet,
        phi_target_rad=phi_t,
        phi_rad=phi,
        energy=energy,
        positions=pos * ANGSTROM_NM,
        hysteresis=hyst,
        converged=converged,
        max_bond_change_pm=float(np.max(np.abs(r_all - r0)) * 100.0),
        params={"step_deg": float(step_deg), "fmax_eV_per_A": float(fmax), "max_steps": int(max_steps)},
    )
