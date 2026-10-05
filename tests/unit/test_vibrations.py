"""T22：validation/vibrations.py（ff_hessian、minimize_ff、kabsch_rmsd）的必测用例。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from ffcompiler.bonded.internal import linear_angles, wilson_b
from ffcompiler.importers import from_amber
from ffcompiler.ir import Angle, Atom, Bond, ForceFieldIR
from ffcompiler.validation.vibrations import (
    _forces_at,
    _reference_context,
    ff_hessian,
    kabsch_rmsd,
    minimize_ff,
)

REPO = Path(__file__).resolve().parents[2]
GAFF2 = REPO / "tests" / "fixtures" / "gaff2"


def _load_ir(name: str) -> tuple[ForceFieldIR, np.ndarray]:
    mdir = GAFF2 / name
    return from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )


def test_ff_hessian_matches_analytic() -> None:
    """用例 1：只有 bond+angle 项的 IR，在 q⁰=q(x₀) 几何上 FD == Σ k b bᵀ（≤1e-6）。"""
    pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.09572, 0.0, 0.0],
            [0.09572 * np.cos(np.radians(104.52)), 0.09572 * np.sin(np.radians(104.52)), 0.0],
        ]
    )
    coords = [("bond", (0, 1)), ("bond", (0, 2)), ("angle", (1, 0, 2))]
    ks = [4.0e5, 5.0e5, 300.0]
    b_matrix = wilson_b(pos, coords)
    h_analytic = sum(k * np.outer(b, b) for k, b in zip(ks, b_matrix))

    atoms = tuple(
        Atom(index=i, name=f"A{i}", element="O" if i == 0 else "H",
             atomic_number=8 if i == 0 else 1, mass=16.0 if i == 0 else 1.008,
             atom_type="O" if i == 0 else "H", formal_charge=0, partial_charge=0.0,
             sigma=0.3, epsilon=0.0)
        for i in range(3)
    )
    ir = ForceFieldIR(
        name="water-min", atoms=atoms,
        bonds=(
            Bond(i=0, j=1, r0=0.09572, k=ks[0], source="t"),
            Bond(i=0, j=2, r0=float(np.linalg.norm(pos[2] - pos[0])), k=ks[1], source="t"),
        ),
        angles=(Angle(i=1, j=0, k_atom=2, theta0=float(np.radians(104.52)),
                      force_constant=ks[2], source="t"),),
        proper_torsions=(), impropers=(),
        exclusions=((0, 1), (0, 2), (1, 2)), pairs14=(), total_charge=0,
        forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0, vdw_14_scale=1.0,
    )
    h_fd = ff_hessian(ir, pos)
    scale = float(np.abs(h_analytic).max())
    assert float(np.abs(h_fd - h_analytic).max()) / scale <= 1e-6


def test_exclude_bonds_angles_keeps_linear_angles() -> None:
    """exclude_terms="bonds_angles" 只去掉待推导项：线性键角留在固定部分。"""
    ir, pos = _load_ir("acetonitrile")
    angles = [(a.i, a.j, a.k_atom) for a in ir.angles]
    linear = set(linear_angles(pos, angles))
    assert linear
    expected_ir = ir.replace_terms(
        bonds=(), angles=tuple(a for a in ir.angles if (a.i, a.j, a.k_atom) in linear)
    )
    np.testing.assert_allclose(
        ff_hessian(ir, pos, exclude_terms="bonds_angles"),
        ff_hessian(expected_ir, pos),
        rtol=0, atol=1e-9,
    )


def test_minimize_ff_reaches_force_tolerance() -> None:
    """minimize_ff 返回的几何上 max|F| ≤ tolerance（未收敛必须报错，不得静默返回）。"""
    ir, positions = _load_ir("ethanol")
    tol = 1e-3
    optimized = minimize_ff(ir, positions, tolerance_kj_nm=tol)
    forces = _forces_at(_reference_context(ir), optimized)
    assert float(np.abs(forces).max()) <= tol


def test_kabsch_rmsd_rotation_and_offset() -> None:
    _, positions = _load_ir("ethanol")
    rotated = positions @ Rotation.from_euler("xyz", [30, -50, 110], degrees=True).as_matrix().T
    assert kabsch_rmsd(positions, rotated + 0.3) <= 1e-12
    displaced = positions.copy()
    displaced[0] += np.array([0.01, 0.0, 0.0])
    rmsd = kabsch_rmsd(positions, displaced)
    # 单原子位移 d：对齐后 RMSD 介于 0 与 d/√N 之间
    assert 0.0 < rmsd <= 0.01 / np.sqrt(len(positions)) + 1e-12
