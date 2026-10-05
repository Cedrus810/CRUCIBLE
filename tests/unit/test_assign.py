"""T21：derive_bonded 驱动函数与 CLI 的必测用例。"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.bonded.assign import BondedReport, TopologyMismatchError, derive_bonded
from ffcompiler.bonded.internal import (
    LinearAngleError,
    linear_angles,
    proper_quartets,
    wilson_b,
)
from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.hessian.gaussian import read_fchk
from ffcompiler.hessian.xtb import read_xtb_hessian
from ffcompiler.importers import from_amber
from ffcompiler.ir import ForceFieldIR
from ffcompiler.ir.topology import bond_graph

REPO = Path(__file__).resolve().parents[2]
GAFF2 = REPO / "tests" / "fixtures" / "gaff2"
XTB_FIXTURE = REPO / "tests" / "fixtures" / "hessian" / "xtb_ethanol_gfn2"
METHODS = ["mseminario", "m2seminario", "drih"]


def _load_ir(name: str) -> tuple[ForceFieldIR, np.ndarray]:
    mdir = GAFF2 / name
    return from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )


def _ethanol_hessian() -> CartesianHessian:
    return read_xtb_hessian(
        XTB_FIXTURE / "hessian", XTB_FIXTURE / "xtbopt.xyz", method="gfn2"
    )


@pytest.mark.parametrize("method", METHODS)
def test_ethanol_three_methods(method: str) -> None:
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    new_ir, hess_pos, report = derive_bonded(
        ir, h, method=method, input_positions_nm=positions
    )
    assert isinstance(report, BondedReport)
    assert new_ir.proper_torsions == ir.proper_torsions
    assert new_ir.impropers == ir.impropers
    assert new_ir.pairs14 == ir.pairs14
    for bond in new_ir.bonds:
        assert bond.source == method
        assert bond.meta["hessian_provider"] == "xtb"
    for angle in new_ir.angles:
        assert angle.source == method
    # 报告字段
    assert report.hessian_provider == "xtb"
    assert report.hessian_sha256 == h.source_sha256
    assert report.diagnostics.frequencies_cm.size == 3 * len(ir.atoms) - 6


def test_shuffled_atom_order_rejected() -> None:
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    perm = np.arange(len(h.atomic_numbers))[::-1]
    order = (3 * perm[:, None] + np.arange(3)).ravel()
    h_shuffled = CartesianHessian(
        atomic_numbers=h.atomic_numbers[perm],
        positions_nm=h.positions_nm[perm],
        hessian=h.hessian[np.ix_(order, order)],
        gradient=None, energy_kj=None, provider="xtb", method="gfn2",
    )
    with pytest.raises(TopologyMismatchError):
        derive_bonded(ir, h_shuffled, method="m2seminario")


def test_stretched_bond_rejected() -> None:
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    # 把 Hessian 几何里 C-O 键（原子 1-2）拉长到 1.5 倍
    pos = h.positions_nm.copy()
    vec = pos[2] - pos[1]
    pos[2] = pos[1] + 1.5 * vec
    h_stretched = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=pos,
        hessian=h.hessian, gradient=None, energy_kj=None,
        provider="xtb", method="gfn2",
    )
    with pytest.raises(TopologyMismatchError, match=r"1, 2"):
        derive_bonded(ir, h_stretched, method="m2seminario",
                      input_positions_nm=positions)


def _synthetic_linear_molecule_hessian(ir: ForceFieldIR, positions: np.ndarray) -> CartesianHessian:
    """在 IR 几何上构造合成 Hessian（正 k 的 b bᵀ 和），用于线性角策略测试。

    线性键角及经过它的二面角按 §2.5 规则剔除（wilson_b 对线性角会抛错）。
    """
    graph = bond_graph(ir)
    coords = [("bond", (b.i, b.j)) for b in ir.bonds]
    coords += [("angle", (a.i, a.j, a.k_atom)) for a in ir.angles]
    coords += [("dihedral", q) for q in proper_quartets(graph)]

    def canonical(t):
        i, j, k = t
        return (min(i, k), j, max(i, k))

    lin = {canonical(a) for a in linear_angles(positions, [c[1] for c in coords if c[0] == "angle"])}
    coords = [
        c for c in coords
        if not (c[0] == "angle" and canonical(c[1]) in lin)
        and not (c[0] == "dihedral" and canonical(c[1][:3]) in lin)
    ]
    b_matrix = wilson_b(positions, coords)
    ks = np.full(len(coords), 1000.0)
    hessian = sum(k * np.outer(b, b) for k, b in zip(ks, b_matrix))
    return CartesianHessian(
        atomic_numbers=np.array([a.atomic_number for a in ir.atoms]),
        positions_nm=positions,
        hessian=hessian,
        gradient=np.zeros(hessian.shape[0]),
        energy_kj=None, provider="test", method="synthetic",
    )


@pytest.mark.parametrize("method", METHODS)
def test_acetonitrile_linear_angle_policies(method: str) -> None:
    ir, positions = _load_ir("acetonitrile")
    h = _synthetic_linear_molecule_hessian(ir, positions)
    with pytest.raises(LinearAngleError):
        derive_bonded(ir, h, method=method)

    new_ir, _, report = derive_bonded(
        ir, h, method=method, linear_angle_policy="library"
    )
    assert len(report.linear_angles_kept) == 1
    i, j, k_at = report.linear_angles_kept[0]
    original = next(
        a for a in ir.angles
        if (a.i, a.j, a.k_atom) == (i, j, k_at)
    )
    kept = next(
        a for a in new_ir.angles
        if (a.i, a.j, a.k_atom) == (i, j, k_at)
    )
    assert kept.theta0 == original.theta0
    assert kept.force_constant == original.force_constant
    assert kept.source == original.source
    assert kept.meta["linear_angle"] is True
    for angle in new_ir.angles:
        if (angle.i, angle.j, angle.k_atom) != (i, j, k_at):
            assert angle.source == method


def test_symmetrize_on_off() -> None:
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    new_on, _, report_on = derive_bonded(
        ir, h, method="m2seminario", symmetrize=True, input_positions_nm=positions
    )
    ch3_h = [a.index for a in ir.atoms if a.element == "H"][:3]  # 不严谨，用图判定
    graph = bond_graph(ir)
    carbon0 = next(i for i in graph if len(graph[i]) == 4)
    ch3_h = sorted(j for j in graph[carbon0] if ir.atoms[j].element == "H")
    ks_on = [
        b.k for b in new_on.bonds
        if {b.i, b.j} <= set(ch3_h) | {carbon0} and carbon0 in (b.i, b.j)
        and ir.atoms[[b.i, b.j][1 if b.i == carbon0 else 0]].element == "H"
    ]
    assert len(ks_on) == 3
    assert len(set(ks_on)) == 1
    assert report_on.class_spread

    new_off, _, report_off = derive_bonded(
        ir, h, method="m2seminario", symmetrize=False, input_positions_nm=positions
    )
    ks_off = [
        b.k for b in new_off.bonds
        if carbon0 in (b.i, b.j)
        and ir.atoms[b.j if b.i == carbon0 else b.i].element == "H"
        and {b.i, b.j} <= set(ch3_h) | {carbon0}
    ]
    assert len(ks_off) == 3
    assert len(set(ks_off)) > 1  # 不平均：Hessian 里三个 C-H 并不严格等价
    ch3_spread = {
        b.meta["class_spread"] for b in new_off.bonds
        if carbon0 in (b.i, b.j) and ir.atoms[b.j if b.i == carbon0 else b.i].element == "H"
    }
    assert len(ch3_spread) == 1
    spread = ch3_spread.pop()
    assert spread == pytest.approx((max(ks_off) - min(ks_off)) / abs(np.mean(ks_off)))
    assert spread in report_off.class_spread.values()  # 不平均也记录离散度


def test_allow_nonstationary_reaches_drih() -> None:
    """allow_nonstationary=True 必须一路传到 drih（drih 内部也做驻点检查）。"""
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    h_grad = CartesianHessian(
        atomic_numbers=h.atomic_numbers, positions_nm=h.positions_nm,
        hessian=h.hessian, gradient=np.full(h.hessian.shape[0], 100.0),
        energy_kj=None, provider="xtb", method="gfn2",
    )
    derive_bonded(ir, h_grad, method="drih", allow_nonstationary=True,
                  input_positions_nm=positions)


def test_diagnostics_use_ir_masses() -> None:
    """§2.2：有 IR 时频率用 IR 的 mass。"""
    from ffcompiler.hessian.base import harmonic_frequencies

    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    _, _, report = derive_bonded(ir, h, method="m2seminario", input_positions_nm=positions)
    # GAFF2 质量（如 O 16.00）与周期表平均原子量（O 15.999）不同，频率差 ~1e-2 cm⁻¹
    expected = harmonic_frequencies(h.hessian, h.positions_nm, np.array([a.mass for a in ir.atoms]))
    np.testing.assert_allclose(report.diagnostics.frequencies_cm, expected, rtol=0, atol=1e-9)


def test_term_meta_class_spread_is_per_class() -> None:
    """T21 规格 3：每项 meta 的 class_spread 是它所在等价类的离散度。"""
    ir, positions = _load_ir("ethanol")
    h = _ethanol_hessian()
    new_ir, _, report = derive_bonded(
        ir, h, method="m2seminario", symmetrize=False, input_positions_nm=positions
    )
    graph = bond_graph(ir)
    oxygen = next(a.index for a in ir.atoms if a.element == "O")
    oh = next(b for b in new_ir.bonds if oxygen in (b.i, b.j)
              and ir.atoms[b.j if b.i == oxygen else b.i].element == "H")
    assert oh.meta["class_spread"] == 0.0  # O-H 单独成类
    carbon0 = next(i for i in graph if len(graph[i]) == 4
                   and sum(ir.atoms[j].element == "H" for j in graph[i]) == 3)
    ch3 = [b for b in new_ir.bonds if carbon0 in (b.i, b.j)
           and ir.atoms[b.j if b.i == carbon0 else b.i].element == "H"]
    spreads = {b.meta["class_spread"] for b in ch3}
    assert len(spreads) == 1 and spreads.pop() > 0.0
    assert max(report.class_spread.values()) > oh.meta["class_spread"]


@pytest.mark.external
def test_cli_derive_bonded_ethanol_gfnff(tmp_path: Path) -> None:
    mdir = GAFF2 / "ethanol"
    cmd = [
        sys.executable, "-m", "ffcompiler.cli", "derive-bonded",
        "--prmtop", str(mdir / "ethanol.prmtop"),
        "--inpcrd", str(mdir / "ethanol.inpcrd"),
        "--sdf", str(mdir / "ethanol.sdf"),
        "--hessian-provider", "gfnff",
        "--method", "m2seminario",
        "--output", "gromacs", "amber", "openmm",
        "--outdir", str(tmp_path / "out"),
    ]
    result = subprocess.run(
        cmd, cwd=REPO, env={"PYTHONPATH": str(REPO / "src"), "PATH": "/usr/bin:/bin"},
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    out = tmp_path / "out"
    assert (out / "gromacs" / "topol.top").is_file()
    assert list((out / "amber").glob("*.prmtop"))
    assert list((out / "openmm").glob("*.xml"))


def test_cli_derive_bonded_rejects_unconverged_xtb(tmp_path: Path) -> None:
    """CLI 必须走 run_xtb_hessian 的成功条件：没有 CONVERGED 行 → 退出码非 0、不写输出。

    伪造的 xtb 退出码 0、给出 hessian/xtbopt.xyz，但不打印收敛行；同时放在
    FFC_XTB 和 PATH 上，确保 CLI 不论怎么选二进制都会用到它。
    """
    fake_dir = tmp_path / "bin"
    fake_dir.mkdir()
    fake = fake_dir / "xtb"
    fake.write_text(
        "#!/bin/bash\n"
        f"cp {XTB_FIXTURE / 'hessian'} hessian\n"
        f"cp {XTB_FIXTURE / 'xtbopt.xyz'} xtbopt.xyz\n"
        "echo 'fake xtb: no convergence line'\n"
    )
    fake.chmod(0o755)
    mdir = GAFF2 / "ethanol"
    env = dict(__import__("os").environ, FFC_XTB=str(fake),
               PATH=f"{fake_dir}:{__import__('os').environ['PATH']}",
               PYTHONPATH=str(REPO / "src"))
    outdir = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, "-m", "ffcompiler.cli", "derive-bonded",
         "--prmtop", str(mdir / "ethanol.prmtop"), "--inpcrd", str(mdir / "ethanol.inpcrd"),
         "--sdf", str(mdir / "ethanol.sdf"), "--hessian-provider", "gfn2",
         "--method", "m2seminario", "--output", "gromacs", "--outdir", str(outdir)],
        capture_output=True, text=True, env=env,
    )
    assert proc.returncode != 0, proc.stdout
    assert "CONVERGED" in proc.stderr
    assert not (outdir / "gromacs").exists()


def test_nonpositive_force_constant_reported_clearly() -> None:
    """决策 0014：方法给出非正 k 时抛 NonPositiveForceConstantError，列出方法与项，
    提示换方法；不自动回退到库参数。"""
    from ffcompiler.bonded.assign import NonPositiveForceConstantError

    fx = REPO / "tests" / "fixtures" / "hessian" / "sobtop_known_limitations" / "cyclopropane_gfnff"
    ir, positions = _load_ir("cyclopropane")
    h = read_xtb_hessian(fx / "hessian", fx / "xtbopt.xyz", method="gfnff")
    with pytest.raises(NonPositiveForceConstantError, match="m2seminario") as info:
        derive_bonded(ir, h, method="m2seminario", input_positions_nm=positions)
    err = info.value
    assert err.method == "m2seminario"
    assert sorted(atoms for atoms, _ in err.terms) == [(0, 1, 2), (0, 2, 1), (1, 0, 2)]
    assert all(k <= 0 for _, k in err.terms)
    assert "mseminario" in str(err) and "drih" in str(err)
    # 其他方法在同一 Hessian 上正常
    derive_bonded(ir, h, method="mseminario", input_positions_nm=positions)
