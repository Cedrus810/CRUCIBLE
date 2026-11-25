"""T13：xtb 运行器与读取器的必测用例。

外部测试（标 external）在本机 xtb 6.7.1 存在时可真实运行；
计划测试用例 5 的奇电子报错由运行器侧校验实现（xtb 6.7.1 在 Fermi
展宽下对奇电子体系不报错，见 run_xtb_hessian 的 _check_charge_and_electrons）。
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.hessian.base import check_hessian, harmonic_frequencies
from ffcompiler.hessian.xtb import (
    XtbError,
    _check_success,
    read_xtb_gradient,
    read_xtb_hessian,
    read_xtb_vibspectrum,
    run_xtb_hessian,
    run_xtb_singlepoint,
)

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "hessian" / "xtb_ethanol_gfn2"
GAFF2 = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _fixture_hessian():
    return read_xtb_hessian(
        FIXTURE / "hessian", FIXTURE / "xtbopt.xyz", method="gfn2", version="6.7.1"
    )


def test_fixture_frequencies_match_vibspectrum() -> None:
    h = _fixture_hessian()
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    masses = np.array(
        [table.GetAtomicWeight(int(z)) for z in h.atomic_numbers], dtype=np.float64
    )
    computed = harmonic_frequencies(h.hessian, h.positions_nm, masses)
    reference = read_xtb_vibspectrum(FIXTURE / "vibspectrum")
    assert computed.size == reference.size
    assert float(np.abs(np.sort(computed) - np.sort(reference)).max()) <= 1.0


def test_hessian_missing_value_raises(tmp_path: Path) -> None:
    text = (FIXTURE / "hessian").read_text().splitlines()
    # 删掉 $hessian 块中间的一个数
    removed = False
    for idx, line in enumerate(text):
        stripped = line.strip()
        if not stripped or stripped.startswith("$"):
            continue
        tokens = stripped.split()
        if len(tokens) > 1:
            tokens.pop()
            text[idx] = "  " + "  ".join(tokens)
            removed = True
            break
    assert removed
    bad = tmp_path / "hessian"
    bad.write_text("\n".join(text) + "\n")
    with pytest.raises(XtbError, match="个数"):
        read_xtb_hessian(bad, FIXTURE / "xtbopt.xyz", method="gfn2")


def test_missing_converged_flag_raises(tmp_path: Path) -> None:
    workdir = tmp_path / "run"
    workdir.mkdir()
    shutil.copy2(FIXTURE / "hessian", workdir / "hessian")
    shutil.copy2(FIXTURE / "xtbopt.xyz", workdir / "xtbopt.xyz")
    tampered = [
        line
        for line in (FIXTURE / "xtb.out").read_text().splitlines()
        if "GEOMETRY OPTIMIZATION CONVERGED" not in line
    ]
    (workdir / "xtb.out").write_text("\n".join(tampered) + "\n")
    with pytest.raises(XtbError, match="CONVERGED"):
        _check_success(workdir, 0)


@pytest.mark.external
def test_run_ethanol_gfnff_and_gfn2(tmp_path: Path) -> None:
    from ffcompiler.importers import from_amber

    mdir = GAFF2 / "ethanol"
    ir, pos = from_amber(
        mdir / "ethanol.prmtop", mdir / "ethanol.inpcrd", mdir / "ethanol.sdf"
    )
    z = np.array([a.atomic_number for a in ir.atoms])
    for method in ("gfnff", "gfn2"):
        h = run_xtb_hessian(
            z, pos, method=method, charge=ir.total_charge,
            workdir=tmp_path / method, threads=2,
        )
        diag = check_hessian(h)
        assert diag.frequencies_cm.min() > -20.0
        assert h.meta["stationary_by"] == "xtb ohess vtight"


def _parse_charge(xtb_out: Path) -> tuple[int, float]:
    text = xtb_out.read_text()
    net = re.findall(r"net charge\s+(-?\d+)", text)
    total = re.findall(r"total charge\s+(-?[\d.]+)\s+e", text)
    assert net and total
    return int(net[0]), float(total[0])


@pytest.mark.external
@pytest.mark.parametrize(
    "name,charge", [("acetate", -1), ("methylammonium", 1)]
)
def test_run_charged_molecules(tmp_path: Path, name: str, charge: int) -> None:
    from ffcompiler.importers import from_amber

    mdir = GAFF2 / name
    ir, pos = from_amber(
        mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
    )
    assert ir.total_charge == charge
    z = np.array([a.atomic_number for a in ir.atoms])
    h = run_xtb_hessian(
        z, pos, method="gfn2", charge=charge, workdir=tmp_path / name, threads=2
    )
    net, total = _parse_charge(Path(h.source_file).parent / "xtb.out")
    assert net == charge
    assert total == pytest.approx(float(charge), abs=1e-6)


@pytest.mark.external
def test_wrong_charge_rejected(tmp_path: Path) -> None:
    from ffcompiler.importers import from_amber

    mdir = GAFF2 / "acetate"
    ir, pos = from_amber(
        mdir / "acetate.prmtop", mdir / "acetate.inpcrd", mdir / "acetate.sdf"
    )
    z = np.array([a.atomic_number for a in ir.atoms])
    # 醋酸根给 charge=0 → 价电子数为奇数，运行器在启动 xtb 前拒绝
    with pytest.raises(XtbError, match="奇数"):
        run_xtb_hessian(
            z, pos, method="gfn2", charge=0, workdir=tmp_path / "odd", threads=2
        )


_VIB_HEADER = (
    "$vibrational spectrum\n"
    "#  mode     symmetry     wave number   IR intensity    selection rules\n"
    "#                         cm**(-1)      (km*mol⁻¹)        IR     \n"
)


def test_vibspectrum_without_symmetry_label(tmp_path: Path) -> None:
    # 无对称标签的行：mode 后面直接就是波数，不能把 IR 强度当成波数
    rows = [f"     {i}                       0.00         0.00000          - " for i in range(1, 7)]
    rows += [
        "     7            225.65         2.67660         YES",
        "     8            340.99       140.84247         YES",
        "     9        a   420.78        30.75955         YES",
    ]
    path = tmp_path / "vibspectrum"
    path.write_text(_VIB_HEADER + "\n".join(rows) + "\n$end\n")
    assert read_xtb_vibspectrum(path).tolist() == pytest.approx([225.65, 340.99, 420.78])


def test_vibspectrum_removes_exactly_six_zero_modes(tmp_path: Path) -> None:
    # 零模未完全投影干净（残余 |ν| ≥ 0.5）：仍应恰好去掉 |ν| 最小的 6 个
    zeros = [-1.20, 0.80, 0.00, 0.55, -0.70, 1.10]
    reals = [-35.0, 225.65, 340.99]  # 虚频（负值）保留
    rows = [f"     {i + 1}        a   {w:10.2f}     0.00000          - " for i, w in enumerate(zeros)]
    rows += [
        f"     {i + 7}        a   {w:10.2f}     1.00000         YES"
        for i, w in enumerate(reals)
    ]
    path = tmp_path / "vibspectrum"
    path.write_text(_VIB_HEADER + "\n".join(rows) + "\n$end\n")
    assert read_xtb_vibspectrum(path).tolist() == pytest.approx(reals)


def test_vibspectrum_linear_removes_five(tmp_path: Path) -> None:
    # 线性分子（CO2，3N = 9）：5 个零模，4 个振动
    zeros = [0.0] * 5
    reals = [650.0, 650.0, 1350.0, 2400.0]
    rows = [f"     {i + 1}        {w:10.2f}     0.00000          - " for i, w in enumerate(zeros)]
    rows += [f"     {i + 6}   {w:10.2f}     1.00000         YES" for i, w in enumerate(reals)]
    path = tmp_path / "vibspectrum"
    path.write_text(_VIB_HEADER + "\n".join(rows) + "\n$end\n")
    assert read_xtb_vibspectrum(path, linear=True).tolist() == pytest.approx(reals)
    # 默认按非线性去掉 6 个
    assert read_xtb_vibspectrum(path).tolist() == pytest.approx(reals[1:])


def test_vibspectrum_bad_mode_count_raises(tmp_path: Path) -> None:
    rows = [f"     {i}        0.00     0.00000          - " for i in range(1, 6)]
    path = tmp_path / "vibspectrum"
    path.write_text(_VIB_HEADER + "\n".join(rows) + "\n$end\n")
    with pytest.raises(ValueError, match="3N"):
        read_xtb_vibspectrum(path)


# ---- xtb 单点（覆盖集诊断 xtb_sp 阶段，2026-11-20 计划 §2.2 第 4 项） ----

# 取自 xtb 6.7.1 `--gfn 2 --grad` 对水分子的真实输出（2026-11-24 本机）
_WATER_GRADIENT = """$grad
  cycle =      1    SCF energy =    -5.07037076185   |dE/dxyz| =  0.007222
    0.00000000000000      0.00000000000000      0.22166487441798      O 
    0.00000000000000      1.43090062156260     -0.88665949767191      H 
    0.00000000000000     -1.43090062156260     -0.88665949767191      H 
  -2.6881125797143E-18  -1.1519142618602E-16   2.9402221161207E-03
   1.6496056751868E-18  -4.4262812903950E-03  -1.4701110580603E-03
   1.0385069045275E-18   4.4262812903951E-03  -1.4701110580603E-03
$end
"""


def test_read_xtb_gradient_water(tmp_path: Path) -> None:
    path = tmp_path / "gradient"
    path.write_text(_WATER_GRADIENT)
    energy, grad = read_xtb_gradient(path, 3)
    assert energy == -5.07037076185
    assert grad.shape == (3, 3)
    assert grad[0, 2] == 2.9402221161207e-03
    assert grad[1, 1] == -4.4262812903950e-03
    # 平移不变：梯度和为 0（量级 1e-18）
    assert np.abs(grad.sum(axis=0)).max() < 1e-15


def test_read_xtb_gradient_takes_last_cycle_and_d_exponent(tmp_path: Path) -> None:
    second = """  cycle =      2    SCF energy =    -5.10000000000D+00   |dE/dxyz| =  0.000001
    0.0 0.0 0.2 O
    0.0 1.4 -0.9 H
    0.0 -1.4 -0.9 H
   1.0D-03  2.0D-03  3.0D-03
   4.0D-03  5.0D-03  6.0D-03
   7.0D-03  8.0D-03  9.0D-03
$end
"""
    path = tmp_path / "gradient"
    path.write_text(_WATER_GRADIENT.replace("$end\n", "") + second)
    energy, grad = read_xtb_gradient(path, 3)
    assert energy == -5.1
    np.testing.assert_array_equal(grad[2], [7.0e-3, 8.0e-3, 9.0e-3])


def test_read_xtb_gradient_wrong_atom_count_raises(tmp_path: Path) -> None:
    path = tmp_path / "gradient"
    path.write_text(_WATER_GRADIENT)
    with pytest.raises(XtbError, match="2×4"):
        read_xtb_gradient(path, 4)


def test_read_xtb_gradient_not_gradient_file_raises(tmp_path: Path) -> None:
    path = tmp_path / "gradient"
    path.write_text("$energy\n 1 -5.0\n$end\n")
    with pytest.raises(XtbError, match="cycle"):
        read_xtb_gradient(path, 3)


def test_singlepoint_odd_electrons_rejected_before_launch(tmp_path: Path) -> None:
    # OH 中性：价电子 7 → 奇数；binary 指向不存在的路径，证明拒绝发生在启动 xtb 之前
    z = np.array([8, 1])
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.097]])
    with pytest.raises(XtbError, match="奇数"):
        run_xtb_singlepoint(z, pos, method="gfn2", charge=0, workdir=tmp_path,
                            binary=tmp_path / "no-such-xtb")


@pytest.mark.external
def test_run_singlepoint_ethanol_gfn2(tmp_path: Path) -> None:
    from ffcompiler.importers import from_amber

    mdir = GAFF2 / "ethanol"
    ir, pos = from_amber(
        mdir / "ethanol.prmtop", mdir / "ethanol.inpcrd", mdir / "ethanol.sdf"
    )
    z = np.array([a.atomic_number for a in ir.atoms])
    out = run_xtb_singlepoint(z, pos, method="gfn2", charge=ir.total_charge,
                              workdir=tmp_path, threads=2)
    assert out["gradient_hartree_per_bohr"].shape == (len(z), 3)
    assert np.isfinite(out["energy_hartree"]) and out["energy_hartree"] < 0
    assert np.abs(out["gradient_hartree_per_bohr"].sum(axis=0)).max() < 1e-6
    assert out["xtb_version"] == "6.7.1"
    # 重跑同一目录：旧 gradient 被清掉，结果逐位相同（不追加 cycle）
    again = run_xtb_singlepoint(z, pos, method="gfn2", charge=ir.total_charge,
                                workdir=tmp_path, threads=2)
    assert again["energy_hartree"] == out["energy_hartree"]


# ---- 鞍点重启（2026-11-25；诊断 experiments/20261125_saddle_diagnosis/） ----

def _fake_ohess_factory(lowest_freqs: list[float], tmp_path: Path):
    """替身 _ohess_once：第 i 次调用写出 vibspectrum（最低频率 lowest_freqs[i]）、xtbopt.xyz、
    以及虚频时的 xtbhess.xyz；返回带标记的 CartesianHessian。"""
    import ffcompiler.hessian.xtb as xtb_mod

    h0 = _fixture_hessian()
    n = len(h0.atomic_numbers)
    calls: list[dict] = []

    def fake(binary, z, pos, *, method, charge, workdir, threads):
        i = len(calls)
        calls.append({"workdir": Path(workdir), "pos": np.array(pos)})
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        rows = ["$vibrational spectrum"] + [f"  {k + 1}  0.00  0.0  -" for k in range(6)]
        freqs = [lowest_freqs[i]] + [100.0 + k for k in range(3 * n - 7)]
        rows += [f"  {k + 7}  a  {f:.2f}  1.0  YES" for k, f in enumerate(freqs)] + ["$end"]
        (workdir / "vibspectrum").write_text("\n".join(rows) + "\n")
        xyz = [str(n), f" energy: {-1.0 - 0.001 * i:.6f} gnorm: 0.0001"] + \
              [f"C {0.1 * k:.4f} {0.2 + i:.4f} 0.0" for k in range(n)]
        (workdir / "xtbopt.xyz").write_text("\n".join(xyz) + "\n")
        if lowest_freqs[i] < -5.0:
            (workdir / "xtbhess.xyz").write_text("\n".join(xyz) + "\n")
        return xtb_mod.replace(h0, source_file=str(workdir / "hessian"))

    return fake, calls


def _run_with_fake(monkeypatch, tmp_path, lowest, **kw):
    import ffcompiler.hessian.xtb as xtb_mod

    fake, calls = _fake_ohess_factory(lowest, tmp_path)
    monkeypatch.setattr(xtb_mod, "_ohess_once", fake)
    h0 = _fixture_hessian()
    h = run_xtb_hessian(h0.atomic_numbers, h0.positions_nm, method="gfnff", charge=0,
                        workdir=tmp_path, binary=tmp_path / "xtb", **kw)
    return h, calls


def test_no_restart_when_no_imaginary(tmp_path, monkeypatch):
    h, calls = _run_with_fake(monkeypatch, tmp_path, [150.0])
    assert len(calls) == 1
    assert dict(h.meta) == {"stationary_by": "xtb ohess vtight"}     # 不加键，逐位同改动前


def test_no_restart_for_small_imaginary_above_gate(tmp_path, monkeypatch):
    # −10 cm⁻¹：xtb 自己会写 xtbhess.xyz，但 check_hessian 只警告不拒 → 不重启、行为不变
    h, calls = _run_with_fake(monkeypatch, tmp_path, [-10.0])
    assert len(calls) == 1 and "imaginary_restarts" not in h.meta


def test_restart_from_xtbhess_then_success(tmp_path, monkeypatch):
    h, calls = _run_with_fake(monkeypatch, tmp_path, [-56.0, 19.7])
    assert len(calls) == 2
    assert calls[1]["workdir"] == tmp_path / "imag_restart_1"
    # 第二次的起点 = 第一次的 xtbhess.xyz（Å → nm）
    assert calls[1]["pos"][0, 1] == pytest.approx(0.02)
    rec = h.meta["imaginary_restarts"]
    assert len(rec) == 1 and rec[0]["round"] == 1
    assert rec[0]["lowest_freq_before_cm"] == -56.0
    assert rec[0]["energy_before_hartree"] == -1.0 and rec[0]["energy_after_hartree"] == -1.001
    assert h.source_file == str(tmp_path / "imag_restart_1" / "hessian")


def test_restart_stops_at_limit_and_returns_last(tmp_path, monkeypatch):
    h, calls = _run_with_fake(monkeypatch, tmp_path, [-60.0, -50.0, -40.0, -30.0],
                              max_imaginary_restarts=2)
    assert len(calls) == 3                                          # 原始 + 2 次重启
    assert [r["round"] for r in h.meta["imaginary_restarts"]] == [1, 2]
    assert h.source_file == str(tmp_path / "imag_restart_2" / "hessian")


def test_restart_disabled_with_zero(tmp_path, monkeypatch):
    h, calls = _run_with_fake(monkeypatch, tmp_path, [-60.0], max_imaginary_restarts=0)
    assert len(calls) == 1 and "imaginary_restarts" not in h.meta


@pytest.mark.external
def test_acetone_gfnff_saddle_restart_passes_quality_gate(tmp_path):
    """覆盖集 acetone：GFN-FF --ohess 停在 C2v 鞍点（2 虚频），重启后无虚频、过 check_hessian。"""
    from rdkit import Chem

    sdf = Path(__file__).resolve().parents[2] / "benchmarks" / "coverage_set" / "sdf" / "acetone.sdf"
    mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
    z = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    pos = mol.GetConformer().GetPositions() * 0.1
    h = run_xtb_hessian(z, pos, method="gfnff", charge=0, workdir=tmp_path, threads=2)
    assert h.meta["imaginary_restarts"][0]["lowest_freq_before_cm"] < -20.0
    e = h.meta["imaginary_restarts"][-1]
    assert e["energy_after_hartree"] < e["energy_before_hartree"]
    diag = check_hessian(h)
    assert diag.n_imaginary == 0
