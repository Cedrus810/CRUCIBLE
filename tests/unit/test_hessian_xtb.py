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
    read_xtb_hessian,
    read_xtb_vibspectrum,
    run_xtb_hessian,
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
