"""ORCA .hess 读取器的必测用例。

主用例是 ORCA 6.1.1 实际写出的文件（``tests/fixtures/hessian/orca/water.hess``，
HF/STO-3G 水，输入见同目录 water.inp；生成记录
``experiments/20261017_orca_hess_format/``）：用读到的 Hessian 和 ORCA 自己的
质量算谐振频率，与文件里 ORCA 给的 ``$vibrational_frequencies`` 比较。另用
SF6 fchk 固件按同一真实布局（按列分块满阵）写出再读回，覆盖 5 原子以上的多块
情形；其余用例覆盖布局错误必须显式报错。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.hessian.base import check_hessian, harmonic_frequencies
from ffcompiler.hessian.gaussian import read_fchk
from ffcompiler.hessian.orca import read_orca_hess
from ffcompiler.units import (
    BOHR_NM,
    HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2 as HA_BOHR2_TO_KJ_NM2,
)

HESS_FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "hessian"
SOBTOP = HESS_FIXTURES / "sobtop_examples"
ORCA_WATER = HESS_FIXTURES / "orca" / "water.hess"


def _blocked(matrix: np.ndarray, width: int = 5) -> list[str]:
    """按 ORCA 的写法输出 $hessian 节体：维度行 + 每 width 列一块（列号表头 + 行号 值…）。"""
    n = matrix.shape[0]
    lines = [str(n)]
    for start in range(0, n, width):
        cols = range(start, min(start + width, n))
        lines.append("  ".join(str(c) for c in cols))
        for r in range(n):
            lines.append(f"{r} " + " ".join(f"{matrix[r, c]:.17e}" for c in cols))
    return lines


def _write_hess(tmp_path: Path, hess_lines: list[str], atom_lines: list[str],
                name: str = "mini.hess", header: str = "$orca_hessian_file") -> Path:
    text = "\n".join(
        [header, "", "$act_atom", "0", "$hessian", *hess_lines,
         "$atoms", str(len(atom_lines)), *atom_lines, "$end"]
    ) + "\n"
    path = tmp_path / name
    path.write_text(text)
    return path


def _orca_frequencies(path: Path) -> np.ndarray:
    lines = [l.strip() for l in path.read_text().splitlines()]
    start = lines.index("$vibrational_frequencies")
    n = int(lines[start + 1])
    return np.array([float(lines[start + 2 + i].split()[1]) for i in range(n)])


def test_real_orca_water_frequencies_match_orca() -> None:
    """真实 ORCA 文件：频率与 ORCA 自己的结果相对差 ≤ 1e-6。

    依据：文件里 Hessian 有 11 位有效数字（相对 1e-10）；我们与 ORCA 的
    Hartree/Bohr 常数不同，实测三个模式统一偏 5.74e-8（常数比，非布局误差）。
    1e-6 约为实测的 17 倍；单位或布局错会给 O(1) 偏差。
    """
    h = read_orca_hess(ORCA_WATER)
    assert h.atomic_numbers.tolist() == [8, 1, 1]
    assert h.provider == "orca-hess"
    assert h.gradient is None and h.energy_kj is None
    assert np.array_equal(h.hessian, h.hessian.T)
    # 输入 H 的 y = ±0.7572 Å（ORCA 把坐标平移到质心，y 不变）
    assert h.positions_nm[1, 1] == pytest.approx(0.07572, abs=1e-9)
    masses = np.array([float(m) for m in h.meta["orca_masses"].split(",")])
    assert masses.tolist() == [15.999, 1.008, 1.008]
    ours = harmonic_frequencies(h.hessian, h.positions_nm, masses)
    ref = _orca_frequencies(ORCA_WATER)[6:]
    assert ours == pytest.approx(ref, rel=1e-6)


def test_sf6_roundtrip_multiblock(tmp_path: Path) -> None:
    """SF6 fchk → 按 ORCA 布局写出（21 维、5 块）→ 读回，逐项一致，check_hessian 通过。"""
    ref = read_fchk(SOBTOP / "SF6.fchk")
    symbols = ["S" if z == 16 else "F" for z in ref.atomic_numbers]
    coords_bohr = ref.positions_nm / BOHR_NM
    atom_lines = [
        f"{s} 1.0 {xyz[0]:.17e} {xyz[1]:.17e} {xyz[2]:.17e}"
        for s, xyz in zip(symbols, coords_bohr)
    ]
    path = _write_hess(tmp_path, _blocked(ref.hessian / HA_BOHR2_TO_KJ_NM2), atom_lines, "SF6.hess")
    h = read_orca_hess(path)
    assert h.atomic_numbers.tolist() == ref.atomic_numbers.tolist()
    assert np.allclose(h.positions_nm, ref.positions_nm, rtol=1e-12, atol=0.0)
    assert np.allclose(h.hessian, ref.hessian, rtol=1e-12, atol=1e-6)
    diag = check_hessian(h)
    assert diag.frequencies_cm.size == 15
    assert diag.frequencies_cm.min() > 0.0


def test_d_exponents_and_atomic_number_tokens(tmp_path: Path) -> None:
    hess = ["3", "0 1 2"] + [f"{r} " + " ".join(["1.0D+00"] * 3) for r in range(3)]
    path = _write_hess(tmp_path, hess, ["8 15.999 0.0 0.0 0.0"])
    h = read_orca_hess(path)
    assert h.atomic_numbers.tolist() == [8]
    assert h.hessian[0, 0] == pytest.approx(HA_BOHR2_TO_KJ_NM2)


def test_not_orca_file_raises(tmp_path: Path) -> None:
    hess = ["3", "0 1 2"] + [f"{r} 1.0 0.0 0.0" for r in range(3)]
    path = _write_hess(tmp_path, hess, ["O 15.999 0 0 0"], header="$something_else")
    with pytest.raises(ValueError, match="orca_hessian_file"):
        read_orca_hess(path)


def test_missing_hessian_section_raises(tmp_path: Path) -> None:
    path = tmp_path / "no_hess.hess"
    path.write_text("$orca_hessian_file\n$atoms\n1\nO 15.999 0 0 0\n$end\n")
    with pytest.raises(ValueError, match="\\$hessian"):
        read_orca_hess(path)


def test_missing_atoms_section_raises(tmp_path: Path) -> None:
    hess = ["3", "0 1 2"] + [f"{r} 1.0 0.0 0.0" for r in range(3)]
    path = tmp_path / "no_atoms.hess"
    path.write_text("\n".join(["$orca_hessian_file", "$hessian", *hess, "$end"]) + "\n")
    with pytest.raises(ValueError, match="\\$atoms"):
        read_orca_hess(path)


def test_truncated_block_raises(tmp_path: Path) -> None:
    hess = ["3", "0 1 2", "0 1.0 0.0 0.0", "1 0.0 1.0 0.0"]  # 少一行
    path = _write_hess(tmp_path, hess, ["O 15.999 0 0 0"])
    with pytest.raises(ValueError, match="不足|格式不符"):
        read_orca_hess(path)


def test_missing_column_block_raises(tmp_path: Path) -> None:
    hess = ["6", "0 1 2 3 4"] + [f"{r} " + " ".join(["0.0"] * 5) for r in range(6)]  # 缺第 5 列
    path = _write_hess(tmp_path, hess, ["H 1.008 0 0 0", "H 1.008 1 0 0"])
    with pytest.raises(ValueError, match="不完整"):
        read_orca_hess(path)


def test_atom_count_mismatch_raises(tmp_path: Path) -> None:
    hess = ["3", "0 1 2"] + [f"{r} 1.0 0.0 0.0" for r in range(3)]
    path = _write_hess(tmp_path, hess, ["O 15.999 0 0 0", "H 1.008 1 0 0"])
    with pytest.raises(ValueError, match="不一致"):
        read_orca_hess(path)


def test_atom_line_without_mass_raises(tmp_path: Path) -> None:
    hess = ["3", "0 1 2"] + [f"{r} 1.0 0.0 0.0" for r in range(3)]
    path = _write_hess(tmp_path, hess, ["O 0.0 0.0 0.0"])  # 缺质量列
    with pytest.raises(ValueError, match="列数不符"):
        read_orca_hess(path)
