"""T12：Gaussian fchk 读取器的必测用例。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.hessian.base import check_hessian
from ffcompiler.hessian.gaussian import read_fchk

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "hessian" / "sobtop_examples"


def test_sf6_fchk_geometry_and_hessian() -> None:
    h = read_fchk(FIXTURES / "SF6.fchk")
    assert h.atomic_numbers.tolist() == [16, 9, 9, 9, 9, 9, 9]
    assert h.method == "RB3LYP"
    # S–F 键长（Å）与 Sobtop param.txt 一致，误差 ≤ 1e-6 Å
    distances_a = [
        np.linalg.norm(h.positions_nm[0] - h.positions_nm[i]) * 10.0
        for i in range(1, 7)
    ]
    for d in distances_a:
        assert d == pytest.approx(1.578599, abs=1e-6)
    # check_hessian 通过；最小频率 > 0（八面体非线性：3·7−6 = 15 个频率）
    diag = check_hessian(h)
    assert diag.frequencies_cm.size == 15
    assert diag.frequencies_cm.min() > 0.0


def test_methyl_benzoate_atom_count_matches_mol2() -> None:
    h = read_fchk(FIXTURES / "Methyl_benzoate.fchk")
    mol2_lines = (FIXTURES / "Methyl_benzoate.mol2").read_text().splitlines()
    idx = mol2_lines.index("@<TRIPOS>MOLECULE")
    # MOLECULE 块：分子名行、计数行（原子数 [键数 …]）
    n_mol2 = int(mol2_lines[idx + 2].split()[0])
    assert h.atomic_numbers.size == n_mol2
    check_hessian(h)  # 通过即无异常


def test_lower_triangle_expansion_and_symmetry(tmp_path: Path) -> None:
    # 手写 2 原子 fchk：3N = 6，下三角 21 个值按 (0,0),(1,0),(1,1),… 展开
    n3 = 6
    values = [f"{i + 1}.000000000000D+00" for i in range(n3 * (n3 + 1) // 2)]
    lines = [
        "Title Card Required",
        "Freq      RB3LYP                                                      test",
        "Number of atoms                            I                2",
        "Atomic numbers                             I   N=           2",
        "   1   8",
        "Current cartesian coordinates              R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00",
        "  1.889726124626E+00  0.000000000000E+00  0.000000000000E+00",
        "Cartesian Gradient                         R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00",
        "  0.000000000000E+00  0.000000000000E+00  1.000000000000D-03",
        "Cartesian Force Constants                  R   N=         21",
        *[f"  {v}" for v in values],
        "Total Energy                               R     -7.600000000000000E+01",
    ]
    path = tmp_path / "mini.fchk"
    path.write_text("\n".join(lines) + "\n")

    h = read_fchk(path)
    assert h.atomic_numbers.tolist() == [1, 8]
    assert h.method == "RB3LYP"
    # D 幂写法解析：gradient 第 6 个分量 = 1e-3 Ha/Bohr
    assert h.gradient is not None
    assert h.gradient[5] == pytest.approx(1e-3 * 2625.499639479 / 0.0529177210544)
    # 每个元素位置正确：h[i,j] = h[j,i] = values[下三角行主序编号]
    rows, cols = np.tril_indices(n3)
    expected = np.zeros((n3, n3))
    for k, (i, j) in enumerate(zip(rows, cols)):
        expected[i, j] = expected[j, i] = k + 1
    conversion = 2625.499639479 / 0.0529177210544**2
    assert np.allclose(h.hessian, expected * conversion, rtol=0, atol=1e-9)
    assert np.array_equal(h.hessian, h.hessian.T)
    assert h.energy_kj == pytest.approx(-76.0 * 2625.499639479)


def test_missing_force_constants_raises(tmp_path: Path) -> None:
    lines = [
        "Title Card Required",
        "Freq      RB3LYP                                                      test",
        "Number of atoms                            I                2",
        "Atomic numbers                             I   N=           2",
        "   1   8",
        "Current cartesian coordinates              R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00",
        "  1.889726124626E+00  0.000000000000E+00  0.000000000000E+00",
    ]
    path = tmp_path / "no_hess.fchk"
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="freq"):
        read_fchk(path)


def _mini_fchk_lines(gradient_block: list[str]) -> list[str]:
    """2 原子小 fchk；gradient_block 原样插在坐标与力常数之间。"""
    n3 = 6
    values = [f"{i + 1}.000000000000D+00" for i in range(n3 * (n3 + 1) // 2)]
    return [
        "Title Card Required",
        "Freq      RB3LYP                                                      test",
        "Number of atoms                            I                2",
        "Atomic numbers                             I   N=           2",
        "   1   8",
        "Current cartesian coordinates              R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00",
        "  1.889726124626E+00  0.000000000000E+00  0.000000000000E+00",
        *gradient_block,
        "Cartesian Force Constants                  R   N=         21",
        *[f"  {v}" for v in values],
        "Total Energy                               R     -7.600000000000000E+01",
    ]


def test_absent_gradient_gives_none(tmp_path: Path) -> None:
    path = tmp_path / "no_grad.fchk"
    path.write_text("\n".join(_mini_fchk_lines([])) + "\n")
    assert read_fchk(path).gradient is None


def test_truncated_gradient_raises(tmp_path: Path) -> None:
    # 声明 6 个、只写了 3 个，紧接下一个字段头：不能当成"没有梯度"
    block = [
        "Cartesian Gradient                         R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00  1.000000000000E-01",
    ]
    path = tmp_path / "short_grad.fchk"
    path.write_text("\n".join(_mini_fchk_lines(block)) + "\n")
    with pytest.raises(ValueError, match="Cartesian Gradient"):
        read_fchk(path)


def test_truncated_gradient_at_eof_raises(tmp_path: Path) -> None:
    lines = _mini_fchk_lines([])
    lines += [
        "Cartesian Gradient                         R   N=          6",
        "  0.000000000000E+00  0.000000000000E+00",
    ]
    path = tmp_path / "eof_grad.fchk"
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="Cartesian Gradient"):
        read_fchk(path)


def test_gradient_wrong_length_raises(tmp_path: Path) -> None:
    # 字段本身完整，但个数 ≠ 3N
    block = [
        "Cartesian Gradient                         R   N=          3",
        "  0.000000000000E+00  0.000000000000E+00  1.000000000000E-01",
    ]
    path = tmp_path / "len_grad.fchk"
    path.write_text("\n".join(_mini_fchk_lines(block)) + "\n")
    with pytest.raises(ValueError, match="gradient|梯度"):
        read_fchk(path)
