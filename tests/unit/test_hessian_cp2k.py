"""CP2K .out Hessian 读取器单元测试（决策 0023）。

夹具取自 Sobtop 官方示例 FeCl2/exact/freq.out（CP2K 前向差分 VIB Hessian，
75 原子 5×5×1 周期晶胞）的裁剪版：GLOBAL 行 + QUICKSTEP 坐标段 + VIB Hessian 段。

已知数据质量注记：该文件打印的是 raw 前向差分 Jacobian——与 CP2K 自己报告的
低频/本征向量不自洽（Rayleigh 商验证，见决策 0023）。Reader 的职责边界是忠实
解析打印值（本测试逐位核对），数值质量记录在 meta 里由上层决定。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.hessian.base import CartesianHessian, symmetrized
from ffcompiler.hessian.cp2k import read_cp2k_hess
from ffcompiler.units import HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "hessian" / "cp2k" / (
    "FeCl2_5x5x1_freq.out"
)


@pytest.fixture(scope="module")
def h() -> CartesianHessian:
    return read_cp2k_hess(FIXTURE)


def test_shape_and_elements(h: CartesianHessian):
    assert len(h.atomic_numbers) == 75
    assert h.hessian.shape == (225, 225)
    # QUICKSTEP 坐标表顺序：Cl Fe Cl Cl Fe …
    assert list(h.atomic_numbers[:5]) == [17, 26, 17, 17, 26]


def test_diagonal_value_bit_exact(h: CartesianHessian):
    # 文件打印的 H[1,1] = 1.136309 E_h/Bohr²，对角不受对称化影响
    assert h.hessian[0, 0] == pytest.approx(1.136309 * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2)


def test_offdiagonal_symmetrized(h: CartesianHessian):
    # 文件第 5 行（Fe x）第 14 列 = -1.548430；第 14 行第 5 列 = -0.197095；
    # 对称化后 = 两者平均
    raw_5_14 = -1.548430
    raw_14_5 = -0.197095
    expect = 0.5 * (raw_5_14 + raw_14_5) * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2
    assert h.hessian[4, 13] == pytest.approx(expect)
    assert h.hessian[13, 4] == pytest.approx(expect)


def test_symmetrized_on_read(h: CartesianHessian):
    assert np.array_equal(h.hessian, h.hessian.T)
    # raw 不对称度记入 meta（E_h/Bohr² 单位）
    assert float(h.meta["raw_max_asymmetry_eh_bohr2"]) == pytest.approx(1.351, abs=1e-3)
    assert "symmetrized" in h.meta["note"]


def test_positions_from_quickstep_table(h: CartesianHessian):
    # 坐标表第 1 个 Cl：1.699928 0.981454 12.455510 Å
    assert h.positions_nm[0] == pytest.approx(
        np.array([1.699928, 0.981454, 12.455510]) * 0.1
    )


def test_positions_override(tmp_path: Path):
    fake = np.zeros((75, 3))
    fake[0] = [1.0, 2.0, 3.0]
    h = read_cp2k_hess(FIXTURE, positions_nm=fake)
    assert h.positions_nm[0] == pytest.approx([1.0, 2.0, 3.0])


def test_provider_metadata(h: CartesianHessian):
    assert h.provider == "cp2k-out"
    assert h.source_file.endswith("FeCl2_5x5x1_freq.out")
    assert len(h.source_sha256) == 64


def test_missing_hessian_section(tmp_path: Path):
    bad = tmp_path / "no_hess.out"
    bad.write_text("MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n")
    with pytest.raises(ValueError, match="VIBRATIONAL_ANALYSIS"):
        read_cp2k_hess(bad)


def test_incomplete_matrix(tmp_path: Path):
    lines = FIXTURE.read_text().splitlines()
    # 截掉最后一页的一部分
    bad = tmp_path / "trunc.out"
    bad.write_text("\n".join(lines[: -100]) + "\n")
    with pytest.raises(ValueError, match="225 行|不完整"):
        read_cp2k_hess(bad)


def test_symmetrized_helper_roundtrip(h: CartesianHessian):
    s = symmetrized(h)
    assert np.array_equal(s.hessian, h.hessian)  # 读入时已对称化
