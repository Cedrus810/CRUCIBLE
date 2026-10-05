import pytest
from ffcompiler.validation.gromacs import parse_xvg, find_gromacs


@pytest.mark.parametrize("text", ["", "0 1\n1 2", "0 nan", "0 inf", "0 1 2", "1 2"])
def test_damaged_energy_output_rejected(text):
    with pytest.raises(ValueError):
        parse_xvg(text, ["total"])


def test_precise_energy_parsing():
    assert parse_xvg('# header\n@ s0 legend "Potential"\n0 -8.448573983764\n', ["total"])["total"] == -8.448573983764


def test_explicit_bad_executable_is_error(monkeypatch):
    monkeypatch.setenv("FFC_GMX", "/does/not/exist/gmx")
    with pytest.raises(ValueError, match="not executable"):
        find_gromacs()


def test_empty_report_and_engine_errors_cannot_pass():
    from ffcompiler.validation.engines import CrossEngineReport
    report = CrossEngineReport()
    assert not report.passed
    report.rows.append({"passed": True})
    report.errors["gmx"] = "failed to execute"
    assert not report.passed


def test_gmx_tolerance_double_build_stays_strict():
    """决策 0013：double 版 GROMACS 仍按原经文件容差 1e-4 + 1e-7|E|。"""
    from ffcompiler.validation.engines import gmx_energy_tolerance
    assert gmx_energy_tolerance("bond", 100.0, 100.0, "double") == pytest.approx(1e-4 + 1e-5)


def test_gmx_tolerance_mixed_build_covers_fp32_rounding():
    """决策 0013：mixed 版 1e-4 + 1e-5|E|；G0 实测最大 fp32 偏差 2.7e-6 相对必须通过。"""
    from ffcompiler.validation.engines import gmx_energy_tolerance
    tol = gmx_energy_tolerance("bond", 61.89, 61.89, "mixed")
    assert tol == pytest.approx(1e-4 + 1e-5 * 61.89)
    assert 1.66e-4 <= tol  # 苯 cfg2 bond 实测偏差


def test_gmx_tolerance_mixed_still_catches_conversion_errors():
    """真实转换错误 ≥1e-3 相对（决策 0004.1）必须超出 mixed 容差。"""
    from ffcompiler.validation.engines import gmx_energy_tolerance
    for e in (0.5, 10.0, 100.0, 1000.0):
        assert 1e-3 * e > gmx_energy_tolerance("bond", e, e, "mixed")


def test_gmx_tolerance_mixed_total_scales_with_largest_component():
    """total 行与 0004.1 一致：相对项按最大分量量级（成键/非键抵消时 |total| 很小）。"""
    from ffcompiler.validation.engines import gmx_energy_tolerance
    assert gmx_energy_tolerance("total", 1.0, 500.0, "mixed") == pytest.approx(1e-4 + 1e-5 * 500.0)


def test_gmx_tolerance_unknown_precision_is_strict():
    from ffcompiler.validation.engines import gmx_energy_tolerance
    assert gmx_energy_tolerance("bond", 100.0, 100.0, None) == pytest.approx(1e-4 + 1e-5)
