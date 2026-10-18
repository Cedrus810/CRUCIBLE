"""fallback（hessian/fallback.py，MVP-5）的必测用例：假 runner 注入，
不调用 xtb；覆盖接受/拒绝/不可用/抛错与全失败的报告。"""

from __future__ import annotations

import numpy as np
import pytest

from ffcompiler.hessian.base import CartesianHessian, HessianQualityError
from ffcompiler.hessian.fallback import (
    FallbackError,
    Tier,
    TierUnavailableError,
    run_with_fallback,
)
from ffcompiler.bonded.internal import wilson_b


def _harmonic(k: float, scale: float = 1.0) -> CartesianHessian:
    """双原子谐振 Hessian（k<0 给虚频；|k|=5 → |ν|≈17 cm⁻¹ 只警告，
    |k|=10 → ≈24 cm⁻¹ 超过 −20 报错阈值）。"""
    pos = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    b = wilson_b(pos, [("bond", (0, 1))])[0]
    h = k * np.outer(b, b)
    return CartesianHessian(
        atomic_numbers=np.array([1, 1]), positions_nm=pos, hessian=h,
        gradient=np.zeros(6), energy_kj=None, provider="test",
        method=f"synthetic-{scale}",
    )


def test_first_tier_accepted_immediately() -> None:
    calls = []

    def good(**kwargs):
        calls.append("good")
        return _harmonic(5.0)

    h, report = run_with_fallback([Tier("t1", good), Tier("t2", good)])
    assert report.accepted_tier == "t1"
    assert [a.status for a in report.attempts] == ["accepted"]
    assert calls == ["good"]
    assert h.method.startswith("synthetic")


def test_quality_rejection_escalates_to_next_tier() -> None:
    def bad(**kwargs):
        return _harmonic(-10.0)   # 虚频 ≈ −24 cm⁻¹ → HessianQualityError

    def good(**kwargs):
        return _harmonic(5.0)

    h, report = run_with_fallback([Tier("t1", bad), Tier("t2", good)])
    assert report.accepted_tier == "t2"
    assert [a.status for a in report.attempts] == ["rejected", "accepted"]
    assert "虚频" in report.attempts[0].detail


def test_runner_exception_recorded_as_error() -> None:
    def boom(**kwargs):
        raise RuntimeError("xtb died")

    def good(**kwargs):
        return _harmonic(5.0)

    h, report = run_with_fallback([Tier("t1", boom), Tier("t2", good)])
    assert report.accepted_tier == "t2"
    assert report.attempts[0].status == "error"
    assert "RuntimeError" in report.attempts[0].detail


def test_unavailable_tier_skipped() -> None:
    def unavailable(**kwargs):
        raise TierUnavailableError("no DFT file")

    def good(**kwargs):
        return _harmonic(5.0)

    h, report = run_with_fallback([Tier("dft", unavailable), Tier("gfn2", good)])
    assert report.accepted_tier == "gfn2"
    assert report.attempts[0].status == "unavailable"


def test_warning_does_not_escalate() -> None:
    # 浅虚频（|ν|≈17 cm⁻¹，−20…0 区间）只警告：接受并记录 warning
    def shallow(**kwargs):
        return _harmonic(-5.0)

    h, report = run_with_fallback([Tier("t1", shallow)])
    assert report.accepted_tier == "t1"
    assert report.attempts[0].status == "accepted"
    assert any("虚频" in w for w in report.attempts[0].warnings)


def test_all_fail_raises_with_report() -> None:
    def bad(**kwargs):
        return _harmonic(-10.0)

    with pytest.raises(FallbackError) as info:
        run_with_fallback([Tier("gfnff", bad), Tier("gfn2", bad)])
    report = info.value.report
    assert report.accepted_tier is None
    assert [a.status for a in report.attempts] == ["rejected", "rejected"]
    assert report.to_dict()["attempts"][0]["tier"] == "gfnff"


def test_default_ladder_dft_tier_requires_file(tmp_path) -> None:
    from ffcompiler.hessian.fallback import default_ladder

    pos = np.zeros((2, 3))
    tiers = default_ladder(np.array([1, 1]), pos, charge=0,
                           workdir=tmp_path, dft_file=None)
    assert [t.name for t in tiers] == ["gfnff", "gfn2", "dft"]
    with pytest.raises(TierUnavailableError, match="DFT"):
        tiers[2].runner()


def test_default_ladder_dft_file_must_match_molecule() -> None:
    from pathlib import Path

    from ffcompiler.hessian.fallback import default_ladder

    sf6 = Path(__file__).resolve().parents[1] / "fixtures" / "hessian" / "sobtop_examples" / "SF6.fchk"
    tiers = default_ladder(np.array([1, 1]), np.zeros((2, 3)), charge=0,
                           workdir=Path("unused"), dft_file=sf6)
    with pytest.raises(ValueError, match="原子序数序列"):
        tiers[2].runner()
