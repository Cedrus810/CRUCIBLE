"""T23：gate_mvp1 报告生成函数的必测用例（假指标行驱动）。

passed 必须逐条由 G1.1–G1.7 推出：任何一条不满足都不得为 True（fail-closed）。
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "gate_mvp1", REPO / "scripts" / "gate_mvp1.py"
)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

GOOD_SOBTOP = {"mseminario": "12/12", "m2seminario": "12/12", "drih": "12/12"}
GOOD_PYTEST = {
    "not external": {"accepted": True, "passed": 10},
    "external": {"accepted": True, "passed": 5},
}
GOOD_FOCUS = {"RF1": [("tests/unit/test_assign.py::test_x", "passed")]}


def _row(**overrides):
    base = {
        "molecule": "ethanol", "provider": "gfn2", "method": "drih",
        "freq_rmse_cm": 12.5, "freq_mae_cm": 9.1, "freq_max_cm": 40.2,
        "freq_rmse_lt1000": 5.0, "freq_rmse_ge1000": 20.0,
        "kabsch_rmsd_pm": 30.0, "bond_mae_pm": 2.1, "angle_mae_deg": 0.9,
        "max_class_spread": 0.01,
        "gf_alpha": None, "gf_lam": None, "gf_cond": None,
        "gf_rel_res_fit": None, "gf_rel_res_prior": None,
        "gf_n_at_bound": None, "gf_max_force": None,
        "compare_passed": True,
    }
    base.update(overrides)
    return base


def _rows():
    return [
        _row(),
        _row(method="gaff2", freq_rmse_cm=0.0),
        _row(provider="gfnff", method="global_fit", gf_alpha=1e-3, gf_lam=2.0,
             gf_cond=100.0, gf_rel_res_fit=0.1, gf_rel_res_prior=0.3,
             gf_n_at_bound=0, gf_max_force=5.0),
    ]


def _build(rows=None, failures=None, **overrides):
    kwargs = dict(
        pytest_summary=GOOD_PYTEST, sobtop_match=GOOD_SOBTOP, m2_naming="m2seminario",
        global_fit_self_consistency=4.3e-9, review_focus=GOOD_FOCUS,
    )
    kwargs.update(overrides)
    return gate.build_report(_rows() if rows is None else rows, failures or [], **kwargs)


def test_report_headers_and_acceptance_fields() -> None:
    rows = _rows()
    report_md, metrics_csv, acceptance = _build(rows)
    # 表头固定
    assert report_md.startswith("# Gate G1 评估报告（MVP-1）")
    assert "| method | 指标 | 中位数 | 最大值 |" in report_md
    assert metrics_csv.splitlines()[0] == gate.CSV_HEADER
    assert len(metrics_csv.splitlines()) == 1 + len(rows)
    # acceptance 字段齐全
    for key in ("gate", "passed", "sobtop_match", "m2_naming",
                "global_fit_self_consistency", "cross_engine_all_passed",
                "failures", "pytest_summary", "criteria"):
        assert key in acceptance
    assert acceptance["gate"] == "G1"
    assert acceptance["global_fit_self_consistency"] == 4.3e-9
    assert acceptance["passed"] is True
    assert all(acceptance["criteria"].values())
    json.dumps(acceptance)  # 可序列化
    # G1.6：报告列出 Review Focus 测试名
    assert "tests/unit/test_assign.py::test_x" in report_md


def test_failure_listed_and_passed_false() -> None:
    failures = [{"molecule": "benzene", "provider": "gfn2", "method": "drih",
                 "reason": "boom"}]
    report_md, _, acceptance = _build(failures=failures)
    assert acceptance["passed"] is False
    assert any(f["molecule"] == "benzene" and f["reason"] == "boom"
               for f in acceptance["failures"])
    assert "benzene/gfn2/drih: boom" in report_md


def test_sobtop_mismatch_fails_g1_1() -> None:
    _, _, acceptance = _build(sobtop_match=dict(GOOD_SOBTOP, mseminario="0/12"))
    assert acceptance["criteria"]["G1.1"] is False
    assert acceptance["passed"] is False


def test_m2_named_but_not_matching_fails_g1_2() -> None:
    _, _, acceptance = _build(sobtop_match=dict(GOOD_SOBTOP, m2seminario="11/12"))
    assert acceptance["criteria"]["G1.2"] is False
    assert acceptance["passed"] is False


def test_self_consistency_missing_or_large_fails_g1_3() -> None:
    for value in (None, 2e-5, float("nan")):
        _, _, acceptance = _build(global_fit_self_consistency=value)
        assert acceptance["criteria"]["G1.3"] is False
        assert acceptance["passed"] is False


def test_skipped_compare_is_not_a_pass() -> None:
    rows = _rows()
    rows[0]["compare_passed"] = None  # --skip-compare
    _, _, acceptance = _build(rows)
    assert acceptance["cross_engine_all_passed"] is False
    assert acceptance["criteria"]["G1.4"] is False
    assert acceptance["passed"] is False


def test_nan_metric_fails_g1_5() -> None:
    rows = _rows()
    rows[2]["gf_max_force"] = math.nan
    _, _, acceptance = _build(rows)
    assert acceptance["criteria"]["G1.5"] is False
    assert acceptance["passed"] is False
    assert any("NaN" in f["reason"] for f in acceptance["failures"])


def test_missing_combination_fails_g1_5() -> None:
    expected = [("ethanol", "gfn2", "drih"), ("ethanol", "gfn2", "gaff2"),
                ("ethanol", "gfnff", "global_fit"), ("benzene", "gfn2", "drih")]
    _, _, acceptance = _build(expected=expected)
    assert acceptance["criteria"]["G1.5"] is False
    assert any(f["molecule"] == "benzene" for f in acceptance["failures"])


def test_review_focus_test_not_passed_fails_g1_6() -> None:
    _, _, acceptance = _build(review_focus={"RF1": [("t::a", "missing")]})
    assert acceptance["criteria"]["G1.6"] is False
    assert acceptance["passed"] is False


def test_pytest_not_accepted_fails_g1_7() -> None:
    summary = dict(GOOD_PYTEST, external={"accepted": False, "returncode": 0,
                                          "unexpected_skips": 1})
    _, _, acceptance = _build(pytest_summary=summary)
    assert acceptance["criteria"]["G1.7"] is False
    assert acceptance["passed"] is False


def test_junit_outcomes(tmp_path: Path) -> None:
    xml = tmp_path / "r.xml"
    xml.write_text(
        '<testsuites><testsuite>'
        '<testcase classname="tests.unit.test_a" name="test_ok"/>'
        '<testcase classname="tests.unit.test_a" name="test_bad"><failure/></testcase>'
        '</testsuite></testsuites>'
    )
    outcomes = gate.junit_outcomes([xml])
    assert outcomes["tests/unit/test_a.py::test_ok"] == "passed"
    assert outcomes["tests/unit/test_a.py::test_bad"] == "failed"


def _limitation(evidence):
    return {"molecule": "cyclopropane", "provider": "gfnff", "method": "m2seminario",
            "reason": "非正 k", "category": "method_limitation", "evidence": evidence}


def test_sobtop_confirmed_method_limitation_does_not_fail_g1_5() -> None:
    """决策 0014：有 Sobtop 同 Hessian 对照证据的方法局限照样列出，但不卡 G1.5。"""
    evidence = {"sobtop_confirmed": True,
                "terms": [{"atoms": [0, 1, 2], "ours": -33.59, "sobtop": -33.59}]}
    expected = [("ethanol", "gfn2", "drih"), ("ethanol", "gfn2", "gaff2"),
                ("ethanol", "gfnff", "global_fit"), ("cyclopropane", "gfnff", "m2seminario")]
    report_md, _, acceptance = _build(failures=[_limitation(evidence)], expected=expected)
    assert acceptance["criteria"]["G1.5"] is True
    assert acceptance["passed"] is True
    assert acceptance["method_limitations"][0]["molecule"] == "cyclopropane"
    assert "已知方法局限" in report_md and "cyclopropane/gfnff/m2seminario" in report_md


def test_unconfirmed_method_limitation_fails_g1_5() -> None:
    for evidence in (None, {"sobtop_confirmed": False, "terms": []}):
        _, _, acceptance = _build(failures=[_limitation(evidence)])
        assert acceptance["criteria"]["G1.5"] is False
        assert acceptance["passed"] is False


import pytest  # noqa: E402


@pytest.mark.external
def test_sobtop_confirmation_on_known_limitation(tmp_path: Path) -> None:
    """Gate 现场用 Sobtop 复算同一 Hessian：环丙烷 GFN-FF 的 m2seminario 非正 k 得到确认。"""
    from ffcompiler.bonded.assign import NonPositiveForceConstantError, derive_bonded
    from ffcompiler.hessian.xtb import read_xtb_hessian
    from ffcompiler.importers import from_amber

    fx = REPO / "tests/fixtures/hessian/sobtop_known_limitations/cyclopropane_gfnff"
    d = REPO / "tests/fixtures/gaff2/cyclopropane/cyclopropane"
    ir, positions = from_amber(Path(f"{d}.prmtop"), Path(f"{d}.inpcrd"), Path(f"{d}.sdf"))
    h = read_xtb_hessian(fx / "hessian", fx / "xtbopt.xyz", method="gfnff")
    with pytest.raises(NonPositiveForceConstantError) as info:
        derive_bonded(ir, h, method="m2seminario", input_positions_nm=positions)
    evidence = gate.sobtop_confirm(ir, h, fx / "hessian", info.value, tmp_path)
    assert evidence["sobtop_confirmed"] is True
    assert len(evidence["terms"]) == 3
    assert all(t["sobtop"] < 0 for t in evidence["terms"])


class _FakeAtom:
    mass = 1.0


class _FakeIR:
    atoms = (_FakeAtom(),)


def test_evaluate_matrix_records_every_exception(tmp_path: Path, monkeypatch) -> None:
    """main 的矩阵循环：xtb 失败记到每个方法；评估异常、非正 k 都进 failures，不丢。"""
    from ffcompiler.bonded.assign import NonPositiveForceConstantError

    def fake_hessian(key):
        if key == ("bad", "gfn2"):
            raise RuntimeError("xtb crashed")
        return "H"

    def fake_evaluate(name, provider, method, ir, positions, h, freqs, run_dir, skip):
        if method == "drih":
            raise ValueError("boom")
        if method == "m2seminario":
            raise NonPositiveForceConstantError("m2seminario", [((0, 1, 2), -1.0)])
        return {"molecule": name, "provider": provider, "method": method}

    monkeypatch.setattr(gate, "check_hessian", lambda h, masses: type("D", (), {"frequencies_cm": []})())
    monkeypatch.setattr(gate, "sobtop_confirm",
                        lambda *a, **k: {"sobtop_confirmed": True, "terms": []})
    keys = [("good", "gfn2"), ("bad", "gfn2")]
    systems = {"good": (_FakeIR(), None), "bad": (_FakeIR(), None)}
    rows, failures = gate.evaluate_matrix(keys, systems, fake_hessian, tmp_path, False,
                                          evaluate=fake_evaluate)
    assert {r["method"] for r in rows} == {"gaff2", "mseminario", "global_fit"}
    by = {(f["molecule"], f["method"]): f for f in failures}
    assert by[("good", "drih")]["reason"] == "ValueError: boom"
    assert by[("good", "m2seminario")]["category"] == "method_limitation"
    assert all(("bad", m) in by and "xtb crashed" in by[("bad", m)]["reason"]
               for m in gate.METHODS)
