"""T37：gate_mvp3 报告生成函数（假结果驱动；判据 fail-closed，MVP-3 计划 §6）。"""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("gate_mvp3", REPO / "scripts" / "gate_mvp3.py")
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

GOOD_PYTEST = {"not external": {"accepted": True}, "external": {"accepted": True}}


def _build(**kw):
    args = dict(
        full_matrix=True,
        g31={"ethanol": 0.0, "acetate": 0.0},
        g32={("am1bcc", "ethanol"): 1e-16, ("nagl", "ethanol"): 0.0},
        compare={("am1bcc", "ethanol"): True, ("nagl", "ethanol"): True},
        sander_notes={},
        g34={"ethanol": {"identical": True, "nagl_statuses": ["mlp", "mlp"]}},
        nagl_vs_am1bcc={"ethanol": {"max_abs": 0.04, "rms": 0.018}},
        profiles=[{"molecule": "ethanol", "bond": "1-2", "charges": "nagl", "profile_rmse": 0.3}],
        pytest_summary=GOOD_PYTEST,
        failures=[],
        expected_g31=["ethanol", "acetate"],
        expected_compare=[("am1bcc", "ethanol"), ("nagl", "ethanol")],
        expected_g34={"ethanol": 2},
    )
    args.update(kw)
    return gate.build_report(**args)


def test_report_header_and_acceptance_fields() -> None:
    report_md, acc = _build()
    assert report_md.startswith(gate.REPORT_HEADER)
    assert acc["gate"] == "G3" and acc["passed"] is True
    assert set(acc["criteria"]) == {"G3.1", "G3.2", "G3.3", "G3.4", "G3.5"}
    for key in ("g31", "g32", "compare", "g34", "nagl_vs_am1bcc", "sander_known_limitation",
                "failures", "pytest_summary"):
        assert key in acc


def test_each_criterion_is_fail_closed() -> None:
    assert _build(g31={"ethanol": 0.0})[1]["criteria"]["G3.1"] is False            # 缺分子
    assert _build(g31={"ethanol": 0.0, "acetate": 2e-10})[1]["criteria"]["G3.1"] is False
    assert _build(g32={("am1bcc", "ethanol"): 1e-16, ("nagl", "ethanol"): 3e-10})[1]["criteria"]["G3.2"] is False
    assert _build(g32={("am1bcc", "ethanol"): 0.0})[1]["criteria"]["G3.2"] is False  # 缺组合
    assert _build(compare={("am1bcc", "ethanol"): True, ("nagl", "ethanol"): None})[1]["criteria"]["G3.3"] is False
    assert _build(g34={"ethanol": {"identical": False, "nagl_statuses": ["mlp", "mlp"]}})[1]["criteria"]["G3.4"] is False
    assert _build(g34={"ethanol": {"identical": True, "nagl_statuses": ["mlp"]}})[1]["criteria"]["G3.4"] is False
    assert _build(g34={})[1]["criteria"]["G3.4"] is False
    bad = {"not external": {"accepted": True}, "external": {"accepted": False}}
    assert _build(pytest_summary=bad)[1]["criteria"]["G3.5"] is False


def test_failures_fail_gate_and_sander_notes_only_reported() -> None:
    _, acc = _build(failures=[{"stage": "g31", "molecule": "nma", "reason": "boom"}])
    assert acc["passed"] is False
    report_md, acc = _build(sander_notes={"nagl/nma": ["sander-ours dihed: delta 3e-4 > tol 2e-4"]})
    assert acc["passed"] is True and "决策 0017" in report_md


def test_subset_run_never_passes() -> None:
    # 决策 0019：各判据全部满足，但只跑了子集 → passed=false、scope=subset
    report_md, *_, acc = _build(full_matrix=False)
    assert acc["scope"] == "subset"
    assert acc["passed"] is False
    assert all(acc["criteria"].values())
    assert "scope=subset" in report_md


def test_full_run_records_scope() -> None:
    acc = _build()[-1]
    assert acc["scope"] == "full" and acc["passed"] is True


def test_scan_cache_is_required_and_must_exist(tmp_path) -> None:
    """--scan-cache 无默认值且必须是目录：G3 必须显式指向本次 G2 的 scans/（2026-10-13 用户定）。"""
    import subprocess
    import sys

    script = str(REPO / "scripts" / "gate_mvp3.py")
    missing = subprocess.run([sys.executable, script], capture_output=True, text=True, timeout=60)
    assert missing.returncode == 2 and "--scan-cache" in missing.stderr
    wrong = subprocess.run([sys.executable, script, "--scan-cache", str(tmp_path / "nope")],
                           capture_output=True, text=True, timeout=60)
    assert wrong.returncode == 2 and "不是目录" in wrong.stderr
