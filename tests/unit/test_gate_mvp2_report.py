"""T30：gate_mvp2 报告生成函数（假指标行驱动；passed 由 G2.1–G2.5 逐条推出，fail-closed）。"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("gate_mvp2", REPO / "scripts" / "gate_mvp2.py")
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

GOOD_PYTEST = {"not external": {"accepted": True}, "external": {"accepted": True}}
GOOD_FOCUS = {"RF1": [("tests/unit/test_torsion_rotors.py::test_rotor_counts[benzene]", "passed")]}
GOOD_SELF = {"G2.1 fourier": "passed", "G2.1 mm": "passed"}


def _row(**kw):
    base = {"molecule": "ethanol", "bond": "1-2", "rotor_class": 0, "status": "mlp", "reason": "",
            "fit_rmse": 0.12, "profile_rmse_fit": 0.3, "profile_rmse_gaff2": 2.5,
            "max_hysteresis": 0.01, "max_bond_change_pm": 0.5}
    base.update(kw)
    return base


def _build(rows=None, failures=None, compare=None, **kw):
    args = dict(full_matrix=True, pytest_summary=GOOD_PYTEST, review_focus=GOOD_FOCUS, self_consistency=GOOD_SELF,
                expected=[("ethanol", "1-2"), ("ethanol", "0-1")],
                compare={"ethanol": True} if compare is None else compare)
    args.update(kw)
    rows = [_row(), _row(bond="0-1", status="fallback", reason="mlp_error", fit_rmse=None,
                         profile_rmse_fit=None)] if rows is None else rows
    return gate.build_report(rows, failures or [], **args)


def test_report_header_and_acceptance_fields() -> None:
    report_md, csv, acc = _build()
    assert report_md.startswith(gate.REPORT_HEADER)
    assert csv.splitlines()[0] == gate.CSV_HEADER
    assert acc["gate"] == "G2" and acc["passed"] is True
    assert set(acc["criteria"]) == {"G2.1", "G2.2", "G2.3", "G2.4", "G2.5"}
    for key in ("failures", "fallbacks", "compare", "pytest_summary", "rotor_counts"):
        assert key in acc
    assert acc["fallbacks"] == [{"molecule": "ethanol", "bond": "0-1", "reason": "mlp_error"}]


def test_exception_goes_to_failures_and_fails_gate() -> None:
    _, _, acc = _build(failures=[{"molecule": "nma", "bond": "1-3", "reason": "worker crashed"}])
    assert acc["passed"] is False and acc["criteria"]["G2.2"] is False
    assert any(f["reason"] == "worker crashed" for f in acc["failures"])


def test_missing_rotor_and_nan_are_failures() -> None:
    _, _, acc = _build(rows=[_row(fit_rmse=math.nan)])
    reasons = " ".join(f["reason"] for f in acc["failures"])
    assert "NaN" in reasons and "0-1" in " ".join(f"{f['bond']}" for f in acc["failures"])
    assert acc["passed"] is False


def test_each_criterion_is_fail_closed() -> None:
    assert _build(compare={"ethanol": False})[2]["criteria"]["G2.3"] is False
    assert _build(compare={})[2]["criteria"]["G2.3"] is False
    assert _build(self_consistency={"G2.1 fourier": "passed", "G2.1 mm": "missing"})[2]["criteria"]["G2.1"] is False
    assert _build(review_focus={"RF1": [("x", "failed")]})[2]["criteria"]["G2.4"] is False
    assert _build(review_focus={})[2]["criteria"]["G2.4"] is False
    bad = {"not external": {"accepted": True}, "external": {"accepted": False}}
    assert _build(pytest_summary=bad)[2]["criteria"]["G2.5"] is False


# 决策 0017：G2.3 以 OpenMM 为准，sander 腿只报告
def _xreport(rows=(), force_rows=(), errors=None):
    def row(leg, passed, category="dihed"):
        return {"leg": leg, "reference": "numpy-ir", "category": category, "E_leg": 1.0,
                "E_ref": 1.0, "delta": 0.0 if passed else 3e-4, "tol": 2e-4, "passed": passed}
    return {"rows": [row(*r) for r in rows],
            "force_rows": [{"leg": leg, "reference": "openmm-ir", "max_df": 0.0, "tol": 1e-6,
                            "passed": ok} for leg, ok in force_rows],
            "skipped": {}, "errors": errors or {}, "metadata": {}}


def test_sander_only_failure_passes_and_is_noted() -> None:
    ok, notes = gate.judge_compare(_xreport(rows=[("openmm-prmtop", True), ("sander-ours", False)]))
    assert ok is True
    assert len(notes) == 1 and "sander-ours" in notes[0] and "dihed" in notes[0]


def test_non_sander_failures_still_fail_g23() -> None:
    assert gate.judge_compare(_xreport(rows=[("openmm-prmtop", False), ("sander-ours", True)]))[0] is False
    assert gate.judge_compare(_xreport(rows=[("openmm-top", True)], force_rows=[("openmm-xml", False)]))[0] is False
    assert gate.judge_compare(_xreport(rows=[("openmm-top", True)], errors={"gmx": "boom"}))[0] is False
    assert gate.judge_compare(_xreport(rows=[("sander-ours", True)]))[0] is False  # 没有非 sander 腿


def test_sander_notes_reported_not_counted() -> None:
    report_md, _, acc = _build(sander_notes={"nma": ["sander-ours dihed: delta 2.3e-4 > tol 2e-4"]})
    assert acc["passed"] is True
    assert acc["sander_known_limitation"] == {"nma": ["sander-ours dihed: delta 2.3e-4 > tol 2e-4"]}
    assert "决策 0017" in report_md and "sander-ours dihed" in report_md


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


def test_seed_scan_cache_copies_without_touching_old(tmp_path) -> None:
    old = tmp_path / "old" / "scans"
    (old / "ethanol").mkdir(parents=True)
    (old / "dmso").mkdir()
    (old / "ethanol" / "1-2.npz").write_bytes(b"scan-a")
    (old / "ethanol" / "1-2.tmp.npz").write_bytes(b"half")
    (old / "dmso" / "0-1.npz").write_bytes(b"scan-b")
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in old.rglob("*") if p.is_file()}
    new = tmp_path / "new" / "scans"
    copied = gate.seed_scan_cache(old, new, ["ethanol"])
    assert copied == ["ethanol/1-2.npz"]
    assert (new / "ethanol" / "1-2.npz").read_bytes() == b"scan-a"
    assert not (new / "ethanol" / "1-2.tmp.npz").exists()
    assert not (new / "dmso").exists()
    after = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in old.rglob("*") if p.is_file()}
    assert after == before


def test_resume_without_scans_dir_errors_before_running(tmp_path) -> None:
    """--resume 路径写错（没有 scans/）在解析参数时报错退出，不静默全重扫（2026-10-13 用户定）。"""
    import subprocess
    import sys

    proc = subprocess.run([sys.executable, str(REPO / "scripts" / "gate_mvp2.py"), "--resume", str(tmp_path / "nope")],
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 2
    assert "没有 scans/" in proc.stderr
