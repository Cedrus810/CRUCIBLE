"""scripts/torsion_multistart_report.py：只读诊断报告（MVP-7，计划 §3.1）。合成输入，不扫描。"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "torsion_multistart_report.py"


def _load():
    spec = importlib.util.spec_from_file_location("torsion_multistart_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


N = 4
QUARTET = [0, 1, 2, 3]


def _leg(path: Path, energy, hysteresis, quartet=QUARTET):
    np.savez(path, key=np.array("k"), quartet=np.array(quartet),
             phi_target_rad=np.radians(np.arange(N) * 90.0 - 180.0),
             phi_rad=np.zeros(N), energy=np.asarray(energy, float),
             positions=np.zeros((N, 4, 3)), hysteresis=np.asarray(hysteresis, float),
             converged=np.ones(N, bool), max_bond_change_pm=np.array(0.1),
             params=np.array("{}"))


def _fixture(tmp_path: Path):
    scans = tmp_path / "scans"
    scans.mkdir()
    x0_e = [10.0, 30.0, 20.0, 12.0]
    _leg(scans / "1-2.npz", x0_e, [0.0, 50.0, 1.0, 0.0])
    _leg(scans / "1-2.etkdg-1.npz", [10.0, 5.0, 20.0, 12.0], [0.0, 0.5, 0.0, 0.0])
    diag = {"threshold_nm": 0.05, "n_points": N, "branch_count_max": 2, "branch_count_mean": 1.25,
            "n_points_multibranch": 1, "branch_gap_max": 25.0, "envelope_gain_max": 25.0,
            "n_points_improved": 1, "unconverged_points": [],
            "branch_counts": [1, 2, 1, 1], "branch_gap_kJ_mol": [0.0, 25.0, 0.0, 0.0],
            "envelope_delta_kJ_mol": [0.0, -25.0, 0.0, 0.0], "start_labels": ["x0", "etkdg-1"]}
    report = {"molecule_fallback": None,
              "rotors": [{"bond": [1, 2], "rotor_class": 0, "scan_quartet": QUARTET,
                          "status": "mlp", "reason": None, "fit_rmse": 0.3,
                          "max_hysteresis": 0.5, "multistart": diag},
                         {"bond": [2, 3], "rotor_class": 1, "scan_quartet": [1, 2, 3, 4],
                          "status": "fallback", "reason": "not_selected"}],
              "multistart": {"spec": {"n_confs": 2}, "identities": [None, {"label": "etkdg-1"}]}}
    rpath = tmp_path / "torsion_report.json"
    rpath.write_text(json.dumps(report))
    base = {"molecule_fallback": None,
            "rotors": [{"bond": [1, 2], "status": "mlp", "fit_rmse": 2.0, "max_hysteresis": 50.0}]}
    bpath = tmp_path / "base.json"
    bpath.write_text(json.dumps(base))
    return rpath, scans, bpath


def test_build_summary_and_points(tmp_path):
    mod = _load()
    rpath, scans, bpath = _fixture(tmp_path)
    summaries, points = mod.build(rpath, scans, bpath)
    assert len(summaries) == 1                       # 只有带 multistart 诊断的键
    s = summaries[0]
    assert s["leg_max_hysteresis_kJ_mol"] == {"x0": 50.0, "etkdg-1": 0.5}
    assert s["fit_rmse_envelope"] == 0.3 and s["fit_rmse_baseline"] == 2.0
    assert s["missing_legs"] == []
    assert len(points) == N
    # 包络 = x0 + delta；相对能各自减自身最小值
    assert points[1]["x0_rel_kJ_mol"] == 20.0
    assert points[1]["envelope_rel_kJ_mol"] == 0.0   # 包络 [10,5,20,12] 最小 5
    assert points[0]["envelope_rel_kJ_mol"] == 5.0
    assert points[1]["phi_target_deg"] == pytest.approx(-90.0)
    assert points[1]["branch_count"] == 2


def test_main_writes_outputs_without_cache(tmp_path):
    mod = _load()
    rpath, _, _ = _fixture(tmp_path)
    out = tmp_path / "out"
    assert mod.main(["--report", str(rpath), "--outdir", str(out)]) == 0
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "## 键 1-2" in text and "未给" in text
    assert (out / "per_point.csv").read_text().count("\n") == N + 1
    assert json.loads((out / "summary.json").read_text())["bonds"][0]["bond"] == [1, 2]


def test_rejects_cache_from_other_run(tmp_path):
    mod = _load()
    rpath, scans, _ = _fixture(tmp_path)
    _leg(scans / "1-2.npz", [0.0] * N, [0.0] * N, quartet=[9, 1, 2, 3])
    with pytest.raises(ValueError, match="四元组"):
        mod.build(rpath, scans, None)


def test_rejects_single_start_report(tmp_path):
    mod = _load()
    path = tmp_path / "r.json"
    path.write_text(json.dumps({"molecule_fallback": None, "rotors": []}))
    with pytest.raises(ValueError, match="multistart"):
        mod.build(path, None, None)
