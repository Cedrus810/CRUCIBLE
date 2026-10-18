"""正式接口 refinement_rounds=3 vs B3 机制实验第 3 轮（jobs/torsion_b3_sequential_20261017），只读缓存。"""
import csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from diagnose_torsion_profiles import _load_fixture
from ffcompiler.ir import io as ir_io
from ffcompiler.torsion.mace import default_model_path
from ffcompiler.torsion.mlp import MLPTorsionProvider

def forbidden(ir):
    raise AssertionError("reference potential must not be called")

g2 = ROOT / "jobs/gate_mvp2_20261013"; b3 = ROOT / "jobs/torsion_b3_sequential_20261017"
ref_rows = {(r["molecule"], r["bond"]): r for r in csv.DictReader(open(b3 / "per_bond.csv")) if r["round"] == "3"}
worst = {"dk": 0.0, "dphase": 0.0, "dbest": 0.0, "dmin": 0.0, "dfit": 0.0}
for name in sorted(p.stem.replace("_round3", "") for p in (b3 / "ir").glob("*_round3.json")):
    base, pos, mol = _load_fixture(name)
    ir3, rep = MLPTorsionProvider(forbidden, supported=frozenset(range(1, 119)), scan_cache=g2 / "scans" / name,
                                  model_path=default_model_path(), cache_only=True,
                                  refinement_rounds=3).parameterize(base, pos, mol)
    ref = ir_io.from_json(b3 / "ir" / f"{name}_round3.json")
    layout = [(t.i, t.j, t.k, t.l, t.periodicity, t.source) for t in ir3.proper_torsions]
    same_layout = layout == [(t.i, t.j, t.k, t.l, t.periodicity, t.source) for t in ref.proper_torsions]
    dk = max(abs(a.force_constant - b.force_constant) for a, b in zip(ir3.proper_torsions, ref.proper_torsions))
    dph = max(abs(np.angle(np.exp(1j * (a.phase - b.phase)))) for a, b in zip(ir3.proper_torsions, ref.proper_torsions))
    worst["dk"] = max(worst["dk"], dk); worst["dphase"] = max(worst["dphase"], dph)
    last = rep.refinement["history"][-1]["bonds"] if rep.refinement else []
    for b in last:
        r = ref_rows[(name, f"{b['bond'][0]}-{b['bond'][1]}")]
        worst["dbest"] = max(worst["dbest"], abs(b["profile_rmse_best_shift"] - float(r["bestshift_rmse"])))
        worst["dmin"] = max(worst["dmin"], abs(b["profile_rmse_min_shift"] - float(r["minshift_rmse"])))
        worst["dfit"] = max(worst["dfit"], abs(b["fit_rmse_at_update"] - float(r["fit_rmse_at_update"])))
    print(f"{name:16s} same_layout={same_layout} max|dk|={dk:.2e} max|dphase|={dph:.2e} "
          f"rounds={rep.refinement['completed_rounds'] if rep.refinement else None}")
print("worst:", {k: f"{v:.2e}" for k, v in worst.items()})
