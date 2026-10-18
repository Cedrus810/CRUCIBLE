"""refinement_rounds=0（默认）vs 10-18 G2 存档 IR（改动前代码，同在 181），只读缓存。"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from diagnose_torsion_profiles import _load_fixture
from ffcompiler.ir import io as ir_io
from ffcompiler.torsion.mace import default_model_path
from ffcompiler.torsion.mlp import MLPTorsionProvider

def forbidden(ir):
    raise AssertionError("reference potential must not be called")

g2 = ROOT / "jobs/gate_mvp2_20261018"
for name in sorted(p.name for p in (g2 / "scans").iterdir() if p.is_dir()):
    base, pos, mol = _load_fixture(name)
    ir0, rep = MLPTorsionProvider(forbidden, supported=frozenset(range(1, 119)), scan_cache=g2 / "scans" / name,
                                  model_path=default_model_path(), cache_only=True).parameterize(base, pos, mol)
    stored = ir_io.from_json(g2 / "ir" / f"{name}.json")
    print(f"{name:16s} torsions_identical={ir0.proper_torsions == stored.proper_torsions} "
          f"refinement={rep.refinement} protocol={ir0.provenance['torsion']['fit_protocol']}")
