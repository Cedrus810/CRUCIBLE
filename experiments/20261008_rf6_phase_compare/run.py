"""RF6 冒烟：任意相位二面角经真实 sander / gmx（mixed）的跨引擎比较。"""
import sys
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO))
from tests.unit.test_export_arbitrary_phase import _phased_ethanol
from ffcompiler.validation.engines import compare
ir, pos = _phased_ethanol()
print(sorted({round(t.phase, 3) for t in ir.proper_torsions if t.source == "mlp"}))
rep = compare(ir, pos, Path(__file__).parent / "work")
print("passed:", rep.passed)
print(rep.to_markdown()[-1500:])
