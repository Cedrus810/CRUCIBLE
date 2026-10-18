"""Gate G2（jobs/gate_mvp2_20261008_3）里 nma 的 compare=False：用存下的 ir/nma.json 重跑，打印逐腿结果。"""
import json
import sys
from pathlib import Path

from ffcompiler.importers import from_amber
from ffcompiler.ir import from_json
from ffcompiler.validation.engines import compare

REPO = Path(__file__).resolve().parents[2]
out = Path(__file__).resolve().parent
mdir = REPO / "tests/fixtures/gaff2/nma"
_, pos = from_amber(*(mdir / f"nma.{e}" for e in ("prmtop", "inpcrd", "sdf")))
which = sys.argv[1] if len(sys.argv) > 1 else "mlp"
if which == "mlp":
    ir = from_json(REPO / "jobs/gate_mvp2_20261008_3/ir/nma.json")
else:
    ir, _ = from_amber(*(mdir / f"nma.{e}" for e in ("prmtop", "inpcrd", "sdf")))
rep = compare(ir, pos, out / f"work_{which}")
print(which, "passed", rep.passed)
print(rep.to_markdown())
(out / f"report_{which}.json").write_text(json.dumps(rep.to_json(), indent=1, default=str))
