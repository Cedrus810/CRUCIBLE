"""冒烟：NAGL openff-gnn-am1bcc-1.0.0 与固件 AM1-BCC 电荷的差（只看量级）、总电荷、单分子耗时。"""
import json, sys, time
from pathlib import Path
import numpy as np
from repro_util import FIX, mols, mol2_charges
from openff.toolkit import Molecule
from openff.nagl import GNNModel
import openff.nagl_models, os
path = Path(openff.nagl_models.__file__).parent / "models" / "am1bcc" / "openff-gnn-am1bcc-1.0.0.pt"
model = GNNModel.load(str(path), eval_mode=True)
rows = []
for name, info in mols.items():
    if not isinstance(info, dict) or "charge" not in info: continue
    m = Molecule.from_file(str(FIX / name / f"{name}.sdf"), allow_undefined_stereo=True)
    t0 = time.time()
    try:
        q = np.asarray(model.compute_property(m, check_domains=True), dtype=float)
    except Exception as exc:
        rows.append({"name": name, "error": repr(exc)[:300]}); print(rows[-1]); continue
    ref = mol2_charges(FIX / name / f"{name}.mol2")
    elem_ok = [a.atomic_number for a in m.atoms]
    rows.append({"name": name, "wall_s": round(time.time() - t0, 3),
                 "max_abs_diff_vs_fixture": float(np.max(np.abs(q - ref))),
                 "rms_diff_vs_fixture": float(np.sqrt(np.mean((q - ref) ** 2))),
                 "sum_minus_Q": float(q.sum() - info["charge"])})
    print(json.dumps(rows[-1]), flush=True)
(Path(__file__).parent / "nagl_results.json").write_text(json.dumps(rows, indent=2) + "\n")
