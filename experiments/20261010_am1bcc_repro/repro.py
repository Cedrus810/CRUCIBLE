"""冒烟：按 make_fixtures.py 的做法在固件 SDF 上重跑 antechamber -c bcc，与固件 mol2 电荷逐原子比较；
顺带算 NAGL（openff-gnn-am1bcc-1.0.0）与固件的差（只看量级）。输出写本目录。"""
import json, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "src"))
from make_fixtures import _equalize_mol2_charges  # 与固件完全相同的均摊规则

FIX = REPO / "tests" / "fixtures" / "gaff2"
manifest = json.loads((FIX / "manifest.json").read_text())
mols = manifest.get("molecules", manifest)
here = Path(__file__).resolve().parent

def mol2_charges(p):
    out, on = [], False
    for line in p.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"): on = True; continue
        if line.startswith("@<TRIPOS>"): on = False; continue
        if on and line.strip(): out.append(float(line.split()[8]))
    return np.array(out)

rows = []
for name, info in mols.items():
    if not isinstance(info, dict) or "charge" not in info: continue
    d = here / "runs" / name; d.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIX / name / f"{name}.sdf", d / f"{name}.sdf")
    t0 = time.time()
    subprocess.run(["antechamber", "-i", f"{name}.sdf", "-fi", "sdf", "-o", f"{name}.mol2", "-fo", "mol2",
                    "-c", "bcc", "-at", "gaff2", "-nc", str(info["charge"]), "-pf", "y"],
                   cwd=d, check=True, capture_output=True)
    raw = mol2_charges(d / f"{name}.mol2")
    _equalize_mol2_charges(d / f"{name}.mol2", info["charge"], name)
    ours, ref = mol2_charges(d / f"{name}.mol2"), mol2_charges(FIX / name / f"{name}.mol2")
    row = {"name": name, "n": len(ref), "wall_s": round(time.time() - t0, 2),
           "max_abs_diff_equalized": float(np.max(np.abs(ours - ref))),
           "max_abs_diff_raw4dp": float(np.max(np.abs(raw - ref)))}
    try:
        from openff.toolkit import Molecule
        from openff.nagl import GNNModel
        from openff.nagl_models import load_nagl_model_path
        m = Molecule.from_file(str(FIX / name / f"{name}.sdf"), allow_undefined_stereo=True)
        q = GNNModel.load(load_nagl_model_path("openff-gnn-am1bcc-1.0.0.pt")).compute_property(m, check_domains=True)
        q = np.asarray(q, dtype=float)
        row["nagl_max_abs_diff"] = float(np.max(np.abs(q - ref)))
        row["nagl_sum_minus_Q"] = float(q.sum() - info["charge"])
    except Exception as exc:
        row["nagl_error"] = repr(exc)[:200]
    rows.append(row); print(json.dumps(row, ensure_ascii=False), flush=True)
(here / "results.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n")
