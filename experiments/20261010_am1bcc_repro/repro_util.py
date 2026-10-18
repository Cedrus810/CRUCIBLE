"""公共部分（摘自 repro.py）：按 make_fixtures.py 的做法在固件 SDF 上重跑 antechamber -c bcc，与固件 mol2 电荷逐原子比较；
顺带算 NAGL（openff-gnn-am1bcc-1.0.0）与固件的差（只看量级）。输出写本目录。"""
import json, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "src"))

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

