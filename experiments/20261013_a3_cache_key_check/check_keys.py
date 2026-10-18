"""A3 前检查：用现行代码重算 G2 缓存键，与 npz 里存的键逐字比较（不扫描）。
用法：PYTHONPATH=src:scripts python experiments/20261013_a3_cache_key_check/check_keys.py [SCANS_DIR]"""
import sys
from pathlib import Path

import numpy as np

from gate_mvp2 import NAMES, _load
from ffcompiler.torsion.mace import MODEL_NAME, default_model_path, molecule_info
from ffcompiler.torsion.mlp import SCAN_DEFAULTS, model_file_identity, scan_cache_key
from ffcompiler.torsion.rotors import find_rotors

root = Path(sys.argv[1] if len(sys.argv) > 1 else "jobs/gate_mvp2_20261008_3/scans")
mf = model_file_identity(default_model_path())
match = mismatch = missing = 0
for name in NAMES:
    ir, pos, mol = _load(name)
    info, numbers, bonds = molecule_info(ir), [a.atomic_number for a in ir.atoms], [(b.i, b.j) for b in ir.bonds]
    for r in find_rotors(ir, mol, pos):
        key = scan_cache_key(MODEL_NAME, r.scan_quartet, info, dict(SCAN_DEFAULTS), model_file=mf,
                             numbers=numbers, bonds=bonds)
        p = root / name / f"{r.bond[0]}-{r.bond[1]}.npz"
        if not p.is_file():
            missing += 1; print("missing ", name, r.bond); continue
        with np.load(p) as d:
            ok = str(d["key"]) == key
        match += ok; mismatch += not ok
        print("match   " if ok else "MISMATCH", name, r.bond)
print(f"{root}: {match} match / {mismatch} mismatch / {missing} missing")
