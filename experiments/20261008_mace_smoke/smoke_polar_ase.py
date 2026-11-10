"""MVP-2 冒烟 4：MACE-POLAR-1-M 经 mace 自带 ASE 计算器（不经 openmm-ml）带总电荷求能。
用法: python smoke_polar_ase.py <fixture> <device> [dtype float32|float64]"""
import sys, time
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from ffcompiler.importers.amber_prmtop import from_amber
from ase import Atoms
from mace.calculators import MACECalculator

name, device = sys.argv[1], sys.argv[2]
dtype = sys.argv[3] if len(sys.argv) > 3 else "float64"
d = REPO / "tests/fixtures/gaff2" / name
ir, pos = from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")
calc = MACECalculator(model_paths="/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model", device=device,
                      default_dtype=dtype, model_type="PolarMACE")
atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos * 10.0)
atoms.calc = calc
for q in (ir.total_charge, 0) if ir.total_charge != 0 else (0,):
    atoms.info["charge"] = q; atoms.info["spin"] = 1
    calc.reset()
    e = atoms.get_potential_energy() * 96.48533212
    t0 = time.time()
    for _ in range(5):
        calc.reset(); atoms.get_forces()
    dt = (time.time() - t0) / 5 * 1000
    ch = calc.results.get("charges")
    chs = "" if ch is None else f" sum(q)={np.sum(ch):+.4f} q[:4]={np.round(ch[:4], 3).tolist()}"
    print(f"{name} charge_in={q:+d} E={e:.3f} kJ/mol  energy+forces {dt:.1f} ms{chs}")
print("result keys:", sorted(calc.results.keys()))
