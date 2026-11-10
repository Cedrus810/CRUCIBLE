"""MVP-2 冒烟 3：MACE-POLAR-1-M 能否经 openmm-ml 1.6 的 modelPath 路径加载、带总电荷求能。
用法: python smoke_polar.py <fixture> <device> [model]"""
import sys, time, contextlib
from pathlib import Path
import numpy as np
import openmm
from openmm import unit
import torch
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from ffcompiler.importers.amber_prmtop import from_amber
from ffcompiler.exporters.openmm import to_openmm_topology

name, device = sys.argv[1], sys.argv[2]
model = sys.argv[3] if len(sys.argv) > 3 else "MACE-POLAR-1-M"
d = REPO / "tests/fixtures/gaff2" / name
ir, pos = from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")

@contextlib.contextmanager
def _load_to(dev):
    orig = torch.load
    def patched(*a, **k):
        obj = orig(*a, **k)
        return obj.to(dev) if isinstance(obj, torch.nn.Module) else obj
    torch.load = patched
    try:
        yield
    finally:
        torch.load = orig

from openmmml import MLPotential
path = f"/home/ruigengji/MLP/mace/{model}.model"
t0 = time.time()
m = torch.load(path, map_location="cpu", weights_only=False)
print(f"model class={type(m).__name__} dtype={next(m.parameters()).dtype} "
      f"elements={len(m.atomic_numbers)} r_max={float(m.r_max):.1f}")
pot = MLPotential("mace", modelPath=path)
with _load_to(device):
    system = pot.createSystem(to_openmm_topology(ir), returnEnergyType="interaction_energy",
                              device=device, charge=ir.total_charge, multiplicity=1)
print(f"load+createSystem {time.time()-t0:.1f}s total_charge={ir.total_charge}")
ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
ctx.setPositions(pos)
e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
t0 = time.time()
for _ in range(5):
    f = ctx.getState(getEnergy=True, getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
print(f"natoms={len(ir.atoms)} E={e:.4f} kJ/mol  energy+forces {(time.time()-t0)/5*1000:.1f} ms  max|F|={np.abs(f).max():.1f}")
