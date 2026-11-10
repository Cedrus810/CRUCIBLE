"""MVP-2 冒烟 2：10° 顺序约束松弛扫描，检查键是否保持；比较 CPU / CUDA。
用法: python smoke2.py <fixture> <device cpu|cuda> <n_steps>"""
import sys, time, math, contextlib
from pathlib import Path
import numpy as np
import openmm
from openmm import unit
import torch
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from ffcompiler.importers.amber_prmtop import from_amber
from ffcompiler.exporters.openmm import to_openmm_topology
from ffcompiler.ir.energy import dihedral_angle

name, device, nsteps = sys.argv[1], sys.argv[2], int(sys.argv[3])
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
pot = MLPotential("mace", modelPath="/home/ruigengji/MLP/mace/MACE-OFF24_medium.model")
with _load_to(device):
    system = pot.createSystem(to_openmm_topology(ir), returnEnergyType="interaction_energy", device=device)
# 重原子四元组：中心键两端都是重原子、两端原子也优先重原子
heavy = {a.index for a in ir.atoms if a.atomic_number > 1}
quads = [t for t in ir.proper_torsions if t.j in heavy and t.k in heavy]
quads.sort(key=lambda t: -((t.i in heavy) + (t.l in heavy)))
q = quads[0]; quad = (q.i, q.j, q.k, q.l)
rest = openmm.CustomTorsionForce("0.5*kr*dt^2; dt=atan2(sin(theta-t0), cos(theta-t0))")
rest.addGlobalParameter("kr", 1e4); rest.addGlobalParameter("t0", 0.0)
rest.addTorsion(*quad, []); rest.setForceGroup(1); system.addForce(rest)
plat = openmm.Platform.getPlatformByName("Reference")
ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), plat)
ctx.setPositions(pos)
def egroup0():
    return ctx.getState(getEnergy=True, groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
t0 = time.time()
for _ in range(5): ctx.getState(getEnergy=True, getForces=True, groups={0})
print(f"device={device} natoms={len(ir.atoms)} quad={quad} energy+forces {(time.time()-t0)/5*1000:.1f} ms")
# 先无约束最小化
ctx.setParameter("kr", 0.0); openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 0); ctx.setParameter("kr", 1e4)
p = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
phi0 = math.degrees(dihedral_angle(*(p[i] for i in quad)))
bonds = np.array([[b.i, b.j] for b in ir.bonds])
r_ref = np.linalg.norm(p[bonds[:, 0]] - p[bonds[:, 1]], axis=1)
e_ref = egroup0()
print(f"relaxed phi0={phi0:.2f} E={e_ref:.3f}")
for s in range(1, nsteps + 1):
    tgt = phi0 + 10.0 * s
    ctx.setParameter("t0", math.radians(tgt))
    t1 = time.time(); openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 0); dt = time.time() - t1
    p = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    phi = math.degrees(dihedral_angle(*(p[i] for i in quad)))
    r = np.linalg.norm(p[bonds[:, 0]] - p[bonds[:, 1]], axis=1)
    print(f"step {s:2d} target {((tgt+180)%360)-180:8.2f} got {phi:8.2f}  dE={egroup0()-e_ref:9.3f} kJ/mol  max|dr|={np.max(np.abs(r-r_ref))*1000:.2f} pm  {dt:5.1f}s", flush=True)
