"""MVP-2 冒烟：MACE-OFF24 经 openmm-ml 在 CUDA 上加载，测单点与带二面角约束的最小化耗时。"""
import sys, time, math
from pathlib import Path
import numpy as np
import openmm, openmm.app
from openmm import unit
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from ffcompiler.importers.amber_prmtop import from_amber
from ffcompiler.exporters.openmm import to_openmm_topology

name = sys.argv[1] if len(sys.argv) > 1 else "ethanol"
d = REPO / "tests/fixtures/gaff2" / name
ir, pos = from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")
t0 = time.time()
from openmmml import MLPotential
pot = MLPotential("mace", modelPath="/home/ruigengji/MLP/mace/MACE-OFF24_medium.model")
top = to_openmm_topology(ir)
import torch, contextlib
@contextlib.contextmanager
def _load_to(device):
    # openmm-ml 1.6 的 modelPath 路径只做 torch.load(map_location=device)，
    # e3nn 编译模块的 4 个 _w3j buffer 会留在 CPU；补一次 .to(device)。
    orig = torch.load
    def patched(*a, **k):
        obj = orig(*a, **k)
        return obj.to(device) if isinstance(obj, torch.nn.Module) else obj
    torch.load = patched
    try:
        yield
    finally:
        torch.load = orig
with _load_to("cuda"):
    system = pot.createSystem(top, returnEnergyType="interaction_energy", device="cuda", **({"precision": sys.argv[2]} if len(sys.argv) > 2 else {}))
print(f"load+createSystem {time.time()-t0:.1f}s; forces={[f.__class__.__name__ for f in system.getForces()]}")
# 二面角约束：取第一条 proper 的四元组
q = ir.proper_torsions[0]
quad = (q.i, q.j, q.k, q.l)
rest = openmm.CustomTorsionForce("0.5*kr*min(dt, 2*pi-dt)^2; dt=abs(theta-t0); pi=3.141592653589793")
rest.addGlobalParameter("kr", 1e4); rest.addGlobalParameter("t0", 0.0)
rest.addTorsion(*quad, [])
rest.setForceGroup(1)
system.addForce(rest)
for plat in ("CUDA", "Reference"):
    try:
        ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName(plat))
        break
    except Exception as e:
        print("platform", plat, "failed:", e)
print("platform:", ctx.getPlatform().getName())
ctx.setPositions(pos)
t0 = time.time(); e = ctx.getState(getEnergy=True, groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole); print(f"first energy {e:.4f} kJ/mol ({time.time()-t0:.2f}s)")
t0 = time.time()
for _ in range(10): ctx.getState(getEnergy=True, getForces=True, groups={0})
print(f"energy+forces x10: {(time.time()-t0)/10*1000:.1f} ms each")
times = []
for ang in (60.0, 120.0, 180.0):
    ctx.setParameter("t0", math.radians(ang))
    t0 = time.time()
    openmm.LocalEnergyMinimizer.minimize(ctx, 1e-2, 0)
    st = ctx.getState(getEnergy=True, getPositions=True, groups={0})
    p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    from ffcompiler.ir.energy import dihedral_angle
    phi = math.degrees(dihedral_angle(*(p[i] for i in quad)))
    times.append(time.time()-t0)
    print(f"target {ang:6.1f} got {phi:8.3f} deg  E_mlp {st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.3f}  minimize {times[-1]:.1f}s")
print(f"natoms={len(ir.atoms)} mean minimize {np.mean(times):.1f}s")
