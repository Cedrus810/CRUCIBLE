"""审阅用只读检查：带 box_nm 的 IR 导出后，引擎能否读入、能量是否与 NumPy 参考一致。"""
import dataclasses, sys
from pathlib import Path
import numpy as np
import openmm, openmm.app as app, openmm.unit as u
from ffcompiler.importers import from_amber
from ffcompiler.exporters.amber import write_amber
from ffcompiler.exporters.openmm import to_openmm_system
from ffcompiler.ir.energy import evaluate

out = Path(sys.argv[1])
fx = Path("tests/fixtures/gaff2/ethanol")
ir, pos = from_amber(fx / "ethanol.prmtop", fx / "ethanol.inpcrd", fx / "ethanol.sdf")
box = ((3.0, 0, 0), (0, 3.0, 0), (0, 0, 3.0))
bir = dataclasses.replace(ir, box_nm=box)
wrapped = np.mod(pos, 3.0)

def omm_energy(system, x):
    ctx = openmm.Context(system, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(x)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)

print("numpy  contiguous :", evaluate(bir, pos)["total"])
print("numpy  wrapped    :", evaluate(bir, wrapped)["total"])
s = to_openmm_system(bir)
print("openmm contiguous :", omm_energy(s, pos))
print("openmm wrapped    :", omm_energy(s, wrapped))
print("openmm forces usesPBC:", [(f.__class__.__name__, f.usesPeriodicBoundaryConditions()) for f in s.getForces()])

files = write_amber(bir, pos, out / "amber")
try:
    p = app.AmberPrmtopFile(str(files["prmtop"]))
    print("OpenMM AmberPrmtopFile read OK; topology box:", p.topology.getPeriodicBoxVectors())
except Exception as e:
    print("OpenMM AmberPrmtopFile FAILED:", type(e).__name__, e)
try:
    import parmed
    st = parmed.load_file(str(files["prmtop"]), str(files["inpcrd"]))
    print("ParmEd read OK; box:", st.box)
except Exception as e:
    print("ParmEd FAILED:", type(e).__name__, e)
