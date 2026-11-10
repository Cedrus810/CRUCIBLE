"""MVP-2 冒烟 5：MACE-POLAR-1-M + ASE FixInternals 硬约束 + LBFGS，10° 顺序扫描。
用法: python smoke_polar_scan.py <fixture> <n_steps> [fmax_eV_per_A]"""
import sys, time
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from ffcompiler.importers.amber_prmtop import from_amber
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import LBFGS
from mace.calculators import MACECalculator

name, nsteps = sys.argv[1], int(sys.argv[2])
fmax = float(sys.argv[3]) if len(sys.argv) > 3 else 0.005
d = REPO / "tests/fixtures/gaff2" / name
ir, pos = from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")
calc = MACECalculator(model_paths="/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model", device="cpu",
                      default_dtype="float64", model_type="PolarMACE")
atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms], positions=pos * 10.0)
atoms.info.update(charge=ir.total_charge, spin=1)
atoms.calc = calc
heavy = {a.index for a in ir.atoms if a.atomic_number > 1}
quads = [t for t in ir.proper_torsions if t.j in heavy and t.k in heavy]
quads.sort(key=lambda t: -((t.i in heavy) + (t.l in heavy)))
q = quads[0]; quad = [q.i, q.j, q.k, q.l]
EV = 96.48533212
t0 = time.time(); opt = LBFGS(atoms, logfile=None); opt.run(fmax=fmax, steps=2000)
print(f"{name} quad={quad} charge={ir.total_charge} relax0 {time.time()-t0:.1f}s nsteps={opt.nsteps}", flush=True)
e0 = atoms.get_potential_energy() * EV
phi0 = atoms.get_dihedral(*quad)
bonds = np.array([[b.i, b.j] for b in ir.bonds])
r0 = np.linalg.norm(atoms.positions[bonds[:, 0]] - atoms.positions[bonds[:, 1]], axis=1)
# 转动的一侧：k 那侧（断开 j-k 后含 k 的连通分量）
adj = {i: set() for i in range(len(atoms))}
for i, j in bonds: adj[i].add(j); adj[j].add(i)
side, stack = {quad[2]}, [quad[2]]
while stack:
    a = stack.pop()
    for b in adj[a]:
        if b not in side and not (a == quad[2] and b == quad[1]):
            side.add(b); stack.append(b)
mask = [i in side for i in range(len(atoms))]
for s in range(1, nsteps + 1):
    tgt = (phi0 + 10.0 * s) % 360.0
    atoms.set_constraint()
    atoms.set_dihedral(*quad, tgt, mask=mask)
    atoms.set_constraint(FixInternals(dihedrals_deg=[[atoms.get_dihedral(*quad), quad]]))
    t1 = time.time(); opt = LBFGS(atoms, logfile=None); conv = opt.run(fmax=fmax, steps=2000); dt = time.time() - t1
    r = np.linalg.norm(atoms.positions[bonds[:, 0]] - atoms.positions[bonds[:, 1]], axis=1)
    print(f"step {s:2d} target {tgt:7.2f} got {atoms.get_dihedral(*quad):7.2f} conv={conv} it={opt.nsteps:3d} "
          f"dE={atoms.get_potential_energy()*EV-e0:8.3f} kJ/mol max|dr|={np.max(np.abs(r-r0))*100:.2f} pm {dt:5.1f}s", flush=True)
