"""假设：sander 读 prmtop 时把 |phase−π| ≤ 4e-3 的相位改成 π。用 numpy 重算 dihed，看能否解释 sander 的 1.8638。"""
import math
from pathlib import Path

from ffcompiler.importers import from_amber
from ffcompiler.ir import from_json
from ffcompiler.ir.energy import dihedral_angle

REPO = Path(__file__).resolve().parents[2]
mdir = REPO / "tests/fixtures/gaff2/nma"
_, pos = from_amber(*(mdir / f"nma.{e}" for e in ("prmtop", "inpcrd", "sdf")))
ir = from_json(REPO / "jobs/gate_mvp2_20261008_3/ir/nma.json")
KCAL = 4.184


def e(terms, snap):
    tot = 0.0
    for t in terms:
        ph = t.phase
        if snap and abs(ph - math.pi) <= 4e-3:
            ph = math.copysign(math.pi, ph)
        phi = dihedral_angle(*(pos[a] for a in (t.i, t.j, t.k, t.l)))
        tot += t.force_constant * (1 + math.cos(t.periodicity * phi - ph))
    return tot / KCAL


terms = list(ir.proper_torsions) + list(ir.impropers)
for t in terms:
    if abs(t.phase - math.pi) <= 4e-3:
        print("snapped:", (t.i, t.j, t.k, t.l), "n", t.periodicity, "phase", t.phase, "k kJ", t.force_constant)
print(f"numpy exact  {e(terms, False):.6f} kcal/mol")
print(f"numpy snap   {e(terms, True):.6f} kcal/mol")
print("sander-ours  1.8638 (4 位小数)")
