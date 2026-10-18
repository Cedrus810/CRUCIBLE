"""MVP-2 Review Focus 6：任意相位（δ ∉ {0, π}）经 AMBER / GROMACS 文件导出后与 OpenMM 内存体系一致。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openmm.app as app
import pytest
from openmm import unit

from ffcompiler.exporters.amber import write_amber
from ffcompiler.exporters.gromacs import write_gromacs
from ffcompiler.exporters.openmm import to_openmm_system
from ffcompiler.importers import from_amber
from ffcompiler.ir import ProperTorsion
from ffcompiler.validation.engines import reference_context

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _energy(system, pos) -> float:
    ctx = reference_context(system)
    ctx.setPositions(pos)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def _phased_ethanol():
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    q = next((t.i, t.j, t.k, t.l) for t in ir.proper_torsions if {t.j, t.k} == {1, 2})
    kept = tuple(t for t in ir.proper_torsions if {t.j, t.k} != {1, 2})
    new = tuple(ProperTorsion(*q, periodicity=n, phase=0.37 * n - 1.1, force_constant=0.8 + 0.1 * n,
                              source="mlp") for n in range(1, 7))
    return ir.replace_terms(proper_torsions=kept + new), pos


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_amber_and_gromacs_files_keep_arbitrary_phases(tmp_path, seed) -> None:
    ir, pos = _phased_ethanol()
    x = pos + np.random.default_rng(seed).normal(scale=0.01, size=pos.shape)
    ref = _energy(to_openmm_system(ir), x)
    tol = 1e-4 + 1e-7 * abs(ref)  # 经文件容差（决策 0004）

    amber = write_amber(ir, pos, tmp_path / "amber")
    prm = app.AmberPrmtopFile(str(amber["prmtop"]))
    assert _energy(prm.createSystem(nonbondedMethod=app.NoCutoff, constraints=None), x) == pytest.approx(ref, abs=tol)

    gmx = write_gromacs(ir, pos, tmp_path / "gmx")
    top = app.GromacsTopFile(str(gmx["top"]))
    assert _energy(top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None), x) == pytest.approx(ref, abs=tol)
