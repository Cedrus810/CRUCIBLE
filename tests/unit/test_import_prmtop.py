"""T5：prmtop 导入器的必测用例。"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import openmm
import openmm.app
import pytest
from openmm import unit

from ffcompiler.formats.amber import read_prmtop, write_prmtop
from ffcompiler.importers import from_amber
from ffcompiler.ir.energy import evaluate
from ffcompiler.validation.engines import reference_context

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
    "acetonitrile",
)


def _fixture(name: str) -> tuple[Path, Path, Path]:
    mdir = FIXTURES / name
    return mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"


@pytest.mark.parametrize("name", NAMES)
def test_all_fixtures_import_and_validate(name: str) -> None:
    prmtop, inpcrd, sdf = _fixture(name)
    ir, positions = from_amber(prmtop, inpcrd, sdf)
    assert positions.shape == (len(ir.atoms), 3)
    assert ir.forcefield_family == "amber"


def _dihedral_entries(prm):
    """(数组名, 起始下标) 列表：每个二面角条目在原始数组中的位置。"""
    out = []
    for flag in ("DIHEDRALS_WITHOUT_HYDROGEN", "DIHEDRALS_INC_HYDROGEN"):
        for k in range(0, len(prm[flag]), 5):
            out.append((flag, k))
    return out


def _is_active_proper(prm, flag, k):
    return prm[flag][k + 2] > 0 and prm[flag][k + 3] > 0


def test_ethanol_counts_and_bond_conversion() -> None:
    prmtop, inpcrd, sdf = _fixture("ethanol")
    ir, _ = from_amber(prmtop, inpcrd, sdf)
    prm = read_prmtop(prmtop)

    assert len(ir.atoms) == prm.pointer("NATOM")
    assert len(ir.bonds) == prm.pointer("NBONH") + prm.pointer("MBONA")
    assert len(ir.angles) == prm.pointer("NTHETH") + prm.pointer("MTHETA")
    n_propers = sum(1 for flag, k in _dihedral_entries(prm) if prm[flag][k + 3] > 0)
    assert len(ir.proper_torsions) == n_propers

    # C0-C1 bond: k_IR = 2 × K_prmtop × 4.184 × 100（1e-12 相对误差）
    entries = prm["BONDS_WITHOUT_HYDROGEN"] + prm["BONDS_INC_HYDROGEN"]
    i, j, t = next(entries[k:k + 3] for k in range(0, len(entries), 3)
                   if {entries[k] // 3, entries[k + 1] // 3} == {0, 1})
    ir_bond = next(b for b in ir.bonds if {b.i, b.j} == {0, 1})
    expected_k = 2.0 * prm["BOND_FORCE_CONSTANT"][t - 1] * 4.184 * 100.0
    assert ir_bond.k == pytest.approx(expected_k, rel=1e-12)
    assert ir_bond.r0 == pytest.approx(prm["BOND_EQUIL_VALUE"][t - 1] * 0.1, rel=1e-12)
    assert ir_bond.source == "gaff2-library"
    assert all("amber_types" in b.meta for b in ir.bonds)


def test_benzene_impropers() -> None:
    prmtop, inpcrd, sdf = _fixture("benzene")
    ir, _ = from_amber(prmtop, inpcrd, sdf)
    assert len(ir.impropers) == 6  # GAFF2 对每个芳香 C 一个
    seen_centers = set()
    for imp in ir.impropers:
        members = (imp.i, imp.j, imp.k, imp.l)
        # 决策 0002：中心原子可以在任何位置，但必须与其余三个成键
        centers = [
            m for m in members
            if all(
                any({b.i, b.j} == {m, o} for b in ir.bonds)
                for o in members if o != m
            )
        ]
        assert len(centers) == 1, f"quartet {members} centers: {centers}"
        seen_centers.add(centers[0])
    # 6 个连通中心 = 6 个芳香 C
    assert {ir.atoms[c].atom_type for c in seen_centers} == {"ca"}
    assert len(seen_centers) == 6


@pytest.mark.parametrize(
    "name,charge,bad_elem",
    (("acetate", -1, "O"), ("methylammonium", 1, "N")),
)
def test_charged_molecules(name: str, charge: int, bad_elem: str) -> None:
    prmtop, inpcrd, sdf = _fixture(name)
    ir, _ = from_amber(prmtop, inpcrd, sdf)
    assert ir.total_charge == charge
    charged = [atom for atom in ir.atoms if atom.formal_charge != 0]
    assert len(charged) == 1
    assert charged[0].element == bad_elem
    assert abs(sum(a.partial_charge for a in ir.atoms) - charge) <= 1e-4


def test_tampered_scee_raises(tmp_path: Path) -> None:
    """一个计算 1-4 的条目改用 SCEE=1.0 的新类型 → 非均匀 1-4 缩放被拒。"""
    prmtop, inpcrd, sdf = _fixture("ethanol")
    prm = read_prmtop(prmtop)
    flag, k = next((f, k) for f, k in _dihedral_entries(prm) if _is_active_proper(prm, f, k))
    t = prm[flag][k + 4] - 1
    for name in ("DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PERIODICITY", "DIHEDRAL_PHASE",
                 "SCNB_SCALE_FACTOR"):
        prm[name].append(prm[name][t])
    prm["SCEE_SCALE_FACTOR"].append(1.0)
    prm["POINTERS"][17] += 1  # NPTRA
    prm[flag][k + 4] = prm.pointer("NPTRA")
    bad = tmp_path / "bad.prmtop"
    write_prmtop(prm, bad)

    with pytest.raises(ValueError, match="SCEE"):
        from_amber(bad, inpcrd, sdf)


@pytest.mark.parametrize("mode", ["all_ignored", "one_missing", "duplicate"])
def test_invalid_active_14_carriers_rejected(tmp_path, mode):
    prmtop, inpcrd, sdf = _fixture("ethanol")
    prm = read_prmtop(prmtop)
    active = [(f, k) for f, k in _dihedral_entries(prm) if _is_active_proper(prm, f, k)]
    if mode == "all_ignored":
        for f, k in active:
            prm[f][k + 2] = -prm[f][k + 2]
    elif mode == "one_missing":
        f, k = active[0]
        prm[f][k + 2] = -prm[f][k + 2]
    else:
        f, k = active[0]
        prm[f].extend(prm[f][k:k + 5])
        prm["POINTERS"][6 if f == "DIHEDRALS_INC_HYDROGEN" else 7] += 1  # NPHIH / MPHIA
        if f == "DIHEDRALS_WITHOUT_HYDROGEN":
            prm["POINTERS"][14] += 1  # NPHIA
    bad = tmp_path / "bad.prmtop"
    write_prmtop(prm, bad)
    with pytest.raises(ValueError, match="active 1-4 carriers"):
        from_amber(bad, inpcrd, sdf)


def test_non_lorentz_berthelot_lj_rejected(tmp_path: Path) -> None:
    """非对角 LJ 系数偏离 LB 组合 → 拒绝（IR 只支持 LB，不默默假设）。"""
    prmtop, inpcrd, sdf = _fixture("ethanol")
    prm = read_prmtop(prmtop)
    ntypes = prm.pointer("NTYPES")
    off_diag = prm["NONBONDED_PARM_INDEX"][ntypes * 0 + 1] - 1
    prm["LENNARD_JONES_ACOEF"][off_diag] *= 1.01
    bad = tmp_path / "bad.prmtop"
    write_prmtop(prm, bad)
    with pytest.raises(ValueError, match="Lorentz"):
        from_amber(bad, inpcrd, sdf)


def test_importer_does_not_load_parmed() -> None:
    """决策 0015：导入器与导出器不依赖 ParmEd。"""
    import subprocess
    import sys

    code = (
        "import sys; from ffcompiler.importers import from_amber; "
        "from ffcompiler.exporters.amber import write_amber; "
        "print('parmed' in sys.modules)"
    )
    # pytest.ini 的 pythonpath 不会传给子进程：显式给出 src
    import os

    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"))
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         env=env, check=True).stdout.strip()
    assert out == "False"


def test_element_mismatch_raises(tmp_path: Path) -> None:
    eth_prmtop, eth_inpcrd, eth_sdf = _fixture("ethanol")
    # 原子数不变，交换第 0/2 个原子的元素符号（C ↔ O）
    lines = eth_sdf.read_text().splitlines()
    n_atoms = int(lines[3][:3])
    for idx in (0, 2):
        row = 4 + idx
        sym = lines[row][31:34]
        lines[row] = lines[row][:31] + ("O " if sym.strip() == "C" else "C ") + lines[row][34:]
    bad_sdf = tmp_path / "swapped.sdf"
    bad_sdf.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="mismatch"):
        from_amber(eth_prmtop, eth_inpcrd, bad_sdf)


@pytest.mark.parametrize("name", NAMES)
def test_import_energy_matches_openmm_amber_prmtop(name: str) -> None:
    """导入器和 numpy 参考能量互相印证：evaluate vs OpenMM 读原 prmtop（经文件容差）。"""
    prmtop, inpcrd, sdf = _fixture(name)
    ir, positions = from_amber(prmtop, inpcrd, sdf)

    amber = openmm.app.AmberPrmtopFile(str(prmtop))
    system = amber.createSystem(
        nonbondedMethod=openmm.app.NoCutoff, constraints=None, rigidWater=False,
        removeCMMotion=False,
    )
    integrator = openmm.VerletIntegrator(0.001)
    context = reference_context(system)
    context.setPositions([openmm.Vec3(*p) for p in positions] * unit.nanometer)
    e_omm = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )

    e_ref = evaluate(ir, positions)["total"]
    tol = 1e-4 + 1e-7 * abs(e_ref)
    assert abs(e_omm - e_ref) <= tol, f"{name}: OpenMM {e_omm} vs IR {e_ref} (tol {tol})"
