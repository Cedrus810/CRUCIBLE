"""§C1/§C2 纯逻辑层测试：留出盆地划分（validation/holdout.py）与替换变体（variants.py）。

盆地划分用合成 ScanResult（两盆地/单盆地）；变体用 ethanol 固件 IR 加合成改动，
验证"只换被点名部分、其余保留 base"与对照矩阵的完备性。不调外部工具、不做 MD。
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.importers import from_amber
from ffcompiler.ir.schema import ForceFieldIR
from ffcompiler.torsion.scan import ScanResult
from ffcompiler.validation.holdout import basin_split, split_summary
from ffcompiler.validation.variants import PARTS, comparison_matrix, substitute_terms

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ETHANOL = FIXTURES / "ethanol"
HEAVY8 = np.ones(8, dtype=bool)
N = 12


def _scan_two_basins(n=N):
    """帧 0..n-1：偶数帧几何 A、奇数帧几何 B（内部位移 0.3 nm → 两个盆地）。"""
    rng = np.random.default_rng(1)
    base = rng.normal(size=(1, 8, 3)) * 0.1
    bent = base.copy()
    bent[0, 4:, 1] += 0.3
    positions = np.repeat(base, n, axis=0)
    positions[1::2] = bent[0]
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    return ScanResult(quartet=(0, 1, 2, 3), phi_target_rad=phi, phi_rad=phi,
                      energy=np.arange(n) * 0.5, positions=positions,
                      hysteresis=np.zeros(n), converged=np.ones(n, dtype=bool),
                      max_bond_change_pm=0.0, params={"step_deg": 30.0})


# ---- basin_split -------------------------------------------------------------

def test_basin_split_whole_basin_in_or_out():
    scan = _scan_two_basins()  # 两盆地各 6 帧；默认 cap=4 放不下任一盆地
    split = basin_split(scan, HEAVY8, max_holdout_frac=0.6)
    assert split["n_basins"] == 2
    hold = set(split["holdout_indices"])
    fit = set(split["fit_indices"])
    assert hold and fit and not (hold & fit)
    assert hold | fit == set(range(N))
    # 整盆地：奇偶各成一个盆地，holdout 必须是其中一整个
    basins = {b: {s for s, bb in enumerate(split["basin_of_frame"]) if bb == b}
              for b in range(2)}
    assert hold in basins.values()


def test_basin_split_respects_fraction_cap():
    scan = _scan_two_basins(n=10)  # 两盆地各 5 帧，cap = floor(0.4×10) = 4 < 5
    split = basin_split(scan, HEAVY8)
    assert split["holdout_indices"] == []
    assert split["note"] and "cap" in split["note"]
    # 放宽上限后可划分
    loose = basin_split(scan, HEAVY8, max_holdout_frac=0.6)
    assert len(loose["holdout_indices"]) == 5


def test_basin_split_single_basin_notes_empty_holdout():
    scan = _scan_two_basins()
    same = ScanResult(quartet=scan.quartet, phi_target_rad=scan.phi_target_rad,
                      phi_rad=scan.phi_rad, energy=scan.energy,
                      positions=np.repeat(scan.positions[:1], N, axis=0),
                      hysteresis=scan.hysteresis, converged=scan.converged,
                      max_bond_change_pm=0.0, params={"step_deg": 30.0})
    split = basin_split(same, HEAVY8)
    assert split["n_basins"] == 1
    assert split["holdout_indices"] == [] and split["note"]
    assert split_summary(split).startswith("basins=1 fit=12 holdout=0")


def test_basin_split_deterministic():
    scan = _scan_two_basins()
    a = basin_split(scan, HEAVY8)
    b = basin_split(scan, HEAVY8)
    assert a == b


def test_basin_split_rejects_mismatched_frames():
    scan = _scan_two_basins()
    broken = ScanResult(quartet=scan.quartet, phi_target_rad=scan.phi_target_rad[:3],
                        phi_rad=scan.phi_rad[:3], energy=scan.energy[:3],
                        positions=scan.positions, hysteresis=scan.hysteresis[:3],
                        converged=scan.converged[:3], max_bond_change_pm=0.0, params={})
    with pytest.raises(ValueError, match="帧数不一致"):
        basin_split(broken, HEAVY8)


# ---- substitute_terms / comparison_matrix ------------------------------------

@pytest.fixture(scope="module")
def ethanol_pair():
    base, pos = from_amber(*(ETHANOL / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    # 合成 final：bonded ×1.1、扭转项清零重写、电荷守恒扰动、ε ×0.9
    final = base.replace_terms(
        bonds=tuple(replace(b, k=b.k * 1.1) for b in base.bonds),
        angles=tuple(replace(a, force_constant=a.force_constant * 1.1) for a in base.angles),
        proper_torsions=tuple(replace(t, force_constant=t.force_constant + 0.05)
                              for t in base.proper_torsions),
        atoms=tuple(
            replace(a,
                    partial_charge=(a.partial_charge + (0.01 if a.index % 2 else -0.01)
                                    * (1 if a.index < len(base.atoms) - 1 else 0.0)),
                    epsilon=a.epsilon * 0.9)
            for a in base.atoms),
    )
    return base, final


def test_substitute_only_named_part_changes(ethanol_pair):
    base, final = ethanol_pair
    bonded = substitute_terms(base, final, {"bonded"})
    assert bonded.bonds == final.bonds and bonded.angles == final.angles
    assert bonded.proper_torsions == base.proper_torsions
    assert bonded.atoms == base.atoms  # 电荷/LJ 未换

    torsion = substitute_terms(base, final, {"torsion"})
    assert torsion.proper_torsions == final.proper_torsions
    assert torsion.bonds == base.bonds

    charge = substitute_terms(base, final, {"charge"})
    assert [a.partial_charge for a in charge.atoms] == \
           [a.partial_charge for a in final.atoms]
    assert [(a.sigma, a.epsilon) for a in charge.atoms] == \
           [(a.sigma, a.epsilon) for a in base.atoms]  # vdw 未换

    vdw = substitute_terms(base, final, {"vdw"})
    assert [a.epsilon for a in vdw.atoms] == [a.epsilon * 0.9 for a in base.atoms]
    assert [a.partial_charge for a in vdw.atoms] == [a.partial_charge for a in base.atoms]

    full = substitute_terms(base, final, PARTS)
    assert full.bonds == final.bonds
    assert [a.partial_charge for a in full.atoms] == [a.partial_charge for a in final.atoms]


def test_substitute_empty_parts_returns_base_terms(ethanol_pair):
    base, final = ethanol_pair
    same = substitute_terms(base, final, set())
    assert same.bonds == base.bonds and same.atoms == base.atoms
    assert same.provenance["variant"]["parts"] == []


def test_substitute_rejects_unknown_and_mismatched_charge(ethanol_pair):
    base, final = ethanol_pair
    with pytest.raises(ValueError, match="未知 parts"):
        substitute_terms(base, final, {"improper_scaling"})
    # 声明总电荷不同（+1）的 IR：replace_terms 不改 total_charge，需整体重建
    wrong_q = final.replace_terms(
        atoms=tuple(replace(a, partial_charge=a.partial_charge + 1.0 / len(final.atoms))
                    for a in final.atoms))
    wrong_q = ForceFieldIR(
        name=wrong_q.name, atoms=wrong_q.atoms, bonds=wrong_q.bonds,
        angles=wrong_q.angles, proper_torsions=wrong_q.proper_torsions,
        impropers=wrong_q.impropers, exclusions=wrong_q.exclusions,
        pairs14=wrong_q.pairs14, total_charge=1,
        forcefield_family=wrong_q.forcefield_family,
        combination_rule=wrong_q.combination_rule,
        coulomb_14_scale=wrong_q.coulomb_14_scale,
        vdw_14_scale=wrong_q.vdw_14_scale, provenance=wrong_q.provenance)
    with pytest.raises(ValueError, match="总电荷不一致"):
        substitute_terms(base, wrong_q, {"charge"})


def test_comparison_matrix_covers_plan_variants(ethanol_pair):
    base, final = ethanol_pair
    matrix = comparison_matrix(base, final)
    assert set(matrix) == {"gaff2", "bonded", "torsion", "charge", "vdw", "full"}
    assert matrix["gaff2"] is base
    assert matrix["full"].bonds == final.bonds
    assert matrix["torsion"].bonds == base.bonds
    for name, ir in matrix.items():
        assert ir.total_charge == base.total_charge, name
