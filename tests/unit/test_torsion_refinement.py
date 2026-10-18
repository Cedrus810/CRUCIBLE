"""Optional sequential sweeps: synthetic coupled profiles, no reference-potential calls.

Written under the user's no-compute instruction. These tests have NOT been run.
The two scan blocks have features (1, 1/2) and (1/2, 1) and known target
amplitudes (1, 1). Each block is fitted on its own scan, as in the B3 experiment.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.bonded.staleness import staleness
from ffcompiler.charge.apply import assign_charges
from ffcompiler.charge.base import ChargeProvider
from ffcompiler.charge.loaded import LoadedCharges
from ffcompiler.importers import from_amber
from ffcompiler.ir import ProperTorsion
from ffcompiler.ir import io as ir_io
from ffcompiler.torsion import mlp
from ffcompiler.torsion.scan import ScanResult

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2" / "ethanol"


@pytest.fixture
def coupled(monkeypatch, tmp_path):
    base, pos = from_amber(*(FIXTURES / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(FIXTURES / "ethanol.sdf"), removeHs=False)
    rotors = mlp.find_rotors(base, mol, pos)
    assert len(rotors) == 2
    base = base.replace_terms(proper_torsions=tuple(
        ProperTorsion(*r.scan_quartet, 1, 0.0, 0.1, "gaff2") for r in rotors))
    base, _ = assign_charges(base, LoadedCharges(), pos, mol)
    bonds = [r.bond for r in rotors]
    phi = np.arange(18) * (2 * np.pi / 18)
    scans = {}
    for idx, rotor in enumerate(rotors):
        frames = np.zeros((len(phi), len(base.atoms), 3))
        # Markers for the analytic evaluator, not molecular geometries.
        frames[:, 0, 0] = idx
        frames[:, 0, 1] = phi
        scans[rotor.bond] = ScanResult(
            quartet=rotor.scan_quartet, phi_target_rad=phi, phi_rad=phi,
            energy=1.5 * (1 + np.cos(phi)), positions=frames,
            hysteresis=np.zeros(len(phi)), converged=np.ones(len(phi), dtype=bool),
            max_bond_change_pm=0.0, params={"step_deg": 20.0})

    def energy(ir, frames, excluded=None):
        theta = np.asarray(frames)[:, 0, 1]
        scan_idx = int(frames[0, 0, 0])
        result = np.zeros_like(theta)
        for term in ir.proper_torsions:
            key = frozenset((term.j, term.k))
            if key == excluded:
                continue
            idx = next(i for i, b in enumerate(bonds) if frozenset(b) == key)
            weight = 1.0 if idx == scan_idx else 0.5
            result += weight * term.force_constant * (
                1 + np.cos(term.periodicity * theta - term.phase))
        return result

    def forbidden(*args, **kwargs):
        raise AssertionError("refinement/refit must not call the reference potential or write scans")

    monkeypatch.setattr(mlp, "find_rotors", lambda *args: rotors)
    monkeypatch.setattr(mlp, "load_scan", lambda path, key: scans[tuple(map(int, path.stem.split("-")))])
    monkeypatch.setattr(mlp, "save_scan", forbidden)
    monkeypatch.setattr(mlp, "mm_energy", energy)
    monkeypatch.setattr(mlp, "mm_energy_without_bond",
                        lambda ir, bond, frames: energy(ir, frames, frozenset(bond)))

    def provider(rounds=0):
        return mlp.MLPTorsionProvider(forbidden, supported=frozenset(range(1, 119)),
                                      scan_cache=tmp_path, cache_only=True,
                                      refinement_rounds=rounds)

    return base, pos, mol, rotors, provider


def _amplitudes(ir, bonds):
    return [sum(t.force_constant * np.cos(t.phase) for t in ir.proper_torsions
                if {t.j, t.k} == set(b) and t.periodicity == 1) for b in bonds]


@pytest.mark.parametrize("value", [-1, 4, 1.5, "3", True, None])
def test_invalid_round_budget_rejected_before_loading_model(value):
    with pytest.raises(ValueError, match="refinement_rounds"):
        mlp.MLPTorsionProvider(refinement_rounds=value)


def test_default_stays_single_pass(coupled):
    base, pos, mol, rotors, provider = coupled
    fitted, report = provider().parameterize(base, pos, mol)
    assert _amplitudes(fitted, [r.bond for r in rotors]) == pytest.approx([1.45, 1.45])
    assert report.refinement is None and "refinement" not in report.to_dict()
    assert fitted.provenance["torsion"]["fit_protocol"]["refinement_rounds"] == 0


def test_three_sweeps_update_current_other_bond_and_keep_original_terms(coupled):
    base, pos, mol, rotors, provider = coupled
    fitted, report = provider(3).parameterize(base, pos, mol)
    bonds = [r.bond for r in rotors]
    # Analytic recurrence: a <- 1.5 - b/2; b <- 1.5 - a/2.
    assert _amplitudes(fitted, bonds) == pytest.approx([0.9859375, 1.00703125], abs=1e-10)
    assert [h["round"] for h in report.refinement["history"]] == [0, 1, 2, 3]
    assert report.refinement["completed_rounds"] == 3
    assert report.refinement["objective"] == "per_bond_scan"
    assert report.refinement["history"][-1]["pooled_profile_rmse"] < 0.02
    assert report.refinement["history"][0]["pooled_profile_rmse"] > 0.4
    for row in report.rotors:
        assert [t["force_constant"] for t in row["replaced_terms"]] == [0.1]
    assert staleness(fitted) == []
    assert base.proper_torsions[0].force_constant == 0.1


def test_update_order_is_canonical_even_if_rotors_arrive_reversed(coupled, monkeypatch):
    base, pos, mol, rotors, provider = coupled
    expected, _ = provider(3).parameterize(base, pos, mol)
    monkeypatch.setattr(mlp, "find_rotors", lambda *args: list(reversed(rotors)))
    actual, report = provider(3).parameterize(base, pos, mol)
    assert actual.proper_torsions == expected.proper_torsions
    assert report.refinement["order"] == [list(r.bond) for r in rotors]


def test_refit_restores_original_base_and_reuses_round_protocol_after_json(coupled, tmp_path):
    base, pos, mol, _, provider = coupled
    fitted, _ = provider(3).parameterize(base, pos, mol)
    ir_io.to_json(fitted, tmp_path / "ir.json")
    restored = ir_io.from_json(tmp_path / "ir.json")
    again, report = assign_charges(restored, LoadedCharges(), pos, mol, scan_cache=tmp_path)
    assert again.proper_torsions == fitted.proper_torsions
    assert report.torsion_refit["refinement"]["completed_rounds"] == 3
    third, _ = assign_charges(again, LoadedCharges(), pos, mol, scan_cache=tmp_path)
    assert third.proper_torsions == fitted.proper_torsions

    class Shifted(ChargeProvider):
        method = "shifted-refinement-test"

        def compute(self, ir, positions_nm, mol):
            q = np.array([a.partial_charge for a in base.atoms])
            q[0] += 0.05
            q[1] -= 0.05
            return q

        def source_info(self):
            return {"kind": "test"}

    shifted = Shifted()
    charged_base, _ = assign_charges(base, shifted, pos, mol)
    expected, _ = provider(3).parameterize(charged_base, pos, mol)
    actual, _ = assign_charges(fitted, shifted, pos, mol, scan_cache=tmp_path)
    assert actual.proper_torsions == expected.proper_torsions
    assert actual.provenance["torsion"]["mm_dependency"]["charges"] == [
        a.partial_charge for a in actual.atoms]


def test_missing_scan_refit_raises_without_reference_call(coupled, tmp_path, monkeypatch):
    base, pos, mol, rotors, provider = coupled
    fitted, _ = provider(3).parameterize(base, pos, mol)
    loader = mlp.load_scan
    missing = f"{rotors[0].bond[0]}-{rotors[0].bond[1]}"
    monkeypatch.setattr(mlp, "load_scan", lambda path, key: None if path.stem == missing else loader(path, key))
    with pytest.raises(mlp.TorsionRefitError, match=missing):
        assign_charges(fitted, LoadedCharges(), pos, mol, scan_cache=tmp_path)


def test_original_fallback_is_not_upgraded_and_fixed_term_changes_go_stale(coupled, tmp_path, monkeypatch):
    base, pos, mol, rotors, provider = coupled
    loader = mlp.load_scan
    failed = rotors[1].bond
    label = f"{failed[0]}-{failed[1]}"
    monkeypatch.setattr(mlp, "load_scan", lambda path, key: None if path.stem == label else loader(path, key))
    fitted, _ = provider(3).parameterize(base, pos, mol)
    original_row = deepcopy(fitted.provenance["torsion"]["rotors"][1])
    monkeypatch.setattr(mlp, "load_scan", loader)
    again, _ = assign_charges(fitted, LoadedCharges(), pos, mol, scan_cache=tmp_path)
    assert again.provenance["torsion"]["rotors"][1] == original_row
    assert again.provenance["torsion"]["rotors"][1]["status"] == "fallback"
    assert staleness(again) == []
    modified = again.replace_terms(proper_torsions=tuple(
        replace(t, force_constant=t.force_constant + 0.2) if {t.j, t.k} == set(failed) else t
        for t in again.proper_torsions))
    assert any("固定 proper" in item for item in staleness(modified))


def test_report_mutation_cannot_rewrite_ir_provenance(coupled):
    base, pos, mol, _, provider = coupled
    fitted, report = provider(3).parameterize(base, pos, mol)
    saved = deepcopy(fitted.provenance)
    report.rotors[0]["replaced_terms"][0]["force_constant"] = 999.0
    report.refinement["history"][0]["bonds"].clear()
    assert fitted.provenance == saved


def test_legacy_ir_without_fit_protocol_stays_single_pass(coupled, tmp_path):
    base, pos, mol, _, provider = coupled
    fitted, _ = provider().parameterize(base, pos, mol)
    prov = deepcopy(fitted.provenance["torsion"])
    prov.pop("fit_protocol")
    legacy = fitted.replace_terms(provenance_update={"torsion": prov})
    again, report = assign_charges(legacy, LoadedCharges(), pos, mol, scan_cache=tmp_path)
    assert again.proper_torsions == fitted.proper_torsions
    assert report.torsion_refit.get("refinement") is None


def test_unknown_or_mismatched_protocol_is_rejected(coupled, tmp_path):
    base, pos, mol, _, provider = coupled
    fitted, _ = provider(3).parameterize(base, pos, mol)
    for patch in ({"name": "unknown"}, {"order": []}):
        prov = deepcopy(fitted.provenance["torsion"])
        prov["fit_protocol"].update(patch)
        bad = fitted.replace_terms(provenance_update={"torsion": prov})
        assert any("fit_protocol" in item for item in staleness(bad))
        with pytest.raises(mlp.TorsionRefitError, match="fit_protocol"):
            assign_charges(bad, LoadedCharges(), pos, mol, scan_cache=tmp_path)


def test_no_successful_bonds_does_not_run_empty_sweeps(coupled, monkeypatch):
    base, pos, mol, _, provider = coupled
    monkeypatch.setattr(mlp, "load_scan", lambda *args: None)
    fitted, report = provider(3).parameterize(base, pos, mol)
    assert fitted.proper_torsions == base.proper_torsions
    assert report.refinement["completed_rounds"] == 0
    assert len(report.refinement["history"]) == 1
    assert report.refinement["history"][0]["pooled_profile_rmse"] is None


def test_non_finite_profile_does_not_return_partial_refinement(coupled, monkeypatch):
    base, pos, mol, _, provider = coupled
    monkeypatch.setattr(mlp, "mm_energy", lambda ir, frames: np.full(len(frames), np.nan))
    with pytest.raises(mlp.TorsionRefitError, match="non-finite refinement profile"):
        provider(3).parameterize(base, pos, mol)
    assert base.proper_torsions[0].force_constant == 0.1
