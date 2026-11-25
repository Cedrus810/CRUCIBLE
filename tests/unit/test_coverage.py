"""覆盖集诊断纯逻辑层（validation/coverage.py）与驱动脚本选择器的单元测试。

只测不依赖外部工具的部分：frcmod 解析、verdict 规则、表加载、汇总渲染
（2026-11-20 计划 §2.2；跑 antechamber/tleap 的验收属 external，不在本文件）。
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ffcompiler.validation.coverage import (
    CATEGORIES,
    VERDICT_BLOCKED,
    VERDICT_GAPS,
    VERDICT_OK,
    MoleculeRecord,
    load_table,
    parse_frcmod_gaps,
    summarize,
    verdict_from_stages,
    write_summary_csv,
)


# ---- parse_frcmod_gaps -------------------------------------------------------

FRCMOD_CLEAN = """\
remark line goes here
MASS

BOND

ANGLE
H1  CT  1   100.0  110.0
DIHE

IMPROPER
ATTN, need revision
NONBON
"""

FRCMOD_GAPS = """\
remark
MASS

BOND
ATTN, need revision
ATTN, need revision
ANGLE

DIHE
ATTN, need revision
IMPROPER

NONBON
"""


def _frcmod(tmp_path: Path, text: str) -> Path:
    f = tmp_path / "mol.frcmod"
    f.write_text(text, encoding="utf-8")
    return f


def test_parse_frcmod_clean_counts_only_improper(tmp_path):
    out = parse_frcmod_gaps(_frcmod(tmp_path, FRCMOD_CLEAN))
    assert out["n_attn"] == 1
    assert out["sections"] == {"IMPROPER": 1}


def test_parse_frcmod_gaps_counts_bond_and_dihe(tmp_path):
    out = parse_frcmod_gaps(_frcmod(tmp_path, FRCMOD_GAPS))
    assert out["n_attn"] == 3
    assert out["sections"] == {"BOND": 2, "DIHE": 1}
    assert len(out["attn_lines"]) == 3


def test_parse_frcmod_empty_raises(tmp_path):
    f = _frcmod(tmp_path, "\n")
    with pytest.raises(ValueError, match="空文件"):
        parse_frcmod_gaps(f)


# ---- verdict_from_stages -----------------------------------------------------

def test_verdict_ok_when_all_light_stages_pass():
    stages = {s: {"status": "ok"} for s in ("entry", "rotors", "frcmod", "export", "refcheck")}
    verdict, blockers = verdict_from_stages(stages)
    assert (verdict, blockers) == (VERDICT_OK, [])


def test_verdict_blocked_on_entry_fail():
    verdict, blockers = verdict_from_stages(
        {"entry": {"status": "fail", "error": "SDF 解析失败：x.sdf"}})
    assert verdict == VERDICT_BLOCKED
    assert "entry" in blockers[0]


def test_verdict_blocked_on_export_not_equivalent():
    stages = {"entry": {"status": "ok"},
              "export": {"status": "fail", "error": "prmtop leg diverged"}}
    verdict, blockers = verdict_from_stages(stages)
    assert verdict == VERDICT_BLOCKED
    assert any("export_not_equivalent" in b for b in blockers)


def test_verdict_gaps_from_frcmod_and_refcheck_but_not_skip():
    stages = {
        "entry": {"status": "ok"}, "rotors": {"status": "ok"},
        "export": {"status": "ok"},
        "frcmod": {"status": "ok", "n_attn": 2, "gaps": ["frcmod_missing_params(2)"]},
        "refcheck": {"status": "skip"},  # 未跑的阶段不产生 gap
    }
    verdict, blockers = verdict_from_stages(stages)
    assert verdict == VERDICT_GAPS
    assert blockers == ["frcmod_missing_params(2)"]


def test_verdict_stage_fail_is_gap_for_non_entry_export():
    stages = {"entry": {"status": "ok"}, "rotors": {"status": "fail", "error": "boom"}}
    verdict, blockers = verdict_from_stages(stages)
    assert verdict == VERDICT_GAPS
    assert blockers == ["rotors: boom"]


# ---- MoleculeRecord ----------------------------------------------------------

def test_record_rejects_unknown_category():
    with pytest.raises(ValueError, match="未知类别"):
        MoleculeRecord(name="x", category="not_a_cat", net_charge=0, n_atoms=1,
                       sdf_sha256="0" * 64)


def test_record_finalize_sets_verdict_and_blockers():
    rec = MoleculeRecord(name="ethanol", category="alcohol_ether", net_charge=0,
                         n_atoms=9, sdf_sha256="a" * 64)
    rec.stages = {"entry": {"status": "ok"},
                  "frcmod": {"status": "ok", "gaps": ["frcmod_missing_params(1)"],
                             "n_attn": 1}}
    rec.finalize()
    assert rec.verdict == VERDICT_GAPS
    assert rec.blockers == ["frcmod_missing_params(1)"]
    assert rec.to_dict()["verdict"] == VERDICT_GAPS


# ---- load_table --------------------------------------------------------------

def _write_table(tmp_path: Path, rows: list[dict]) -> Path:
    table = tmp_path / "molecules.json"
    table.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    return table


def test_load_table_resolves_relative_sdf(tmp_path):
    (tmp_path / "sdfs").mkdir()
    table = _write_table(tmp_path, [
        {"name": "ethanol", "category": "alcohol_ether", "sdf": "sdfs/ethanol.sdf"},
        {"name": "acetate", "category": "ion_zwitterion", "sdf": "sdfs/acetate.sdf",
         "net_charge": -1},
    ])
    rows = load_table(table)
    assert [r["name"] for r in rows] == ["ethanol", "acetate"]
    assert rows[0]["sdf"] == (tmp_path / "sdfs" / "ethanol.sdf").resolve()
    assert rows[0]["net_charge"] == 0
    assert rows[1]["net_charge"] == -1


def test_load_table_rejects_duplicate_and_missing(tmp_path):
    dup = _write_table(tmp_path, [
        {"name": "a", "category": "amine", "sdf": "a.sdf"},
        {"name": "a", "category": "amine", "sdf": "b.sdf"},
    ])
    with pytest.raises(ValueError, match="重复"):
        load_table(dup)
    missing = _write_table(tmp_path, [{"name": "b", "category": "amine"}])
    with pytest.raises(ValueError, match="缺字段"):
        load_table(missing)


def test_load_table_rejects_unknown_category_and_empty(tmp_path):
    bad_cat = _write_table(tmp_path, [{"name": "x", "category": "nope", "sdf": "x.sdf"}])
    with pytest.raises(ValueError, match="未知类别"):
        load_table(bad_cat)
    empty = _write_table(tmp_path, [])
    with pytest.raises(ValueError, match="至少一条"):
        load_table(empty)


# ---- summarize / write_summary_csv ------------------------------------------

def _sample_records() -> list[MoleculeRecord]:
    ok = MoleculeRecord(name="ethanol", category="alcohol_ether", net_charge=0,
                        n_atoms=9, sdf_sha256="a" * 64)
    ok.stages = {"entry": {"status": "ok"}}
    ok.finalize()
    gaps = MoleculeRecord(name="acetate", category="ion_zwitterion", net_charge=-1,
                          n_atoms=7, sdf_sha256="b" * 64)
    gaps.stages = {"entry": {"status": "ok"},
                   "refcheck": {"status": "ok", "gaps": ["charged_species(q=-1; mace_applicability_open)"]}}
    gaps.finalize()
    return [ok, gaps]


def test_summarize_has_category_rows_and_molecule_lines():
    text = summarize(_sample_records())
    assert "醇/醚 (alcohol_ether)" in text
    assert "离子/两性离子 (ion_zwitterion)" in text
    assert "ethanol" in text and "acetate" in text
    assert VERDICT_GAPS in text


def test_write_summary_csv_columns(tmp_path):
    out = tmp_path / "summary.csv"
    write_summary_csv(_sample_records(), out)
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].split(",")[:5] == ["name", "category", "net_charge", "n_atoms", "verdict"]
    ethanol = [ln for ln in lines[1:] if ln.startswith("ethanol,")][0]
    assert "ok" in ethanol
    acetate = [ln for ln in lines[1:] if ln.startswith("acetate,")][0]
    assert "ok_with_gaps" in acetate
    assert "charged_species" in acetate


def test_categories_match_plan_table():
    # 计划 §2.1 的八类一个不少
    assert set(CATEGORIES) == {
        "alcohol_ether", "carbonyl_ester_amide", "amine", "heteroaromatic",
        "aromatic_conjugated_fused", "small_ring", "heteroatom_halogen",
        "ion_zwitterion",
    }


# ---- 驱动脚本的 xtb_sp 阶段（不调 xtb：替换运行器） ----

def _load_driver():
    import importlib.util

    path = Path(__file__).resolve().parents[2] / "scripts" / "coverage_diagnose.py"
    spec = importlib.util.spec_from_file_location("coverage_diagnose", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _water_ir():
    from types import SimpleNamespace

    atoms = [SimpleNamespace(atomic_number=z) for z in (8, 1, 1)]
    return SimpleNamespace(atoms=atoms, total_charge=0)


def test_stage_xtb_sp_records_energy_and_gradient(tmp_path, monkeypatch):
    import numpy as np

    import ffcompiler.hessian.xtb as xtb_mod

    driver = _load_driver()
    seen = {}

    def fake(z, pos, *, method, charge, workdir, threads):
        seen.update(method=method, charge=charge, threads=threads, workdir=workdir)
        return {"energy_hartree": -5.07,
                "gradient_hartree_per_bohr": np.array([[0, 0, 3e-3], [0, 4e-3, 0], [0, 0, 0]]),
                "xtb_version": "6.7.1", "command": []}

    monkeypatch.setattr(xtb_mod, "run_xtb_singlepoint", fake)
    out = driver.stage_xtb_sp(_water_ir(), np.zeros((3, 3)), tmp_path / "xtb_sp", threads=2)
    assert out["status"] == "ok"
    assert seen == {"method": "gfn2", "charge": 0, "threads": 2, "workdir": tmp_path / "xtb_sp"}
    assert out["max_atom_gradient_hartree_per_bohr"] == 4e-3
    assert out["energy_hartree"] == -5.07


def test_stage_xtb_sp_failure_is_gap_not_blocker(tmp_path, monkeypatch):
    import numpy as np

    import ffcompiler.hessian.xtb as xtb_mod

    driver = _load_driver()

    def fake(*args, **kwargs):
        raise xtb_mod.XtbError("价电子数为奇数")

    monkeypatch.setattr(xtb_mod, "run_xtb_singlepoint", fake)
    out = driver.stage_xtb_sp(_water_ir(), np.zeros((3, 3)), tmp_path)
    assert out["status"] == "fail" and "XtbError" in out["error"]
    verdict, notes = verdict_from_stages({"entry": {"status": "ok"}, "xtb_sp": out})
    assert verdict == VERDICT_GAPS and notes[0].startswith("xtb_sp:")


# ---- 驱动脚本的 fallback / scan_probe 阶段（计划 §2.2 第 3、6 项；不调 xtb/MACE） ----

_FIX = Path(__file__).resolve().parents[1] / "fixtures"


def _ethanol_entry():
    from ffcompiler.importers import from_amber

    d = _FIX / "gaff2" / "ethanol"
    ir, pos = from_amber(d / "ethanol.prmtop", d / "ethanol.inpcrd", d / "ethanol.sdf")
    return ir, pos, d / "ethanol.sdf"


def test_stage_fallback_real_bonded_fit_no_linear(tmp_path, monkeypatch):
    import numpy as np

    import ffcompiler.hessian.xtb as xtb_mod

    driver = _load_driver()
    ir, pos, _ = _ethanol_entry()
    h = xtb_mod.read_xtb_hessian(_FIX / "hessian" / "xtb_ethanol_gfn2" / "hessian",
                                 _FIX / "hessian" / "xtb_ethanol_gfn2" / "xtbopt.xyz",
                                 method="gfn2", version="6.7.1")
    seen = {}

    def fake(z, p, *, method, charge, workdir, threads):
        seen.update(method=method, charge=charge, threads=threads)
        return h

    monkeypatch.setattr(xtb_mod, "run_xtb_hessian", fake)
    out = driver.stage_fallback(ir, pos, tmp_path, threads=3)
    assert seen == {"method": "gfnff", "charge": 0, "threads": 3}
    assert out["status"] == "ok", out
    assert out["n_bonds"] == len(ir.bonds) and out["n_angles"] == len(ir.angles)
    assert out["n_linear_angles_library"] == 0 and out["gaps"] == []
    assert np.isfinite(out["max_bond_ratio_dev"])


def test_stage_fallback_linear_angles_become_gap(tmp_path, monkeypatch):
    from types import SimpleNamespace

    import ffcompiler.bonded.assign as assign_mod
    import ffcompiler.hessian.xtb as xtb_mod

    driver = _load_driver()
    ir, pos, _ = _ethanol_entry()
    monkeypatch.setattr(xtb_mod, "run_xtb_hessian", lambda *a, **k: object())
    seen = {}

    def fake_derive(ir_, h, *, method, linear_angle_policy, input_positions_nm):
        seen.update(method=method, policy=linear_angle_policy)
        report = SimpleNamespace(
            linear_angles_kept=((0, 1, 2),),
            diagnostics=SimpleNamespace(n_imaginary=0, warnings=()),
            max_bond_ratio=1.0, class_spread={})
        return ir_, pos, report

    monkeypatch.setattr(assign_mod, "derive_bonded", fake_derive)
    out = driver.stage_fallback(ir, pos, tmp_path)
    assert seen == {"method": "mseminario", "policy": "library"}
    assert out["n_linear_angles_library"] == 1
    assert out["gaps"][0].startswith("bonded_linear_angles_library(n=1")
    verdict, notes = verdict_from_stages({"entry": {"status": "ok"}, "fallback": out})
    assert verdict == VERDICT_GAPS and notes[0].startswith("bonded_linear_angles_library")


def test_stage_fallback_nonpositive_k_is_fail(tmp_path, monkeypatch):
    import ffcompiler.bonded.assign as assign_mod
    import ffcompiler.hessian.xtb as xtb_mod

    driver = _load_driver()
    ir, pos, _ = _ethanol_entry()
    monkeypatch.setattr(xtb_mod, "run_xtb_hessian", lambda *a, **k: object())

    def boom(*a, **k):
        raise assign_mod.NonPositiveForceConstantError("mseminario", [((0, 1, 2), -1.0)])

    monkeypatch.setattr(assign_mod, "derive_bonded", boom)
    out = driver.stage_fallback(ir, pos, tmp_path)
    assert out["status"] == "fail" and "NonPositiveForceConstantError" in out["error"]


class _FakeProvider:
    calls: list = []
    row_status = "mlp"

    def __init__(self, factory, *, scan_kwargs, model_path, scan_cache, fit_bonds, **kwargs):
        _FakeProvider.calls.append({"scan_kwargs": scan_kwargs, "scan_cache": scan_cache,
                                    "fit_bonds": fit_bonds, "model_path": model_path, **kwargs})

    def parameterize(self, ir, pos, mol):
        from types import SimpleNamespace

        bond = list(_FakeProvider.calls[-1]["fit_bonds"][0])
        row = {"bond": bond, "status": _FakeProvider.row_status,
               "reason": None if _FakeProvider.row_status == "mlp" else "scan_not_converged",
               "fit_rmse": 0.1, "max_hysteresis": 0.01, "max_bond_change_pm": 1.0}
        return ir, SimpleNamespace(rotors=[row], molecule_fallback=None)


def _patch_mace(monkeypatch):
    import ffcompiler.torsion.mace as mace_mod
    import ffcompiler.torsion.mlp as mlp_mod

    _FakeProvider.calls = []
    monkeypatch.setattr(mlp_mod, "MLPTorsionProvider", _FakeProvider)
    monkeypatch.setattr(mace_mod, "default_model_path", lambda: Path("/nonexistent.model"))


def test_stage_scan_probe_picks_first_rotor_and_times(tmp_path, monkeypatch):
    from rdkit import Chem

    from ffcompiler.torsion.rotors import find_rotors

    driver = _load_driver()
    _patch_mace(monkeypatch)
    _FakeProvider.row_status = "mlp"
    ir, pos, sdf = _ethanol_entry()
    first = find_rotors(ir, Chem.MolFromMolFile(str(sdf), removeHs=False), pos)[0].bond
    out = driver.stage_scan_probe(sdf, ir, pos, tmp_path, threads=8)
    assert out["status"] == "ok" and out["gaps"] == []
    assert out["bond"] == list(first) and out["threads"] == 8
    assert out["wall_s"] >= 0.0
    call = _FakeProvider.calls[-1]
    assert call["fit_bonds"] == [tuple(first)]
    assert call["scan_kwargs"] == {"step_deg": 10.0}       # 与 G2 同网格
    assert call["scan_cache"] == tmp_path / "scans"


def test_stage_scan_probe_fallback_is_gap(tmp_path, monkeypatch):
    driver = _load_driver()
    _patch_mace(monkeypatch)
    _FakeProvider.row_status = "fallback"
    ir, pos, sdf = _ethanol_entry()
    out = driver.stage_scan_probe(sdf, ir, pos, tmp_path)
    assert out["gaps"] == ["scan_probe_fallback(scan_not_converged)"]
    verdict, _ = verdict_from_stages({"entry": {"status": "ok"}, "scan_probe": out})
    assert verdict == VERDICT_GAPS


def test_stage_scan_probe_skips_without_rotors(tmp_path, monkeypatch):
    from ffcompiler.importers import from_amber

    driver = _load_driver()
    _patch_mace(monkeypatch)
    d = _FIX / "gaff2" / "benzene"
    ir, pos = from_amber(d / "benzene.prmtop", d / "benzene.inpcrd", d / "benzene.sdf")
    out = driver.stage_scan_probe(d / "benzene.sdf", ir, pos, tmp_path)
    assert out["status"] == "skip"
    assert _FakeProvider.calls == []                     # 没有可转键就不建 provider
    verdict, _ = verdict_from_stages({"entry": {"status": "ok"}, "scan_probe": out})
    assert verdict == VERDICT_OK


@pytest.mark.external
def test_driver_light_stages_on_real_ethanol(tmp_path):
    """入口真跑 antechamber/parmchk2/tleap；下游阶段必须能读回 entry 写出的 ir.json。"""
    driver = _load_driver()
    row = {"name": "ethanol", "category": "alcohol_ether",
           "sdf": _FIX / "gaff2" / "ethanol" / "ethanol.sdf", "net_charge": 0}
    record = driver.diagnose_one(row, tmp_path,
                                 ("entry", "rotors", "frcmod", "refcheck", "xtb_sp"),
                                 force=False, xtb_threads=2)
    assert record.stages["entry"]["status"] == "ok", record.stages["entry"]
    assert (tmp_path / "ethanol" / "entry" / "ir.json").is_file()
    assert (tmp_path / "ethanol" / "entry" / "coords.json").is_file()
    assert record.stages["rotors"]["n_rotors"] == 2
    # 下游拿到的是 nm 坐标：入口几何上 GFN2 单点能收敛（2026-11-24 曾把 Å 当 nm，几何放大
    # 10 倍，SCC 不收敛）；坐标 sidecar 与 inpcrd 只差 Å→nm
    import numpy as np

    from ffcompiler.formats.amber import read_rst7
    from ffcompiler.ir import coords as ir_coords

    assert record.stages["xtb_sp"]["status"] == "ok", record.stages["xtb_sp"]
    pos_nm = ir_coords.from_json(tmp_path / "ethanol" / "entry" / "coords.json").positions_nm
    np.testing.assert_allclose(pos_nm * 10.0,
                               read_rst7(tmp_path / "ethanol" / "entry" / "ethanol.inpcrd"),
                               rtol=0, atol=1e-12)
    assert record.verdict == VERDICT_OK
    # 续跑：record.json 存在即读回、不重跑
    again = driver.diagnose_one(row, tmp_path, ("entry",), force=False)
    assert again.stages.keys() == record.stages.keys()


# ---- 导出阶段：混合精度 GROMACS 超差只记缺项（用户决定 2026-11-24，同决策 0013） ----

def _fake_compare_report(failed_legs, precision="mixed", errors=None, failed_force=False):
    from types import SimpleNamespace

    rows = [{"leg": "openmm-ir", "reference": "numpy-ir", "category": "angle",
             "E_leg": 0.4745, "E_ref": 0.4745, "delta": 1e-9, "tol": 1e-4, "passed": True}]
    for leg in failed_legs:
        rows.append({"leg": leg, "reference": "numpy-ir", "category": "angle",
                     "E_leg": 0.4734, "E_ref": 0.4745, "delta": 1.1e-3, "tol": 1.05e-4,
                     "passed": False})
    force_rows = [{"leg": "openmm-xml", "reference": "openmm-ir", "max_df": 1e-6,
                   "tol": 1e-3, "passed": not failed_force}]
    return SimpleNamespace(rows=rows, force_rows=force_rows, errors=errors or {},
                           skipped={}, metadata={"gmx": {"precision": precision}})


def _run_export(monkeypatch, tmp_path, report):
    import numpy as np

    import ffcompiler.validation.engines as engines_mod

    monkeypatch.setattr(engines_mod, "compare", lambda *a, **k: report)
    driver = _load_driver()
    return driver.stage_export(object(), np.zeros((1, 3)), tmp_path, prmtop="", inpcrd="")


def test_export_mixed_gmx_only_is_gap(tmp_path, monkeypatch):
    out = _run_export(monkeypatch, tmp_path, _fake_compare_report(["gmx"]))
    assert out["status"] == "ok"
    assert out["gaps"][0].startswith("gmx_mixed_precision(angle;")
    assert len(out["gmx_mixed_precision_failed_rows"]) == 1
    verdict, notes = verdict_from_stages({"entry": {"status": "ok"}, "export": out})
    assert verdict == VERDICT_GAPS and notes[0].startswith("gmx_mixed_precision")


@pytest.mark.parametrize("kwargs", [
    {"failed_legs": ["gmx"], "precision": "double"},          # 非混合精度：照旧 blocked
    {"failed_legs": ["gmx"], "precision": None},              # 精度未知：照旧 blocked
    {"failed_legs": ["gmx", "openmm-top"]},                   # 还有别的腿超差
    {"failed_legs": ["gmx"], "errors": {"sander": "crash"}},  # 有腿出错
    {"failed_legs": ["gmx"], "failed_force": True},           # 力超差
])
def test_export_other_failures_still_block(tmp_path, monkeypatch, kwargs):
    out = _run_export(monkeypatch, tmp_path, _fake_compare_report(**kwargs))
    assert out["status"] == "fail"
    verdict, _ = verdict_from_stages({"entry": {"status": "ok"}, "export": out})
    assert verdict == VERDICT_BLOCKED


# ---- 决策 0028：离子/两性离子参考势用 xtb GFN2 ----

def _fixture_entry(name):
    from ffcompiler.importers import from_amber

    d = _FIX / "gaff2" / name
    ir, pos = from_amber(d / f"{name}.prmtop", d / f"{name}.inpcrd", d / f"{name}.sdf")
    return ir, pos, d / f"{name}.sdf"


def test_refcheck_reference_follows_decision_0028():
    driver = _load_driver()
    neutral = driver.stage_refcheck(_fixture_entry("ethanol")[0])
    assert neutral["reference"] == "mace" and neutral["gaps"] == []
    for ion in ("acetate", "methylammonium"):
        out = driver.stage_refcheck(_fixture_entry(ion)[0])
        assert out["reference"] == "xtb-gfn2" and out["gaps"] == [], ion


def test_refcheck_zwitterion_uses_gfn2():
    from dataclasses import replace as dc_replace

    driver = _load_driver()
    ir, _, _ = _fixture_entry("ethanol")
    # 净电荷 0 但有形式电荷 ±1 → 两性离子判定（只改形式电荷，测判定逻辑）
    atoms = list(ir.atoms)
    atoms[0] = dc_replace(atoms[0], formal_charge=1)
    atoms[2] = dc_replace(atoms[2], formal_charge=-1)
    zw = dc_replace(ir, atoms=tuple(atoms))
    assert driver.stage_refcheck(zw)["reference"] == "xtb-gfn2"


def test_scan_probe_ion_uses_gfn2_provider(tmp_path, monkeypatch):
    import ffcompiler.torsion.xtb_ref as xtb_ref_mod

    driver = _load_driver()
    _patch_mace(monkeypatch)
    _FakeProvider.row_status = "mlp"
    monkeypatch.setattr(xtb_ref_mod, "xtb_binary", lambda: Path("/fake/xtb"))
    ir, pos, sdf = _fixture_entry("methylammonium")
    out = driver.stage_scan_probe(sdf, ir, pos, tmp_path)
    assert out["status"] == "ok" and out["reference"] == "xtb-gfn2"
    call = _FakeProvider.calls[-1]
    assert call["model_name"] == "xtb-gfn2" and call["model_path"] == Path("/fake/xtb")
    assert 7 in call["supported"] and 87 not in call["supported"]
