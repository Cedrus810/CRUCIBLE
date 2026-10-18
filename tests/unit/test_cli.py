"""T10：CLI 的必测用例。"""

from __future__ import annotations

from pathlib import Path

import pytest

from ffcompiler.cli import main
from ffcompiler.validation.engines import CrossEngineReport

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"


def _ethanol_args() -> list[str]:
    mdir = FIXTURES / "ethanol"
    return [
        "--prmtop", str(mdir / "ethanol.prmtop"),
        "--inpcrd", str(mdir / "ethanol.inpcrd"),
        "--sdf", str(mdir / "ethanol.sdf"),
    ]


def test_convert_all_three_formats(tmp_path: Path) -> None:
    outdir = tmp_path / "out"
    rc = main([
        "convert", *_ethanol_args(),
        "--output", "gromacs", "amber", "openmm",
        "--outdir", str(outdir),
    ])
    assert rc == 0
    assert (outdir / "openmm" / "system.xml").is_file()
    assert (outdir / "openmm" / "molecule.pdb").is_file()
    assert (outdir / "amber" / "molecule.prmtop").is_file()
    assert (outdir / "amber" / "molecule.inpcrd").is_file()
    assert (outdir / "gromacs" / "molecule.itp").is_file()
    assert (outdir / "gromacs" / "topol.top").is_file()
    assert (outdir / "gromacs" / "conf.gro").is_file()


def test_convert_single_output_and_ir_json(tmp_path: Path) -> None:
    outdir = tmp_path / "out"
    ir_json = tmp_path / "ir.json"
    rc = main([
        "convert", *_ethanol_args(),
        "--output", "amber", "--outdir", str(outdir), "--ir-json", str(ir_json),
    ])
    assert rc == 0
    assert (outdir / "amber" / "molecule.prmtop").is_file()
    assert not (outdir / "gromacs").exists()  # 只要 amber 就不生成 GROMACS 文件
    assert not (outdir / "openmm").exists()
    assert ir_json.is_file()


@pytest.mark.external
def test_check_equivalence_ethanol_passes(tmp_path: Path) -> None:
    rc = main(["check-equivalence", *_ethanol_args(), "--workdir", str(tmp_path / "eq")])
    assert rc == 0


def test_check_equivalence_uses_source_reference(tmp_path, monkeypatch):
    import ffcompiler.validation.engines as engines_mod
    captured = {}
    def compare(ir, positions, workdir, **kwargs):
        captured.update(kwargs)
        result = CrossEngineReport()
        result.rows.append({"leg": "source", "reference": "ir", "category": "total",
                            "E_leg": 1.0, "E_ref": 2.0, "delta": 1.0, "tol": 1e-4, "passed": False})
        return result
    monkeypatch.setattr(engines_mod, "compare", compare)
    assert main(["check-equivalence", *_ethanol_args(), "--workdir", str(tmp_path)]) == 1
    assert captured["tleap_prmtop"] == FIXTURES / "ethanol/ethanol.prmtop"
    assert captured["tleap_inpcrd"] == FIXTURES / "ethanol/ethanol.inpcrd"


def test_check_equivalence_fails_with_bad_ir(tmp_path: Path, monkeypatch) -> None:
    """坏 IR 场景 → 退出码 1（monkeypatch compare 返回失败报告）。"""
    failed_report = CrossEngineReport()
    failed_report.rows.append({
        "leg": "openmm-prmtop", "reference": "numpy-ir", "category": "bond",
        "E_leg": 1.0, "E_ref": 2.0, "delta": 1.0, "tol": 1e-4, "passed": False,
    })
    import ffcompiler.validation.engines as engines_mod

    monkeypatch.setattr(engines_mod, "compare", lambda *a, **k: failed_report)
    rc = main(["check-equivalence", *_ethanol_args(), "--workdir", str(tmp_path / "eq")])
    assert rc == 1


def test_check_equivalence_error_exit_1(tmp_path: Path) -> None:
    """prmtop 与 sdf 不配套 → 导入报错 → 退出码 1。"""
    mdir = FIXTURES / "acetate"
    rc = main([
        "check-equivalence",
        "--prmtop", str(FIXTURES / "ethanol" / "ethanol.prmtop"),
        "--inpcrd", str(FIXTURES / "ethanol" / "ethanol.inpcrd"),
        "--sdf", str(mdir / "acetate.sdf"),
        "--workdir", str(tmp_path / "eq"),
    ])
    assert rc == 1


def test_missing_arguments_exit_2(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["convert", "--prmtop", "x"])  # 缺 --inpcrd/--sdf/--output/--outdir
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "usage" in err.lower()


def test_version_flag(capsys) -> None:
    import pytest

    from ffcompiler import __version__
    from ffcompiler.cli import main

    with pytest.raises(SystemExit) as info:
        main(["--version"])
    assert info.value.code == 0
    assert capsys.readouterr().out.strip() == f"ffcompiler {__version__}"


def test_derive_torsions_with_mm_reference(tmp_path: Path) -> None:
    import json

    from ffcompiler.ir import io as ir_io

    outdir = tmp_path / "out"
    rc = main([
        "derive-torsions", *_ethanol_args(), "--reference", "mm", "--step-deg", "20",
        "--output", "amber", "openmm", "--outdir", str(outdir),
    ])
    assert rc == 0
    report = json.loads((outdir / "torsion_report.json").read_text())
    assert [r["status"] for r in report["rotors"]] == ["mlp", "mlp"]
    assert (outdir / "amber" / "molecule.prmtop").is_file()
    assert (outdir / "openmm" / "system.xml").is_file()
    ir = ir_io.from_json(outdir / "ir.json")
    assert sum(t.source == "mlp" for t in ir.proper_torsions) == 12


@pytest.mark.parametrize("rounds", ["-1", "4", "1.5"])
def test_derive_torsions_rejects_bad_refinement_budget_before_loading(rounds, tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["derive-torsions", *_ethanol_args(), "--reference", "mm",
              "--refinement-rounds", rounds, "--output", "openmm", "--outdir", str(tmp_path / "out")])
    assert exc.value.code == 2
    assert not (tmp_path / "out").exists()


@pytest.mark.external
@pytest.mark.mlp
def test_derive_torsions_mace_ethanol(tmp_path: Path) -> None:
    import json

    outdir = tmp_path / "out"
    rc = main(["derive-torsions", *_ethanol_args(), "--step-deg", "20",
               "--output", "openmm", "--outdir", str(outdir)])
    assert rc == 0
    report = json.loads((outdir / "torsion_report.json").read_text())
    assert [r["status"] for r in report["rotors"]] == ["mlp", "mlp"]


# ---- MVP-3 T36：assign-charges；derive-torsions --scan-cache / --ir-in ----

def test_assign_charges_loaded_text_file(tmp_path: Path) -> None:
    import json

    import numpy as np

    from ffcompiler.ir import io as ir_io

    q = np.linspace(-0.4, 0.4, 9)
    qfile = tmp_path / "q.txt"
    qfile.write_text("\n".join(repr(float(x)) for x in q))
    outdir = tmp_path / "out"
    rc = main(["assign-charges", *_ethanol_args(), "--method", "loaded", "--charges-file", str(qfile),
               "--output", "openmm", "--outdir", str(outdir)])
    assert rc == 0
    report = json.loads((outdir / "charge_report.json").read_text())
    ir = ir_io.from_json(outdir / "ir.json")
    assert report["method"] == "loaded" and report["source"] == {"source": str(qfile)}
    assert max(abs(a.partial_charge - x) for a, x in zip(ir.atoms, q)) <= 1e-12
    assert (outdir / "openmm" / "system.xml").is_file()


def test_derive_torsions_scan_cache_then_assign_charges_refit(tmp_path: Path) -> None:
    from ffcompiler.ir import io as ir_io

    cache = tmp_path / "scans"
    tors = tmp_path / "tors"
    # 先均摊电荷（prmtop 电荷之和与 Q 差 ~1e-9），再派生二面角，换回同一电荷时重拟合应逐项相等
    rc = main(["assign-charges", *_ethanol_args(), "--method", "loaded",
               "--output", "openmm", "--outdir", str(tmp_path / "q0")])
    assert rc == 0
    rc = main(["derive-torsions", *_ethanol_args(), "--ir-in", str(tmp_path / "q0" / "ir.json"),
               "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
               "--output", "openmm", "--outdir", str(tors)])
    assert rc == 0
    assert len(list(cache.glob("*.npz"))) == 2
    mlp_ir = ir_io.from_json(tors / "ir.json")
    rc = main(["assign-charges", *_ethanol_args(), "--ir-in", str(tors / "ir.json"), "--method", "loaded",
               "--scan-cache", str(cache), "--output", "openmm", "--outdir", str(tmp_path / "q1")])
    assert rc == 0
    assert ir_io.from_json(tmp_path / "q1" / "ir.json").proper_torsions == mlp_ir.proper_torsions


def test_assign_charges_on_mlp_ir_without_cache_fails(tmp_path: Path, capsys) -> None:
    tors = tmp_path / "tors"
    assert main(["derive-torsions", *_ethanol_args(), "--reference", "mm", "--step-deg", "20",
                 "--output", "openmm", "--outdir", str(tors)]) == 0
    assert len(list((tors / "scans").glob("*.npz"))) == 2  # 默认缓存目录 <outdir>/scans
    rc = main(["assign-charges", *_ethanol_args(), "--ir-in", str(tors / "ir.json"), "--method", "loaded",
               "--output", "openmm", "--outdir", str(tmp_path / "q")])
    assert rc == 1
    assert "scan_cache" in capsys.readouterr().err


def test_assign_charges_ir_in_must_match_prmtop(tmp_path: Path, capsys) -> None:
    # 审阅 Important 2：--ir-in 的 IR 与 --prmtop/--sdf 拓扑不一致时拒绝
    from ffcompiler.importers import from_amber
    from ffcompiler.ir import io as ir_io

    mdir = FIXTURES / "ethanol"
    ir, _ = from_amber(mdir / "ethanol.prmtop", mdir / "ethanol.inpcrd", mdir / "ethanol.sdf")
    swapped = ir.replace_terms(atoms=tuple(
        a if a.index != 3 else type(a)(**{**a.__dict__, "name": "HX"}) for a in ir.atoms))
    ir_io.to_json(swapped, tmp_path / "other.json")
    rc = main(["assign-charges", *_ethanol_args(), "--ir-in", str(tmp_path / "other.json"),
               "--method", "loaded", "--output", "openmm", "--outdir", str(tmp_path / "q")])
    assert rc == 1
    assert "--ir-in" in capsys.readouterr().err


def test_mm_reference_cache_is_frozen_for_refit_and_rescanned_for_other_charges(tmp_path: Path) -> None:
    import json

    import numpy as np

    from ffcompiler.ir import io as ir_io

    cache = tmp_path / "scans"
    assert main(["assign-charges", *_ethanol_args(), "--method", "loaded",
                 "--output", "openmm", "--outdir", str(tmp_path / "q0")]) == 0
    assert main(["derive-torsions", *_ethanol_args(), "--ir-in", str(tmp_path / "q0" / "ir.json"),
                 "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "t0")]) == 0
    q0 = [a.partial_charge for a in ir_io.from_json(tmp_path / "q0" / "ir.json").atoms]
    q2 = np.array(q0)
    q2[0] += 0.05
    q2[1] -= 0.05
    qfile = tmp_path / "q2.txt"
    qfile.write_text("\n".join(repr(float(x)) for x in q2))
    # 1) 换电荷重拟合：参考势在扫描时冻结，仍命中缓存
    assert main(["assign-charges", *_ethanol_args(), "--ir-in", str(tmp_path / "t0" / "ir.json"),
                 "--method", "loaded", "--charges-file", str(qfile), "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "r2")]) == 0
    # 2) 对另一套电荷的 IR 重新 derive：参考势不同 → 重扫，键里记新电荷
    assert main(["assign-charges", *_ethanol_args(), "--method", "loaded", "--charges-file", str(qfile),
                 "--output", "openmm", "--outdir", str(tmp_path / "q2")]) == 0
    assert main(["derive-torsions", *_ethanol_args(), "--ir-in", str(tmp_path / "q2" / "ir.json"),
                 "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "t2")]) == 0
    q2_ir = [a.partial_charge for a in ir_io.from_json(tmp_path / "q2" / "ir.json").atoms]
    for path in cache.glob("*.npz"):
        with np.load(path) as d:
            assert json.loads(str(d["key"]))["reference"] == {"kind": "mm-self", "charges": q2_ir}
    prov = ir_io.from_json(tmp_path / "t2" / "ir.json").provenance["torsion"]["reference"]
    assert prov == {"kind": "mm-self", "charges": q2_ir}
