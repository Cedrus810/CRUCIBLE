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
