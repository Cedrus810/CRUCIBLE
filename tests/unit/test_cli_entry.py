"""CLI：enter 子命令与 §22 单命令形式的解析/分发（不执行配方）。"""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ETHANOL_SDF = FIXTURES / "ethanol" / "ethanol.sdf"


def test_bare_invocation_returns_2(capsys):
    from ffcompiler import cli

    assert cli.main([]) == 2


def test_single_command_dispatch(monkeypatch):
    from ffcompiler import cli

    seen = {}

    def fake_run(args):
        seen.update(vars(args))
        return 7

    monkeypatch.setattr(cli, "_do_run", fake_run)
    rc = cli.main(["fake.sdf", "--preset", "fast", "--output", "gromacs", "openmm"])
    assert rc == 7
    assert seen["sdf"] == Path("fake.sdf")
    assert seen["preset"] == "fast"
    assert seen["output"] == ["gromacs", "openmm"]


def test_subcommand_still_dispatches(monkeypatch):
    from ffcompiler import cli

    monkeypatch.setattr(cli, "_do_check_env", lambda args: 5)
    assert cli.main(["check-env"]) == 5


def test_single_command_requires_output_and_preset(capsys):
    from ffcompiler import cli

    assert cli.main(["fake.sdf", "--preset", "fast"]) == 2
    assert "--output" in capsys.readouterr().err
    assert cli.main(["fake.sdf", "--output", "gromacs"]) == 2
    assert "--preset" in capsys.readouterr().err


def test_new_charge_errors_with_guidance(capsys):
    """preset 默认 new_charge（MVP-4 未立项）必须显式报错，不静默替换。"""
    from ffcompiler import cli

    rc = cli.main(["fake.sdf", "--preset", "fast", "--output", "gromacs"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "MVP-4" in err
    assert "--charge" in err


def test_high_preset_requires_hessian_file(tmp_path, capsys):
    from ffcompiler import cli

    rc = cli.main(["fake.sdf", "--preset", "high", "--output", "gromacs",
                   "--charge", "entry", "--workdir", str(tmp_path / "w")])
    assert rc == 2
    assert "--hessian-file" in capsys.readouterr().err


def test_existing_workdir_requires_force(tmp_path, capsys):
    from ffcompiler import cli

    wd = tmp_path / "w"
    wd.mkdir()
    rc = cli.main([str(ETHANOL_SDF), "--preset", "fast", "--output", "gromacs",
                   "--charge", "entry", "--workdir", str(wd)])
    assert rc == 2
    assert "--force" in capsys.readouterr().err


def test_derive_torsions_rejects_coarse_step_deg(tmp_path, capsys):
    """A4 提前检查：点数 < 13 在扫描开始前拒绝，不必扫完一根键才发现。"""
    from ffcompiler import cli

    mdir = FIXTURES / "ethanol"
    rc = cli.main(["derive-torsions", "--prmtop", str(mdir / "ethanol.prmtop"),
                   "--inpcrd", str(mdir / "ethanol.inpcrd"), "--sdf", str(mdir / "ethanol.sdf"),
                   "--step-deg", "120", "--output", "gromacs", "--outdir", str(tmp_path)])
    assert rc == 2
    assert "13" in capsys.readouterr().err


def test_enter_round_trip(tmp_path, monkeypatch):
    """enter 子命令：monkeypatch 配方层，验证 ir/coords/manifest 落盘与合同。"""
    from ffcompiler import cli
    from ffcompiler.importers import from_amber
    from ffcompiler.typing.gaff2 import EntryResult

    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    fake = EntryResult(ir=ir, positions_nm=pos, workdir=tmp_path, artifacts={},
                       commands=[["antechamber"]], ambertools_version="test")

    def fake_run_entry(sdf, workdir, *, name, net_charge=None):
        assert name == "ethanol"
        assert net_charge is None  # None = 交给 read_sdf 从 SDF formal charge 推导（真实现）
        return fake

    monkeypatch.setattr("ffcompiler.typing.gaff2.run_entry", fake_run_entry)
    outdir = tmp_path / "out"
    assert cli.main(["enter", str(ETHANOL_SDF), "--outdir", str(outdir)]) == 0
    assert (outdir / "ir.json").is_file()
    assert (outdir / "coords.json").is_file()
    assert (outdir / "entry_manifest.json").is_file()
    from ffcompiler.ir import io as ir_io

    ir2 = ir_io.from_json(outdir / "ir.json")
    assert [a.atom_type for a in ir2.atoms] == [a.atom_type for a in ir.atoms]


def test_run_passes_every_attribute_derive_torsions_reads(tmp_path, monkeypatch):
    """`run` 构造给 derive-torsions 的 Namespace 必须含该 handler 读取的全部 args.X
    （2026-11-25：11-22 给 derive-torsions 加了 --multistart 却没同步 run，run 走到二面角步会
    AttributeError）。另核对 run 用 reference=auto（决策 0028）。"""
    import inspect
    import re
    from types import SimpleNamespace

    from ffcompiler import cli
    from ffcompiler.importers import from_amber
    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    needed = set(re.findall(r"\bargs\.(\w+)", inspect.getsource(cli._do_derive_torsions)))
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(mdir / "ethanol.prmtop", mdir / "ethanol.inpcrd", mdir / "ethanol.sdf")
    def fake_entry(sdf, workdir, **kwargs):  # 真 run_entry 会建 workdir
        Path(workdir).mkdir(parents=True, exist_ok=True)
        return SimpleNamespace(ir=ir, positions_nm=pos, ambertools_version="test")

    monkeypatch.setattr("ffcompiler.typing.gaff2.run_entry", fake_entry)

    def passthrough(args):  # derive-bonded / assign-charges 替身：原样拷贝 IR 与坐标
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        ir_io.to_json(ir_io.from_json(args.ir_in), out / "ir.json")
        ir_coords.to_json(ir_coords.from_json(args.coords_in), out / "coords.json")
        return 0

    seen = {}

    def fake_torsions(args):
        seen["args"] = args
        return passthrough(args)

    monkeypatch.setattr(cli, "_do_derive_bonded", passthrough)
    monkeypatch.setattr(cli, "_do_assign_charges", passthrough)
    monkeypatch.setattr(cli, "_do_derive_torsions", fake_torsions)
    rc = cli.main([str(ETHANOL_SDF), "--preset", "fast", "--charge", "am1bcc",
                   "--output", "openmm", "--workdir", str(tmp_path / "w")])
    assert rc == 0
    missing = sorted(a for a in needed if not hasattr(seen["args"], a))
    assert missing == [], f"run 没传这些 derive-torsions 参数：{missing}"
    assert seen["args"].reference == "auto"
