"""外部工具位置解析（决策 0026）。"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from ffcompiler import toolconfig as tc


def _exe(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\necho fake\n")
    path.chmod(0o755)
    return path


def _config(tmp_path: Path, monkeypatch, text: str) -> Path:
    cfg = tmp_path / "cfg" / "tools.toml"
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text(text, encoding="utf-8")
    monkeypatch.setenv("FFC_CONFIG", str(cfg))
    return cfg


def test_env_beats_config_beats_path(tmp_path, monkeypatch):
    env_gmx = _exe(tmp_path / "env" / "gmx")
    cfg_gmx = _exe(tmp_path / "cfgbin" / "gmx")
    path_gmx = _exe(tmp_path / "pathbin" / "gmx")
    _config(tmp_path, monkeypatch, f'[tools]\ngmx = "{cfg_gmx}"\n')
    monkeypatch.setenv("PATH", str(path_gmx.parent))
    monkeypatch.setenv("FFC_GMX", str(env_gmx))
    assert tc.resolve("gmx") == tc.Resolved("gmx", env_gmx, "env")
    monkeypatch.delenv("FFC_GMX")
    assert tc.resolve("gmx") == tc.Resolved("gmx", cfg_gmx, "config")
    _config(tmp_path, monkeypatch, "[tools]\n")
    assert tc.resolve("gmx") == tc.Resolved("gmx", path_gmx, "path")


def test_xtb_never_searches_path(tmp_path, monkeypatch):
    # MVP-1 计划 §1 固定 xtb 6.7.1：PATH 上的 xtb 不算
    _exe(tmp_path / "pathbin" / "xtb")
    monkeypatch.setenv("PATH", str(tmp_path / "pathbin"))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    monkeypatch.delenv("FFC_XTB", raising=False)
    with pytest.raises(tc.ToolNotFoundError) as info:
        tc.resolve("xtb")
    assert not any("PATH" in t for t in info.value.tried)


def test_legacy_home_default_still_works(tmp_path, monkeypatch):
    # 不写配置的老用户：~/software/xtb-dist/bin/xtb 仍然生效
    home = tmp_path / "home"
    xtb = _exe(home / "software" / "xtb-dist" / "bin" / "xtb")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("FFC_XTB", raising=False)
    assert tc.resolve("xtb") == tc.Resolved("xtb", xtb, "default")


def test_explicit_bad_path_is_error_not_fallback(tmp_path, monkeypatch):
    path_gmx = _exe(tmp_path / "pathbin" / "gmx")
    monkeypatch.setenv("PATH", str(path_gmx.parent))
    monkeypatch.delenv("FFC_GMX", raising=False)
    _config(tmp_path, monkeypatch, '[tools]\ngmx = "/does/not/exist/gmx"\n')
    with pytest.raises(tc.ToolConfigError, match=r"\[tools\]\.gmx"):
        tc.resolve("gmx")
    monkeypatch.setenv("FFC_GMX", "/does/not/exist/gmx")
    with pytest.raises(tc.ToolConfigError, match="FFC_GMX"):
        tc.resolve("gmx")


def test_not_found_lists_every_place_tried(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    monkeypatch.delenv("FFC_SANDER", raising=False)
    with pytest.raises(tc.ToolNotFoundError) as info:
        tc.resolve("sander")
    msg = str(info.value)
    assert "FFC_SANDER" in msg and "[tools].sander" in msg and "PATH" in msg
    assert "ffcompiler init-config" in msg
    assert tc.find("sander") is None


def test_config_syntax_error_names_file(tmp_path, monkeypatch):
    cfg = _config(tmp_path, monkeypatch, '[tools]\ngmx = "unterminated\n')
    with pytest.raises(tc.ToolConfigError, match=re.escape(str(cfg))):
        tc.resolve("gmx")


def test_config_relative_and_tilde_paths(tmp_path, monkeypatch):
    cfg = _config(tmp_path, monkeypatch, "")
    rel = _exe(cfg.parent / "bin" / "gmx")
    home = tmp_path / "home"
    tilde = _exe(home / "opt" / "sander")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("FFC_GMX", raising=False)
    monkeypatch.delenv("FFC_SANDER", raising=False)
    cfg.write_text('[tools]\ngmx = "bin/gmx"\nsander = "~/opt/sander"\n', encoding="utf-8")
    monkeypatch.chdir(tmp_path)  # 相对路径按配置文件目录解析，与工作目录无关
    assert tc.resolve("gmx").path == rel
    assert tc.resolve("sander").path == tilde


def test_empty_env_var_is_unset(tmp_path, monkeypatch):
    cfg_gmx = _exe(tmp_path / "cfgbin" / "gmx")
    _config(tmp_path, monkeypatch, f'[tools]\ngmx = "{cfg_gmx}"\n')
    monkeypatch.setenv("FFC_GMX", "  ")
    assert tc.resolve("gmx").source == "config"


def test_missing_ffc_config_file_is_error(tmp_path, monkeypatch):
    monkeypatch.setenv("FFC_CONFIG", str(tmp_path / "nope.toml"))
    with pytest.raises(tc.ToolConfigError, match="FFC_CONFIG"):
        tc.load_config()


def test_default_config_location_follows_xdg(tmp_path, monkeypatch):
    monkeypatch.delenv("FFC_CONFIG")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    assert tc.config_path() == tmp_path / "xdg" / "crucible" / "tools.toml"
    assert tc.load_config() == ({}, None)


def test_gmx_library_dirs_from_config(tmp_path, monkeypatch):
    lib = tmp_path / "cfg" / "hdf5" / "lib"
    lib.mkdir(parents=True)
    _config(tmp_path, monkeypatch, '[runtime]\ngmx_ld_library_path = ["hdf5/lib"]\n')
    assert tc.gmx_library_dirs() == [lib]


def test_registry_keeps_env_names_and_legacy_defaults():
    # 已有环境变量名不变；旧默认位置保留作最后兜底
    assert tc.TOOLS["xtb"].legacy_default == "~/software/xtb-dist/bin/xtb"
    assert tc.TOOLS["mace_model"].legacy_default == "~/MLP/mace/MACE-POLAR-1-M.model"
    assert tc.env_var("mace_model") == "FFC_MACE_MODEL"
    assert tc.TOOLS["xtb"].path_name is None


def test_find_gromacs_missing_returns_none(tmp_path, monkeypatch):
    from ffcompiler.validation.gromacs import find_gromacs

    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    monkeypatch.delenv("FFC_GMX", raising=False)
    assert find_gromacs() is None


def test_default_model_path_from_config(tmp_path, monkeypatch):
    from ffcompiler.torsion.mace import default_model_path

    model = tmp_path / "m" / "MACE-POLAR-1-M.model"
    model.parent.mkdir()
    model.write_bytes(b"x")
    monkeypatch.delenv("FFC_MACE_MODEL", raising=False)
    _config(tmp_path, monkeypatch, f'[tools]\nmace_model = "{model}"\n')
    assert default_model_path() == model


def test_run_xtb_missing_binary_is_xtb_error(tmp_path, monkeypatch):
    import numpy as np
    from ffcompiler.hessian.xtb import XtbError, run_xtb_hessian

    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    monkeypatch.delenv("FFC_XTB", raising=False)
    with pytest.raises(XtbError, match="找不到 xtb"):
        run_xtb_hessian(np.array([1, 1]), np.array([[0, 0, 0], [0.074, 0, 0]]),
                        method="gfn2", charge=0, workdir=tmp_path / "w")


def test_cli_subprocess_uses_ffc_config(tmp_path, monkeypatch):
    # 子进程只拿到 PATH/PYTHONPATH/FFC_CONFIG 时，也能按配置找到工具
    import subprocess
    import sys

    fake = _exe(tmp_path / "bin" / "sander")
    cfg = _config(tmp_path, monkeypatch, f'[tools]\nsander = "{fake}"\n')
    repo = Path(__file__).resolve().parents[2]
    out = subprocess.run(
        [sys.executable, "-c",
         "from ffcompiler.toolconfig import resolve; print(resolve('sander').path)"],
        env={"PATH": "/usr/bin:/bin", "PYTHONPATH": str(repo / "src"), "FFC_CONFIG": str(cfg)},
        capture_output=True, text=True, check=True,
    )
    assert out.stdout.strip() == str(fake)


def test_compare_skips_missing_sander_only_when_asked(tmp_path, monkeypatch):
    from ffcompiler.importers import from_amber
    from ffcompiler.validation.engines import compare

    fx = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2" / "ethanol"
    ir, pos = from_amber(fx / "ethanol.prmtop", fx / "ethanol.inpcrd", fx / "ethanol.sdf")
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    monkeypatch.delenv("FFC_SANDER", raising=False)
    report = compare(ir, pos, tmp_path / "w", include_gromacs=False, skip_missing_sander=True)
    assert "sander" in report.skipped
    assert not any(r["leg"].startswith("sander") for r in report.rows)
    assert report.passed                    # OpenMM 各路照常比对并通过
    with pytest.raises(tc.ToolNotFoundError):
        compare(ir, pos, tmp_path / "w2", include_gromacs=False)   # 默认（Gate 用法）仍报错


def test_init_config_writes_detected_and_placeholders(tmp_path, monkeypatch, capsys):
    from ffcompiler.cli import main

    gmx = _exe(tmp_path / "pathbin" / "gmx")
    monkeypatch.setenv("PATH", str(gmx.parent))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    for k in tc.TOOLS:
        monkeypatch.delenv(tc.env_var(k), raising=False)
    out = tmp_path / "new" / "tools.toml"
    assert main(["init-config", "--output", str(out)]) == 0
    text = out.read_text(encoding="utf-8")
    assert f'gmx = "{gmx}"' in text
    assert '# xtb = "/path/to/' in text           # 没找到 → 注释掉的占位行
    import tomllib
    assert tomllib.loads(text)["tools"] == {"gmx": str(gmx)}
    # 已存在且没加 --force → 不覆盖
    assert main(["init-config", "--output", str(out)]) == 2
    assert "--force" in capsys.readouterr().err


def test_check_env_reports_without_failing_on_missing(tmp_path, monkeypatch, capsys):
    from ffcompiler.cli import main

    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    for k in tc.TOOLS:
        monkeypatch.delenv(tc.env_var(k), raising=False)
    assert main(["check-env"]) == 0        # 缺工具只报告
    out = capsys.readouterr().out
    assert "xtb" in out and "missing" in out
    assert "| 功能 |" in out
    _config(tmp_path, monkeypatch, '[tools]\ngmx = "/does/not/exist"\n')
    assert main(["check-env"]) == 1        # 配置写错才返回 1
    assert "error" in capsys.readouterr().out


def test_example_config_parses_and_has_no_real_paths():
    import tomllib

    repo = Path(__file__).resolve().parents[2]
    text = (repo / "crucible.example.toml").read_text(encoding="utf-8")
    data = tomllib.loads(text)
    assert set(data["tools"]) <= set(tc.TOOLS)
    assert all(v.startswith("/path/to/") for v in data["tools"].values())
    assert all(p.startswith("/path/to/") for p in data["runtime"]["gmx_ld_library_path"])
    assert "/home/" not in text


def test_check_features_names_missing_tools():
    rows = [{"key": k, "status": "ok" if k != "gmx" else "missing"} for k in tc.TOOLS]
    feats = {f["feature"]: f for f in tc.check_features(rows)}
    assert feats["check-equivalence：GROMACS 一路"]["missing"] == ["gmx"]
    assert feats["check-equivalence：OpenMM 各路"]["available"]
    assert not feats["Gate G0–G3（全套验收）"]["available"]
    assert set(tc.feature_tools()) <= set(tc.TOOLS)


def test_init_config_creates_missing_ffc_config_target(tmp_path, monkeypatch):
    # FFC_CONFIG 指向尚不存在的文件：init-config 应生成它，而不是被"文件不存在"挡住
    from ffcompiler.cli import main

    target = tmp_path / "proj" / "tools.toml"
    monkeypatch.setenv("FFC_CONFIG", str(target))
    assert main(["init-config"]) == 0
    assert target.is_file()


def test_init_config_force_repairs_broken_config(tmp_path, monkeypatch):
    # 已有配置语法错误、或某项写错：--force 仍能重新生成，坏项写成注释并说明原因
    import tomllib

    from ffcompiler.cli import main

    cfg = _config(tmp_path, monkeypatch, '[tools]\ngmx = "unterminated\n')
    assert main(["init-config", "--force"]) == 0
    tomllib.loads(cfg.read_text(encoding="utf-8"))
    cfg.write_text('[tools]\nsander = "/does/not/exist"\n', encoding="utf-8")
    assert main(["init-config", "--force"]) == 0
    text = cfg.read_text(encoding="utf-8")
    assert "# sander" in text and "不可用" in text


def test_relative_ffc_config_resolves_to_absolute(tmp_path, monkeypatch):
    # 相对的 FFC_CONFIG 解析出的工具路径必须是绝对路径（外部程序在别的 cwd 下启动）
    cfgdir = tmp_path / "rel"
    sander = _exe(cfgdir / "bin" / "sander")
    (cfgdir / "tools.toml").write_text('[tools]\nsander = "bin/sander"\n', encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("FFC_CONFIG", "rel/tools.toml")
    monkeypatch.delenv("FFC_SANDER", raising=False)
    got = tc.resolve("sander").path
    assert got.is_absolute() and got == sander
    monkeypatch.setenv("FFC_SANDER", "rel/bin/sander")
    assert tc.resolve("sander").path.is_absolute()
