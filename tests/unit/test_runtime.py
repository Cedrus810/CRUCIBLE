"""T0：RunContext 与 run_logged.sh 的必测用例。"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from ffcompiler.runtime import RunContext, sha256_file

REPO_ROOT = Path(__file__).resolve().parents[2]
LOG_LINE_RE = re.compile(r"^\[\d{4}-\d\d-\d\d \d\d:\d\d:\d\d \+\d+s\] ")


def _read_manifest(run_dir: Path) -> dict:
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def test_run_dir_suffix_and_no_overwrite(tmp_path: Path) -> None:
    with RunContext("demo", root=tmp_path) as ctx:
        first = ctx.dir
        marker = first / "keep.txt"
        marker.write_text("first", encoding="utf-8")
    assert first.name.startswith("demo_") and first.is_dir()

    with RunContext("demo", root=tmp_path) as ctx:
        second = ctx.dir
    assert second.name == first.name + "_2"
    assert marker.read_text(encoding="utf-8") == "first"


def test_manifest_on_normal_and_exception_exit(tmp_path: Path) -> None:
    with RunContext("ok", root=tmp_path) as ctx:
        run_dir = ctx.dir
    manifest = _read_manifest(run_dir)
    assert manifest["exit_status"] == 0
    assert manifest["wall_s"] >= 0
    assert manifest["ended"] is not None
    for field in ("argv", "host", "env_name", "python", "packages", "seed", "started"):
        assert field in manifest

    with pytest.raises(RuntimeError, match="boom"):
        with RunContext("bad", root=tmp_path) as ctx:
            bad_dir = ctx.dir
            raise RuntimeError("boom")
    bad = _read_manifest(bad_dir)
    assert bad["exit_status"] != 0
    assert bad["wall_s"] >= 0
    assert "git_commit" in bad and "git_dirty" in bad


def test_explicit_failure_and_system_exit(tmp_path):
    with RunContext("failure", root=tmp_path) as ctx:
        ctx.finish(3)
    assert _read_manifest(ctx.dir)["exit_status"] == 3
    with pytest.raises(SystemExit):
        with RunContext("exit", root=tmp_path) as ctx:
            raise SystemExit(7)
    assert _read_manifest(ctx.dir)["exit_status"] == 7


def test_log_lines_all_timestamped(tmp_path: Path) -> None:
    with RunContext("logs", root=tmp_path) as ctx:
        ctx.log.info("hello")
        ctx.log.info("second line")
        run_dir = ctx.dir
    lines = (run_dir / "run.log").read_text(encoding="utf-8").splitlines()
    assert len(lines) >= 4  # START + 2 条 + END
    for line in lines:
        assert LOG_LINE_RE.match(line), f"untimestamped line: {line!r}"
    assert lines[0].startswith("[") and " START " in lines[0]
    assert " END " in lines[-1]


def test_inputs_sha256(tmp_path: Path) -> None:
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"ffcompiler fixture bytes")
    expected = hashlib.sha256(payload.read_bytes()).hexdigest()

    with RunContext("hash", root=tmp_path, inputs=[payload]) as ctx:
        run_dir = ctx.dir
    manifest = _read_manifest(run_dir)
    assert manifest["inputs"][str(payload)] == expected
    assert sha256_file(payload) == expected


def test_run_logged_script_exit_code_and_end_line() -> None:
    script = REPO_ROOT / "scripts" / "run_logged.sh"
    proc = subprocess.run(
        ["bash", str(script), "false"],
        capture_output=True, text=True, cwd=REPO_ROOT, timeout=60,
    )
    assert proc.returncode == 1
    last = proc.stdout.strip().splitlines()[-1]
    assert "END exit=1" in last
    assert LOG_LINE_RE.match(last), f"END line not timestamped: {last!r}"


def test_git_state_never_takes_the_index_lock(monkeypatch) -> None:
    """只读探测不得拿 index.lock：git status 可能因超时被杀，残留锁会挡住之后的提交。"""
    import subprocess as sp

    from ffcompiler import runtime

    calls = []

    def fake_run(args, **kwargs):
        calls.append(list(args))
        out = "true\n" if "--is-inside-work-tree" in args else ""
        return sp.CompletedProcess(args, 0, stdout=out, stderr="")

    monkeypatch.setattr(runtime.subprocess, "run", fake_run)
    runtime._git_state()
    assert calls and all("--no-optional-locks" in args for args in calls)


def test_version_is_consistent_and_marks_mvp1() -> None:
    """pyproject 与 ffcompiler.__version__ 一致；0.2.0 = MVP-0 + MVP-1（CHANGELOG）。"""
    import tomllib

    import ffcompiler

    root = Path(__file__).resolve().parents[2]
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    assert project["version"] == ffcompiler.__version__ == "0.2.0"
    assert "## [0.2.0]" in (root / "CHANGELOG.md").read_text()
