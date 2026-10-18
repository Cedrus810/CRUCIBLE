"""RunContext：运行目录、带时间戳的日志与 manifest（MVP-0 计划 T0）。

每个正式运行（Gate、benchmark）在 ``jobs/<name>_<YYYYMMDD>/`` 下进行：
``run.log`` 每行带 ``[YYYY-mm-dd HH:MM:SS +elapsed]`` 前缀，开头 START 行、
结尾 END 行（atexit 兜底）；退出时写 ``manifest.json``。
不探测 Git：本 NFS 工作区禁止运行 Git，提交与 dirty 状态记 null；
源码和输入身份仍由清单记录。
"""

from __future__ import annotations

import atexit
import hashlib
import json
import logging
import os
import platform
import sys
import time
from datetime import datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Sequence

SEED_DEFAULT = 20261004

_LOG_PREFIX_RE = r"^\[\d{4}-\d\d-\d\d \d\d:\d\d:\d\d \+\d+s\] "

_PACKAGES = ("numpy", "scipy", "openmm", "rdkit")


def _package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in _PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except Exception:
            versions[name] = "n/a"
    return versions


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _ElapsedFormatter(logging.Formatter):
    """每行加 ``[YYYY-mm-dd HH:MM:SS +Ns] `` 前缀（多行消息逐行加）。"""

    def __init__(self, t0_monotonic: float) -> None:
        super().__init__()
        self._t0 = t0_monotonic

    def format(self, record: logging.LogRecord) -> str:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        elapsed = max(0, int(time.monotonic() - self._t0))
        prefix = f"[{ts} +{elapsed}s] "
        message = record.getMessage()
        if record.exc_info:
            message += "\n" + self.formatException(record.exc_info)
        return "\n".join(prefix + line for line in message.splitlines())


class RunContext:
    """上下文管理器：进入建目录配日志，退出写 manifest。

    同名同日期再次运行得到 ``_2``、``_3`` 后缀，不覆盖已有目录。
    异常退出时 ``manifest.json`` 照样写出，``exit_status`` 非 0。
    """

    def __init__(
        self,
        name: str,
        root: Path = Path("jobs"),
        seed: int = SEED_DEFAULT,
        inputs: Sequence[Path] = (),
    ) -> None:
        self.name = name
        self.root = Path(root)
        self.seed = int(seed)
        self.inputs = [Path(p) for p in inputs]
        self.dir: Path | None = None
        self.log: logging.Logger
        self.manifest: dict[str, Any] = {}
        self._t0 = time.monotonic()
        self._handler: logging.FileHandler | None = None
        self._atexit_registered = False
        self._ended = False

    # ------------------------------------------------------------------
    def _next_dir(self) -> Path:
        base = f"{self.name}_{datetime.now().strftime('%Y%m%d')}"
        candidate = self.root / base
        n = 1
        while candidate.exists():
            n += 1
            candidate = self.root / f"{base}_{n}"
        return candidate

    def __enter__(self) -> "RunContext":
        self.dir = self._next_dir()
        self.dir.mkdir(parents=True, exist_ok=False)

        self._t0 = time.monotonic()
        logger = logging.getLogger(f"ffcompiler.run.{self.name}.{id(self)}")
        logger.setLevel(logging.DEBUG)
        logger.propagate = False
        handler = logging.FileHandler(self.dir / "run.log", mode="a", encoding="utf-8")
        handler.setFormatter(_ElapsedFormatter(self._t0))
        logger.addHandler(handler)
        self._handler = handler
        self.log = logger

        self.manifest = {
            "git_commit": None,
            "git_dirty": None,
            "argv": list(sys.argv),
            "host": platform.node(),
            "env_name": os.environ.get("CONDA_DEFAULT_ENV", ""),
            "python": platform.python_version(),
            "packages": _package_versions(),
            "seed": self.seed,
            "inputs": {str(p): sha256_file(p) for p in self.inputs},
            "started": datetime.now().isoformat(timespec="seconds"),
            "ended": None,
            "wall_s": None,
            "exit_status": None,
            "source_sha256": {
                str(p.relative_to(Path(__file__).resolve().parents[2])): sha256_file(p)
                for directory in ("src", "scripts", "tests")
                for p in sorted((Path(__file__).resolve().parents[2] / directory).rglob("*.py"))
            },
        }
        self.log.info(
            "START host=%s commit=null argv=%s",
            platform.node(), " ".join(sys.argv),
        )
        if not self._atexit_registered:
            atexit.register(self._atexit_end)
            self._atexit_registered = True
        return self

    # ------------------------------------------------------------------
    def _finish(self, exit_status: int) -> None:
        if self._ended or self.dir is None:
            return
        self._ended = True
        wall_s = time.monotonic() - self._t0
        self.manifest["ended"] = datetime.now().isoformat(timespec="seconds")
        self.manifest["wall_s"] = wall_s
        self.manifest["exit_status"] = exit_status
        self.log.info("END exit=%d wall=%.1fs", exit_status, wall_s)
        with open(self.dir / "manifest.json", "w", encoding="utf-8") as fh:
            json.dump(self.manifest, fh, indent=2, ensure_ascii=False)
        if self._handler is not None:
            self.log.removeHandler(self._handler)
            self._handler.close()
            self._handler = None
        logging.getLogger(f"ffcompiler.run.{self.name}.{id(self)}").handlers.clear()

    def _atexit_end(self) -> None:
        # 进程未经 with 正常收尾（os._exit、未捕获到 __exit__ 等）时兜底写出。
        self._finish(-1)

    def finish(self, exit_status: int = 0) -> None:
        """Record an explicit command result before a normal context exit."""
        self._finish(exit_status)

    def __exit__(self, exc_type, exc, tb) -> bool:
        status = 0 if exc is None else 1
        if isinstance(exc, SystemExit):
            status = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
        self._finish(status)
        return False  # 不吞异常
