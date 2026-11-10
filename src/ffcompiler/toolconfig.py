"""外部工具位置的统一解析（决策 0026）。

每个工具的查找顺序：
1. 环境变量 ``FFC_<KEY>``（去掉首尾空白后非空才算）；
2. 配置文件的 ``[tools]`` 表：``FFC_CONFIG`` 指定的文件，否则
   ``$XDG_CONFIG_HOME/crucible/tools.toml``（缺省 ``~/.config/crucible/tools.toml``）；
3. PATH（``path_name`` 非 None 的工具；xtb 固定版本，不查 PATH）；
4. 旧的相对 HOME 默认位置（向后兼容，存在才用）。

1、2 显式给出但不可用 → ToolConfigError，不往下回退。配置文件不入库；
仓库里只有 ``crucible.example.toml`` 模板（2026-10-28 用户决定：交付时附一份 config，
使用者自己填写位置）。
"""

from __future__ import annotations

import os
import shutil
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

__all__ = [
    "TOOLS", "ToolSpec", "Resolved", "ToolConfigError", "ToolNotFoundError",
    "config_path", "load_config", "env_var", "resolve", "find", "gmx_library_dirs",
    "check_rows", "render_config", "FEATURES", "check_features", "feature_tools",
]

Kind = Literal["executable", "file", "dir"]
Source = Literal["env", "config", "path", "default"]
_KIND_ZH = {"executable": "可执行文件", "file": "文件", "dir": "目录"}


@dataclass(frozen=True)
class ToolSpec:
    key: str
    kind: Kind
    path_name: str | None      # PATH 上的程序名；None = 不查 PATH
    legacy_default: str | None  # 旧的相对 HOME 默认位置；None = 没有
    used_by: str


TOOLS: dict[str, ToolSpec] = {s.key: s for s in (
    ToolSpec("xtb", "executable", None, "~/software/xtb-dist/bin/xtb",
             "derive-bonded / global-fit 的 gfnff、gfn2；Gate G1"),
    ToolSpec("gmx", "executable", "gmx", "~/gmx26.3/bin/gmx",
             "check-equivalence 的 GROMACS 腿（可选）"),
    ToolSpec("antechamber", "executable", "antechamber", None,
             "assign-charges --method am1bcc；scripts/make_fixtures.py"),
    ToolSpec("parmchk2", "executable", "parmchk2", None, "scripts/make_fixtures.py"),
    ToolSpec("tleap", "executable", "tleap", None, "scripts/make_fixtures.py"),
    ToolSpec("sander", "executable", "sander", None, "check-equivalence 的 sander 腿"),
    ToolSpec("mace_model", "file", None, "~/MLP/mace/MACE-POLAR-1-M.model",
             "derive-torsions --reference mace；Gate G2"),
    ToolSpec("sobtop_source", "dir", None, "~/software/sobtop_2026.1.16",
             "scripts/make_sobtop_refs.py（复制到 tools/sobtop/）"),
)}

_LEGACY_GMX_LIBS = ("~/library/HDF5/lib",)


class ToolConfigError(ValueError):
    """配置有错：TOML 语法错误、显式给出的路径不可用、FFC_CONFIG 指向的文件不存在。"""


class ToolNotFoundError(ToolConfigError):
    """所有位置都找不到这个工具。"""

    def __init__(self, key: str, tried: list[str]):
        self.key = key
        self.tried = tried
        lines = "\n".join(f"  - {t}" for t in tried)
        super().__init__(
            f"找不到 {key}（用途：{TOOLS[key].used_by}）。试过：\n{lines}\n"
            f"在配置文件 [tools] 里写上 {key} = \"/path/to/...\"，"
            "或运行 `ffcompiler init-config` 生成配置，再用 `ffcompiler check-env` 检查。"
        )


@dataclass(frozen=True)
class Resolved:
    key: str
    path: Path
    source: Source


def env_var(key: str) -> str:
    return f"FFC_{key.upper()}"


def config_path() -> Path:
    explicit = os.environ.get("FFC_CONFIG", "").strip()
    if explicit:
        return Path(explicit).expanduser().absolute()
    base = os.environ.get("XDG_CONFIG_HOME", "").strip() or "~/.config"
    return (Path(base).expanduser() / "crucible" / "tools.toml").absolute()


def load_config() -> tuple[dict, Path | None]:
    path = config_path()
    if not path.is_file():
        if os.environ.get("FFC_CONFIG", "").strip():
            raise ToolConfigError(f"FFC_CONFIG 指向的配置文件不存在：{path}")
        return {}, None
    try:
        return tomllib.loads(path.read_text(encoding="utf-8")), path
    except tomllib.TOMLDecodeError as exc:
        raise ToolConfigError(f"{path}: TOML 语法错误：{exc}") from exc


def _config_path_value(raw: object, where: str, cfg: Path) -> Path:
    if not isinstance(raw, str) or not raw.strip():
        raise ToolConfigError(f"{cfg}: {where} 应为非空字符串路径，得到 {raw!r}")
    p = Path(raw.strip()).expanduser()
    return p if p.is_absolute() else cfg.parent / p


def _usable(kind: Kind, p: Path) -> bool:
    if kind == "executable":
        return p.is_file() and os.access(p, os.X_OK)
    return p.is_file() if kind == "file" else p.is_dir()


def resolve(key: str, *, use_config: bool = True) -> Resolved:
    """按查找顺序解析；``use_config=False`` 跳过配置文件层（init-config 重新生成配置时用）。

    返回的路径总是绝对路径（外部程序在各自的 workdir 下启动）。
    """
    spec = TOOLS[key]
    name = env_var(key)
    tried: list[str] = []

    raw_env = os.environ.get(name, "").strip()
    if raw_env:
        p = Path(raw_env).expanduser().absolute()
        if not _usable(spec.kind, p):
            raise ToolConfigError(f"{name}={raw_env} 不是可用的{_KIND_ZH[spec.kind]}")
        return Resolved(key, p, "env")
    tried.append(f"环境变量 {name}（未设）")

    data, cfg = load_config() if use_config else ({}, None)
    if cfg is not None and key in data.get("tools", {}):
        p = _config_path_value(data["tools"][key], f"[tools].{key}", cfg)
        if not _usable(spec.kind, p):
            raise ToolConfigError(f"{cfg}: [tools].{key} = {p} 不是可用的{_KIND_ZH[spec.kind]}")
        return Resolved(key, p, "config")
    if use_config:
        tried.append(f"配置文件 {cfg or config_path()} 的 [tools].{key}（"
                     + ("没写" if cfg is not None else "文件不存在") + "）")

    if spec.path_name:
        found = shutil.which(spec.path_name)
        if found:
            return Resolved(key, Path(found).absolute(), "path")
        tried.append(f"PATH 上的 {spec.path_name}（没找到）")

    if spec.legacy_default:
        p = Path(spec.legacy_default).expanduser()
        if _usable(spec.kind, p):
            return Resolved(key, p, "default")
        tried.append(f"默认位置 {p}（不存在）")

    raise ToolNotFoundError(key, tried)


def find(key: str, *, use_config: bool = True) -> Resolved | None:
    """可选工具：找不到返回 None；显式给错的路径照样报错。"""
    try:
        return resolve(key, use_config=use_config)
    except ToolNotFoundError:
        return None


def gmx_library_dirs() -> list[Path]:
    """gmx 运行时额外的库目录：配置 [runtime].gmx_ld_library_path，否则旧默认（存在才用）。"""
    data, cfg = load_config()
    runtime = data.get("runtime", {})
    if cfg is not None and "gmx_ld_library_path" in runtime:
        raw = runtime["gmx_ld_library_path"]
        if not isinstance(raw, list):
            raise ToolConfigError(f"{cfg}: [runtime].gmx_ld_library_path 应为字符串列表")
        return [_config_path_value(v, "[runtime].gmx_ld_library_path", cfg) for v in raw]
    return [p for p in (Path(s).expanduser() for s in _LEGACY_GMX_LIBS) if p.is_dir()]


_VERSION_ARGS = {"xtb": ["--version"], "gmx": ["--version"]}
_VERSION_MARK = {"xtb": "xtb version", "gmx": "GROMACS version"}


def _version(key: str, path: Path) -> str:
    """只对 xtb / gmx 取版本（各自 --version 里的那一行）；其余返回 "—"。"""
    import subprocess

    if key not in _VERSION_ARGS:
        return "—"
    try:
        proc = subprocess.run([str(path), *_VERSION_ARGS[key]], capture_output=True,
                              text=True, timeout=20, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"取版本失败：{exc}"
    for line in (proc.stdout + proc.stderr).splitlines():
        if _VERSION_MARK[key] in line:
            return line.strip()
    return "（输出里没有版本行）"


def check_rows() -> list[dict]:
    rows = []
    for key, spec in TOOLS.items():
        row = {"key": key, "used_by": spec.used_by, "source": "", "path": "", "version": ""}
        try:
            r = resolve(key)
            row.update(status="ok", source=r.source, path=str(r.path), version=_version(key, r.path))
        except ToolNotFoundError:
            row.update(status="missing")
        except ToolConfigError as exc:
            row.update(status="error", path=str(exc))
        rows.append(row)
    return rows


def render_config(detected: dict[str, Resolved | ToolConfigError | None]) -> str:
    """生成配置文本：找到的写成条目，找不到或不可用的写成注释掉的占位行。"""
    lines = [
        "# CRUCIBLE 外部工具位置（ffcompiler init-config 生成；本文件不入库）。",
        "# 查找顺序：环境变量 FFC_<KEY> → 本文件 → PATH → 旧默认位置。",
        "# 改完运行 `ffcompiler check-env` 检查。",
        "",
        "[tools]",
    ]
    for key, spec in TOOLS.items():
        found = detected.get(key)
        lines.append(f"# {key}：{spec.used_by}")
        if isinstance(found, ToolConfigError):
            reason = str(found).splitlines()[0]
            lines.append(f'# {key} = "/path/to/{key}"  # 不可用：{reason}')
        elif found is not None:
            lines.append(f'{key} = "{found.path}"  # 检测来源：{found.source}')
        else:
            lines.append(f'# {key} = "/path/to/{key}"  # 未找到，按需填写')
    lines += ["", "[runtime]", "# gmx 运行时额外的库目录（例如 HDF5）"]
    try:
        libs = gmx_library_dirs()
    except ToolConfigError:
        libs = []  # 配置本身有错时不沿用其中的库目录，留占位行
    if libs:
        lines.append("gmx_ld_library_path = [" + ", ".join(f'"{p}"' for p in libs) + "]")
    else:
        lines.append('# gmx_ld_library_path = ["/path/to/hdf5/lib"]')
    return "\n".join(lines) + "\n"


FEATURES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("derive-bonded / global-fit：--hessian-provider gfnff、gfn2", ("xtb",)),
    ("derive-bonded / global-fit：fchk、orca-hess、cp2k-out 文件", ()),
    ("assign-charges --method am1bcc", ("antechamber",)),
    ("assign-charges --method loaded、nagl", ()),
    ("derive-torsions --reference mace", ("mace_model",)),
    ("derive-torsions --reference mm", ()),
    ("convert（三引擎导出）", ()),
    ("check-equivalence：OpenMM 各路", ()),
    ("check-equivalence：sander 两路", ("sander",)),
    ("check-equivalence：GROMACS 一路", ("gmx",)),
    ("scripts/make_fixtures.py（重建固件）", ("antechamber", "parmchk2", "tleap")),
    ("Gate G0–G3（全套验收）", ("xtb", "gmx", "antechamber", "parmchk2", "tleap", "sander", "mace_model")),
)


def check_features(rows: list[dict]) -> list[dict]:
    ok = {r["key"] for r in rows if r["status"] == "ok"}
    return [{"feature": name, "available": set(need) <= ok,
             "missing": [k for k in need if k not in ok]} for name, need in FEATURES]


def feature_tools() -> set[str]:
    return {k for _, need in FEATURES for k in need}
