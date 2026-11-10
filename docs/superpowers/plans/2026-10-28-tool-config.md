# 外部工具位置配置（tools.toml + init-config + check-env）实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 所有外部工具（xtb、gmx、AmberTools、MACE 模型、Sobtop 源目录）的位置统一由一份不入库的配置文件给出，仓库里只放模板，使用者用 `init-config` 生成、用 `check-env` 自检。

**Architecture:** 新增 `ffcompiler/toolconfig.py`：一张工具登记表，加一个 `resolve(key)`。每个工具的查找顺序都是 环境变量 `FFC_<KEY>` → 配置文件 `[tools]` → PATH（只对允许查 PATH 的工具）→ 旧的相对 HOME 默认位置，前三处都没有、旧默认位置也不存在就报错，报错时列出试过的每一处。现有调用点（xtb、gmx、antechamber、sander、MACE、make_fixtures、make_sobtop_refs）改成调用 `resolve`，CLI 增加 `init-config`、`check-env` 两个子命令。

**Tech Stack:** Python 3.12 标准库（`tomllib`、`shutil`、`dataclasses`），不加依赖；pytest。

**Spec:** 没有对应的设计稿章节（设计稿 §22 只定义了 CLI 形态，没涉及工具位置）。依据是 2026-10-28 本对话里用户的原话，见下节。执行者需要读：本计划、`docs/decisions/0025-review-fixes-periodic-box-semantics.md`（今天的上下文）、项目 `CLAUDE.md` 的环境表。

## 用户原话（2026-10-28）

- 关于默认路径相对 HOME：「是这样 本质是为了脱敏。。。。」
- 「但是 这么说 给别人交付的时候 应该给人家一个config 然后他们用户可以自己写位置啊」
- 「所以搞一个 config 就很好用啊。。。叫其他人自己写？」
- 对 init-config / check-env / 模板的方案：「对的。」

- 「然后我们允许用户缺项然后他自己可以缺功能」；对"sander 缺失时跳过、check-env 输出可用功能表"两处补充：「可以 走起」

## 补充（2026-10-28，用户确认后加入）

- **缺项只影响用到它的功能**：CLI `check-equivalence` 缺 sander 时，sander 两路记进 `report.skipped`，其余各路照常；为此 `compare()` 新增参数 `skip_missing_sander: bool = False`，默认值保持 Gate 的行为不变（Gate 是全套验收，缺 sander 仍然不通过）。→ Task 2 Step 4b。
- **`check-env` 额外输出"可用功能"表**：每个功能需要哪些工具、现在能不能用、缺的是哪个。→ Task 3 的 `FEATURES` / `check_features()`。

## 未确认（Claude 推断，执行前请用户过目）

以下是我补出来的设计细节，用户只确认了大方向，没有逐条确认：

1. **配置文件放在哪**：`FFC_CONFIG` 指定的路径；没指定时用 `$XDG_CONFIG_HOME/crucible/tools.toml`（没设 XDG_CONFIG_HOME 就用 `~/.config/crucible/tools.toml`）。**不支持放在仓库目录里**：放仓库里就得往 `.gitignore` 加一行，而我按规矩不改 `.gitignore`，一旦漏加，真实路径就会被提交进去。
2. **保留旧的相对 HOME 默认位置，作为最后一级兜底**：这样原来那套目录布局不写配置也照常能跑，10-18、10-28 两轮 Gate 的运行方式不受影响。
3. **显式给出的路径不能用时直接报错，不往下回退**：环境变量或配置文件写了一个路径，但它不存在或不可执行。这和现在 `FFC_GMX` 的行为一致，因为写了路径就说明用户想用它，悄悄换成 PATH 上的另一个版本反而更糟。
4. **xtb 不查 PATH**：MVP-1 计划 §1 固定用 xtb 6.7.1，现在的代码也是"FFC_XTB 只用于显式替换，不查 PATH"，这一点保留。gmx、antechamber、parmchk2、tleap、sander 照旧查 PATH。
5. **`check-env` 只在配置有错时返回 1**（语法错误、显式路径不可用），工具缺失只是报告出来、返回 0。原因是 `check-env` 只用来查看状态，不该挡住任何流程。
6. **运行记录里的绝对路径不在本计划范围内**：manifest 和 provenance 会记下实际用的绝对路径，入库前怎么脱敏另外决定，见"不做的事"。

## Global Constraints

- Python ≥ 3.12（`pyproject.toml` 的 `requires-python = ">=3.12"`）；只用标准库，不加依赖。
- openmm_dev 是生产环境：不装、不升级、不卸载任何包。
- ffcompiler 不装进环境，一律 `PYTHONPATH=src`。
- 已有的环境变量 `FFC_XTB`、`FFC_GMX`、`FFC_ANTECHAMBER`、`FFC_MACE_MODEL` 名字和优先级都不变（仍是最高优先级）。
- 不写配置的用户行为不变：旧默认位置仍然生效。验收时 G0–G3 的 report.md 必须与 `jobs/gate_mvp*_20261028_2/` 逐字相同。
- 入库文件（源码、测试、模板、文档）里不能出现真实的用户名或家目录，模板里的路径一律写 `/path/to/...`。要查哪些用户名，看不入库的项目 CLAUDE.md。
- 不加 hash、指纹之类的门（全局规则 §8）；配置文件不做版本号校验。
- agent 不执行任何 git 命令，最后只写提交草稿。

## Review Focus

1. **配置文件 TOML 语法错误**：报错要指出是哪个文件、哪一行，不能抛出一个不知所云的 `TOMLDecodeError`，也不能静默忽略。→ Task 1 `test_config_syntax_error_names_file`。
2. **配置文件里写 `~` 或相对路径**：`~` 要展开；相对路径按配置文件所在目录解析，不按当前工作目录（否则换个目录运行，结果就不一样）。→ Task 1 `test_config_relative_and_tilde_paths`。
3. **环境变量设成了空字符串**（例如 `FFC_XTB=`）：当作没设处理。现在 `Path("")` 会变成 `.`，报出莫名其妙的错。→ Task 1 `test_empty_env_var_is_unset`。
4. **可选工具没装**（例如没有 gmx）：`check-equivalence` 照旧跳过 GROMACS 腿，`check-env` 显示"未找到"，都不崩。→ Task 2 `test_find_gromacs_missing_returns_none`。
5. **CLI 子进程拿不到完整环境**（10-28 实际踩过：测试只传了 PATH 和 PYTHONPATH）：只要能拿到 `FFC_CONFIG` 或 HOME，就能找到配置文件。→ Task 2 `test_cli_subprocess_uses_ffc_config`。

---

## 文件结构

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/ffcompiler/toolconfig.py` | 新建 | 工具登记表、配置文件读取、`resolve`/`find`、gmx 额外库目录、生成配置文本 |
| `src/ffcompiler/hessian/xtb.py` | 改 | `run_xtb_hessian(binary=None)` 走 `resolve("xtb")`；删 `XTB_BINARY` |
| `src/ffcompiler/hessian/__init__.py` | 改 | 去掉 `XTB_BINARY` 导出 |
| `src/ffcompiler/hessian/fallback.py` | 改 | 不再引用 `XTB_BINARY` |
| `src/ffcompiler/cli.py` | 改 | 两处 `FFC_XTB` 读取删掉；新增 `init-config`、`check-env` 子命令 |
| `src/ffcompiler/validation/gromacs.py` | 改 | `find_gromacs` 与 `_environment` 走 toolconfig |
| `src/ffcompiler/validation/engines.py` | 改 | sander 走 `resolve("sander")` |
| `src/ffcompiler/torsion/mace.py` | 改 | `default_model_path` 走 `resolve("mace_model")`；删 `DEFAULT_MODEL` |
| `src/ffcompiler/charge/am1bcc.py` | 改 | `antechamber_exe` 走 `resolve("antechamber")` |
| `scripts/make_fixtures.py`、`scripts/make_sobtop_refs.py` | 改 | antechamber/parmchk2/tleap、Sobtop 源目录走 toolconfig |
| `scripts/env.sh` | 改 | conda 根目录、环境名、modulefiles 目录可由环境变量覆盖，默认值不变 |
| `crucible.example.toml` | 新建 | 入库模板，全部占位路径 |
| `tests/unit/test_toolconfig.py` | 新建 | toolconfig 单测 |
| `tests/unit/test_gromacs_adapter.py`、`test_torsion_mace.py`、`test_cli.py` | 改 | 跟着新接口改 |
| `tests/conftest.py` | 改 | 单测默认把 `FFC_CONFIG` 指到一个不存在的位置之外的空配置，防止开发者本机的配置影响单测（见 Task 1 Step 1） |
| 文档 | 改 | README×3 快速开始、决策 0026、CHANGELOG、HANDOVER、项目 CLAUDE.md |

---

### Task 1：toolconfig 核心

**Files:**
- Create: `src/ffcompiler/toolconfig.py`
- Create: `tests/unit/test_toolconfig.py`
- Modify: `tests/conftest.py`

**Interfaces:**
- Produces:
  - `TOOLS: dict[str, ToolSpec]`；`ToolSpec(key, kind, path_name, legacy_default, used_by)`，其中 `kind ∈ {"executable","file","dir"}`
  - `class ToolConfigError(ValueError)`；`class ToolNotFoundError(ToolConfigError)`，带属性 `key`、`tried: list[str]`
  - `Resolved(key: str, path: Path, source: Literal["env","config","path","default"])`
  - `config_path() -> Path`
  - `load_config() -> tuple[dict, Path | None]`（文件不存在时返回 `({}, None)`）
  - `resolve(key: str) -> Resolved`（找不到时抛 `ToolNotFoundError`，显式路径不可用时抛 `ToolConfigError`）
  - `find(key: str) -> Resolved | None`（找不到返回 None；显式路径错误照样抛）
  - `gmx_library_dirs() -> list[Path]`

- [ ] **Step 1：conftest 隔离开发者本机配置**

`tests/conftest.py` 末尾加：

```python
@pytest.fixture(autouse=True)
def _isolate_tool_config(monkeypatch, tmp_path_factory, request):
    """单测不读开发者本机的 ~/.config/crucible/tools.toml（external 测试除外：它们要真工具）。"""
    if request.node.get_closest_marker("external") is None:
        empty = tmp_path_factory.mktemp("toolcfg") / "tools.toml"
        empty.write_text("[tools]\n", encoding="utf-8")
        monkeypatch.setenv("FFC_CONFIG", str(empty))
```

（若 conftest 里还没有 `import pytest`，在文件头加上。）

- [ ] **Step 2：写失败的测试** `tests/unit/test_toolconfig.py`

```python
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
```

- [ ] **Step 3：跑测试确认失败**

Run：`PYTHONPATH=src /home/user/miniforge3/envs/openmm_dev/bin/python -m pytest tests/unit/test_toolconfig.py -q`（执行时把 python 换成本机 openmm_dev 的实际路径）
Expected：`ModuleNotFoundError: No module named 'ffcompiler.toolconfig'`

- [ ] **Step 4：实现** `src/ffcompiler/toolconfig.py`

```python
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
        return Path(explicit).expanduser()
    base = os.environ.get("XDG_CONFIG_HOME", "").strip() or "~/.config"
    return Path(base).expanduser() / "crucible" / "tools.toml"


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


def resolve(key: str) -> Resolved:
    spec = TOOLS[key]
    name = env_var(key)
    tried: list[str] = []

    raw_env = os.environ.get(name, "").strip()
    if raw_env:
        p = Path(raw_env).expanduser()
        if not _usable(spec.kind, p):
            raise ToolConfigError(f"{name}={raw_env} 不是可用的{_KIND_ZH[spec.kind]}")
        return Resolved(key, p, "env")
    tried.append(f"环境变量 {name}（未设）")

    data, cfg = load_config()
    if cfg is not None and key in data.get("tools", {}):
        p = _config_path_value(data["tools"][key], f"[tools].{key}", cfg)
        if not _usable(spec.kind, p):
            raise ToolConfigError(f"{cfg}: [tools].{key} = {p} 不是可用的{_KIND_ZH[spec.kind]}")
        return Resolved(key, p, "config")
    tried.append(f"配置文件 {cfg or config_path()} 的 [tools].{key}（"
                 + ("没写" if cfg is not None else "文件不存在") + "）")

    if spec.path_name:
        found = shutil.which(spec.path_name)
        if found:
            return Resolved(key, Path(found), "path")
        tried.append(f"PATH 上的 {spec.path_name}（没找到）")

    if spec.legacy_default:
        p = Path(spec.legacy_default).expanduser()
        if _usable(spec.kind, p):
            return Resolved(key, p, "default")
        tried.append(f"默认位置 {p}（不存在）")

    raise ToolNotFoundError(key, tried)


def find(key: str) -> Resolved | None:
    """可选工具：找不到返回 None；显式给错的路径照样报错。"""
    try:
        return resolve(key)
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
```

- [ ] **Step 5：跑测试确认通过**

Run：`PYTHONPATH=src <openmm_dev python> -m pytest tests/unit/test_toolconfig.py -q`
Expected：12 passed

---

### Task 2：现有调用点改走 toolconfig

**Files:**
- Modify: `src/ffcompiler/hessian/xtb.py:23,222,263 附近`、`src/ffcompiler/hessian/__init__.py:14,26`、`src/ffcompiler/hessian/fallback.py:139-141`
- Modify: `src/ffcompiler/cli.py:101,129,363,392`
- Modify: `src/ffcompiler/validation/gromacs.py:26-49`、`src/ffcompiler/validation/engines.py:160`
- Modify: `src/ffcompiler/torsion/mace.py:15-23`、`src/ffcompiler/charge/am1bcc.py:22-27`
- Modify: `tests/unit/test_gromacs_adapter.py:15-18`、`tests/unit/test_torsion_mace.py:21-40`、`tests/unit/test_toolconfig.py`（追加）、`tests/unit/test_assign.py:257-262`

**Interfaces:**
- Consumes: Task 1 的 `resolve`、`find`、`gmx_library_dirs`、`ToolNotFoundError`、`ToolConfigError`
- Produces: `run_xtb_hessian(..., binary: Path | None = None)`；`find_gromacs() -> str | None`（签名不变）；`default_model_path() -> Path`（签名不变，找不到仍抛 FileNotFoundError）；`antechamber_exe() -> str`（签名不变，找不到仍抛 ChargeError）

- [ ] **Step 1：写失败的测试**（追加到 `tests/unit/test_toolconfig.py`）

```python
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
```

把 `tests/unit/test_gromacs_adapter.py` 的

```python
    with pytest.raises(ValueError, match="not executable"):
```

改成

```python
    with pytest.raises(ValueError, match="FFC_GMX"):
```

把 `tests/unit/test_torsion_mace.py` 的 `test_default_model_path_is_polar_m` 换成：

```python
def test_default_model_path_is_polar_m(monkeypatch, tmp_path) -> None:
    # 不写配置、不设 FFC_MACE_MODEL 时，回落到旧的相对 HOME 默认位置
    model = tmp_path / "home" / "MLP" / "mace" / "MACE-POLAR-1-M.model"
    model.parent.mkdir(parents=True)
    model.write_bytes(b"x")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("FFC_MACE_MODEL", raising=False)
    assert default_model_path() == model
```

（同时删掉文件头的 `from ffcompiler.torsion import mace`。）

`tests/unit/test_assign.py` 第 257 行附近传给子进程的环境变量键，从 `("HOME", "FFC_XTB")` 改为 `("HOME", "FFC_XTB", "FFC_CONFIG", "XDG_CONFIG_HOME")`。

- [ ] **Step 2：跑测试确认失败**

Run：`PYTHONPATH=src <openmm_dev python> -m pytest tests/unit/test_toolconfig.py tests/unit/test_gromacs_adapter.py tests/unit/test_torsion_mace.py -q -m "not external"`
Expected：新增的 4 个测试失败（例如 `test_run_xtb_missing_binary_is_xtb_error` 报 "无法启动 xtb" 而不是 "找不到 xtb"）

- [ ] **Step 3：改 xtb**

`src/ffcompiler/hessian/xtb.py`：删掉第 23 行 `XTB_BINARY = ...`；`run_xtb_hessian` 的参数改为 `binary: Path | None = None,`，函数体开头（`workdir = Path(workdir)` 之前）加：

```python
    if binary is None:
        from ffcompiler.toolconfig import ToolNotFoundError, resolve

        try:
            binary = resolve("xtb").path
        except ToolNotFoundError as exc:
            raise XtbError(str(exc)) from exc
    binary = Path(binary)
```

`src/ffcompiler/hessian/__init__.py`：从 import 行和 `__all__` 里删掉 `XTB_BINARY`。

`src/ffcompiler/hessian/fallback.py:139-141`：

```python
    from ffcompiler.hessian.xtb import run_xtb_hessian
    ...
    binary = xtb_binary  # None → run_xtb_hessian 内部按 toolconfig 解析
```

`src/ffcompiler/cli.py` 两处（101/129、363/392）：import 改为 `from ffcompiler.hessian.xtb import run_xtb_hessian`，删掉 `binary=Path(os.environ.get("FFC_XTB", str(XTB_BINARY))),` 这一行（`FFC_XTB` 由 toolconfig 读取，优先级不变）。如果因此 `import os` 在函数里不再用到，一并删掉。

- [ ] **Step 4：改 gmx、sander**

`src/ffcompiler/validation/gromacs.py`：

```python
def find_gromacs() -> str | None:
    """FFC_GMX → 配置 [tools].gmx → PATH → 旧默认 ~/gmx26.3（决策 0026）。显式路径不可用时报错。"""
    from ffcompiler.toolconfig import find

    found = find("gmx")
    return str(found.path.resolve()) if found else None


def _environment(executable: str) -> dict[str, str]:
    from ffcompiler.toolconfig import gmx_library_dirs

    env = dict(os.environ)
    base = Path(executable).resolve().parents[1]
    roots = [base / "lib64", base / "lib", *gmx_library_dirs()]
    libraries = [str(p) for p in roots if p.is_dir()]
    if env.get("LD_LIBRARY_PATH"):
        libraries.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = ":".join(libraries)
    env["OMP_NUM_THREADS"] = "1"
    return env
```

（如果 `shutil` 不再被用到，删掉它的 import。）

`src/ffcompiler/validation/engines.py:160`：`["sander", "-O", ...]` 改为

```python
    from ffcompiler.toolconfig import resolve

    sander = str(resolve("sander").path)
    ...
        [sander, "-O", "-i", str(mdin), "-o", str(mdout),
```

- [ ] **Step 4b：CLI 的 check-equivalence 缺 sander 时跳过**

测试（追加到 `tests/unit/test_toolconfig.py`）：

```python
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
```

实现：`compare()` 签名加 `skip_missing_sander: bool = False`；在 `# ---- sander 两条腿` 之前：

```python
    from ffcompiler.toolconfig import find

    if skip_missing_sander and find("sander") is None:
        report.skipped["sander"] = "sander: not found（在配置文件 [tools].sander 里写上位置即可启用）"
    else:
        ...（原 sander 两条腿的代码整体缩进进 else）
```

`src/ffcompiler/cli.py` 的 `check-equivalence` 调用 `compare(...)` 时传 `skip_missing_sander=True`。

- [ ] **Step 5：改 MACE、antechamber**

`src/ffcompiler/torsion/mace.py`：删掉 `DEFAULT_MODEL`，`default_model_path` 改为

```python
def default_model_path() -> Path:
    """FFC_MACE_MODEL → 配置 [tools].mace_model → 旧默认 ~/MLP/mace/…（决策 0026）。"""
    from ffcompiler.toolconfig import ToolConfigError, resolve

    try:
        return resolve("mace_model").path
    except ToolConfigError as exc:
        raise FileNotFoundError(f"MACE 模型不可用：{exc}") from exc
```

（`import os` 若不再用到就删。`test_default_model_path_env_override` 仍期望 FileNotFoundError，满足。）

`src/ffcompiler/charge/am1bcc.py`：

```python
def antechamber_exe() -> str:
    """FFC_ANTECHAMBER → 配置 [tools].antechamber → PATH（决策 0026）。"""
    from ffcompiler.toolconfig import ToolConfigError, resolve

    try:
        return str(resolve("antechamber").path)
    except ToolConfigError as exc:
        raise ChargeError(str(exc)) from exc
```

（`shutil` 若不再用到就删。）

- [ ] **Step 6：跑测试确认通过**

Run：`PYTHONPATH=src <openmm_dev python> -m pytest -m "not external" -q`
Expected：全部通过（基线是 579 个，加上 Task 1 的 12 个和本任务新增的 4 个）

---

### Task 3：`init-config`、`check-env` 与入库模板

**Files:**
- Modify: `src/ffcompiler/toolconfig.py`（加 `render_config`、`check_rows`）
- Modify: `src/ffcompiler/cli.py`（两个子命令）
- Create: `crucible.example.toml`
- Test: `tests/unit/test_toolconfig.py`（追加）

**Interfaces:**
- Consumes: Task 1 的 `TOOLS`、`resolve`、`find`、`config_path`、`gmx_library_dirs`
- Produces:
  - `render_config(detected: dict[str, Resolved | None]) -> str`
  - `check_rows() -> list[dict]`，每行键为 `key, status, source, path, version, used_by`，其中 `status ∈ {"ok","missing","error"}`
  - CLI：`ffcompiler init-config [--output PATH] [--force]`，`ffcompiler check-env`

- [ ] **Step 1：写失败的测试**

```python
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
```

- [ ] **Step 2：跑测试确认失败**

Run：`PYTHONPATH=src <openmm_dev python> -m pytest tests/unit/test_toolconfig.py -q -k "init_config or check_env or example"`
Expected：FAIL（`invalid choice: 'init-config'`，模板文件不存在）

- [ ] **Step 3：toolconfig 增加两个函数**（追加到 `toolconfig.py`，并加进 `__all__`）

```python
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


def render_config(detected: dict[str, Resolved | None]) -> str:
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
        if found is not None:
            lines.append(f'{key} = "{found.path}"  # 检测来源：{found.source}')
        else:
            lines.append(f'# {key} = "/path/to/{key}"  # 未找到，按需填写')
    lines += ["", "[runtime]", "# gmx 运行时额外的库目录（例如 HDF5）"]
    libs = gmx_library_dirs()
    if libs:
        lines.append("gmx_ld_library_path = [" + ", ".join(f'"{p}"' for p in libs) + "]")
    else:
        lines.append('# gmx_ld_library_path = ["/path/to/hdf5/lib"]')
    return "\n".join(lines) + "\n"
```

注意：`render_config` 里检测来源为 `config` 的项也照样写出，这样 `--force` 重新生成时不会丢掉原来手写的条目。

- [ ] **Step 3b：可用功能表**（追加到 `toolconfig.py`，加进 `__all__`）

```python
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
```

测试：

```python
def test_check_features_names_missing_tools():
    rows = [{"key": k, "status": "ok" if k != "gmx" else "missing"} for k in tc.TOOLS]
    feats = {f["feature"]: f for f in tc.check_features(rows)}
    assert feats["check-equivalence：GROMACS 一路"]["missing"] == ["gmx"]
    assert feats["check-equivalence：OpenMM 各路"]["available"]
    assert not feats["Gate G0–G3（全套验收）"]["available"]
    assert set(tc.feature_tools()) <= set(tc.TOOLS)
```

其中 `feature_tools()` 返回 FEATURES 里出现过的全部工具键（用来保证功能表里没有写错的工具名）：

```python
def feature_tools() -> set[str]:
    return {k for _, need in FEATURES for k in need}
```

`_do_check_env` 在工具表之后再打印：

```python
    print()
    print("| 功能 | 可用 | 缺 |")
    print("|---|---|---|")
    for f in tc.check_features(rows):
        print(f"| {f['feature']} | {'yes' if f['available'] else 'no'} | {', '.join(f['missing'])} |")
```

`test_check_env_reports_without_failing_on_missing` 追加断言：`assert "| 功能 |" in out`。

- [ ] **Step 4：CLI 子命令**（`src/ffcompiler/cli.py`）

在 `main` 的子命令注册处（`sub = parser.add_subparsers(...)` 之后）加：

```python
    p_init = sub.add_parser("init-config", help="生成外部工具位置配置（不入库），能自动找到的先填上")
    p_init.add_argument("--output", type=Path, default=None,
                        help="默认写到 FFC_CONFIG 或 ~/.config/crucible/tools.toml")
    p_init.add_argument("--force", action="store_true", help="覆盖已存在的文件")
    sub.add_parser("check-env", help="逐项显示外部工具解析到哪里、来源与版本")
```

并加两个处理函数，接到现有的分派逻辑里（照 `_do_derive_bonded` 的写法）：

```python
def _do_init_config(args: argparse.Namespace) -> int:
    from ffcompiler import toolconfig as tc

    out = args.output or tc.config_path()
    if out.exists() and not args.force:
        print(f"error: {out} 已存在；要覆盖请加 --force", file=sys.stderr)
        return 2
    detected = {key: tc.find(key) for key in tc.TOOLS}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tc.render_config(detected), encoding="utf-8")
    missing = [k for k, v in detected.items() if v is None]
    print(f"config: {out}")
    if missing:
        print("未找到、需要手填（不用的功能可以不填）：" + ", ".join(missing))
    return 0


def _do_check_env(args: argparse.Namespace) -> int:
    from ffcompiler import toolconfig as tc

    rows = tc.check_rows()
    print(f"config file: {tc.config_path()}")
    print("| tool | status | source | path | version | used by |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['key']} | {r['status']} | {r['source']} | {r['path']} | {r['version']} | {r['used_by']} |")
    return 1 if any(r["status"] == "error" for r in rows) else 0
```

- [ ] **Step 5：入库模板** `crucible.example.toml`

```toml
# CRUCIBLE 外部工具位置模板。
# 复制到 ~/.config/crucible/tools.toml（或用 FFC_CONFIG 指定的位置）后改成本机路径；
# 也可以运行 `ffcompiler init-config` 自动生成，再运行 `ffcompiler check-env` 检查。
# 查找顺序：环境变量 FFC_<KEY> → 配置文件 → PATH → 旧默认位置。不用的功能可以不填。

[tools]
xtb = "/path/to/xtb-6.7.1/bin/xtb"                  # 必填（gfnff/gfn2 Hessian）；不查 PATH
antechamber = "/path/to/amber/bin/antechamber"      # AM1-BCC 电荷；在 PATH 上可省略
sander = "/path/to/amber/bin/sander"                # 跨引擎 sander 腿；在 PATH 上可省略
tleap = "/path/to/amber/bin/tleap"                  # 重建固件；在 PATH 上可省略
parmchk2 = "/path/to/amber/bin/parmchk2"            # 重建固件；在 PATH 上可省略
gmx = "/path/to/gromacs/bin/gmx"                    # 可选：跨引擎 GROMACS 腿
mace_model = "/path/to/MACE-POLAR-1-M.model"        # 可选：MLP 二面角扫描
sobtop_source = "/path/to/sobtop_2026.1.16"         # 可选：重新生成 Sobtop 参考

[runtime]
gmx_ld_library_path = ["/path/to/hdf5/lib"]         # 可选：gmx 需要的额外库目录
```

- [ ] **Step 6：跑测试确认通过**

Run：`PYTHONPATH=src <openmm_dev python> -m pytest tests/unit/test_toolconfig.py -q`
Expected：全部通过

---

### Task 4：脚本与 env.sh

**Files:**
- Modify: `scripts/make_fixtures.py:158,182` 及 parmchk2 调用处、`scripts/make_sobtop_refs.py:35`
- Modify: `scripts/env.sh`

**Interfaces:**
- Consumes: Task 1 的 `resolve`

- [ ] **Step 1：make_fixtures 里的程序名改成解析后的路径**

`scripts/make_fixtures.py` 里 `["antechamber", ...]`、`["parmchk2", ...]`、`["tleap", ...]` 的第一个元素，分别改成 `str(resolve("antechamber").path)`、`str(resolve("parmchk2").path)`、`str(resolve("tleap").path)`；文件头加 `from ffcompiler.toolconfig import resolve  # noqa: E402`（放在已有的 ffcompiler import 旁边）。

- [ ] **Step 2：make_sobtop_refs 的 Sobtop 源目录**

```python
def _sobtop_source() -> Path:
    from ffcompiler.toolconfig import resolve

    return resolve("sobtop_source").path
```

把第 35 行常量 `SOBTOP_SOURCE` 删掉，用到它的两处（`prepare_sobtop` 里的 copytree 和存在性检查）改为调用 `_sobtop_source()`。只有 `tools/sobtop/sobtop` 不存在、需要复制时才会解析，已复制过的机器不受影响。

- [ ] **Step 3：env.sh 可覆盖、默认不变**

```bash
# 激活 openmm_dev（MVP-0 计划 §1.1）
# 用法: source scripts/env.sh
# 可覆盖：FFC_MAMBA_ROOT（默认 ~/miniforge3）、FFC_ENV_NAME（默认 openmm_dev）、
#         FFC_MODULEFILES（默认 ~/modulefiles）、FFC_GMX_MODULE（默认 gromacs2026.3）
_root="${FFC_MAMBA_ROOT:-${HOME}/miniforge3}"
export MAMBA_EXE="${_root}/bin/mamba"
export MAMBA_ROOT_PREFIX="${_root}"
source "${_root}/etc/profile.d/mamba.sh"
mamba activate "${FFC_ENV_NAME:-openmm_dev}"

_mod_dir="${FFC_MODULEFILES:-${HOME}/modulefiles}"
_mod="${FFC_GMX_MODULE:-gromacs2026.3}"
if [ -r "${_mod_dir}/${_mod}" ]; then
    source /etc/profile.d/env-modules.sh
    module use "${_mod_dir}"
    module load "${_mod}"
fi
unset _root _mod_dir _mod
```

- [ ] **Step 4：冒烟**

Run（用一个工具不在旧默认位置的账号，**不改 HOME**，先写好该账号的 `~/.config/crucible/tools.toml`）：
`PYTHONPATH=src <openmm_dev python> -m ffcompiler.cli check-env`
Expected：xtb、antechamber、sander、tleap、parmchk2、gmx、mace_model 都是 ok；sobtop_source 是 missing 也可以接受（tools/sobtop/ 已经复制过了）。

---

### Task 5：文档、验证、提交草稿

**Files:**
- Modify: `README.md`、`README.zh-CN.md`、`README.ja.md`（快速开始加一节"配置外部工具"）
- Create: `docs/decisions/0026-tool-config.md`
- Modify: `docs/decisions/README.md`、`CHANGELOG.md`、`HANDOVER.md`、项目 `CLAUDE.md`（不入库）

- [ ] **Step 1：README 三语"配置外部工具"一节**，内容一致：复制模板或运行 `init-config` → 补缺项 → `check-env`；写明查找顺序与"配置文件不入库"。路径一律用 `/path/to/...`。

- [ ] **Step 2：决策 0026**：记录查找顺序、配置文件位置、不支持仓库内配置的理由、显式错误不回退、xtb 不查 PATH、旧默认位置保留；把"用户原话"一节原样搬过去，并注明以上细节是否经用户确认。

- [ ] **Step 3：项目 CLAUDE.md（不入库）**：非原作者账号的运行方式改为"写好该账号的 `~/.config/crucible/tools.toml` 后直接运行"，删掉"临时改 HOME"的说法。

- [ ] **Step 4：验证**
  - `pytest -m "not external"`：全部通过；
  - 用非原作者账号、只靠配置文件（不改 HOME、不设 `FFC_*`）跑 `pytest -m external`：全部通过（这一步约 10 分钟，按全局规则先报资源计划）；
  - G0–G3：同样只用配置文件、不改 HOME 跑一轮（约 30 分钟，先报计划，用户说"跑"再启动），四个 report.md 与 `jobs/gate_mvp*_20261028_2/` 逐字相同。

- [ ] **Step 5：入库前检查**：对本计划改动和新增的所有入库文件，用不入库的项目 CLAUDE.md 里列出的用户名做 grep，结果必须为空。

- [ ] **Step 6：提交草稿** `.commit_drafts/<日期>_tool-config.md`（commit message、文件清单、可直接执行的 `git add … && git commit …`）。

---

## 不做的事

- **运行记录里的绝对路径**：manifest、provenance、Gate 报告都会写下实际工具的绝对路径，这是可追溯性要求的。这些文件入库前怎么脱敏（提交前统一把家目录替换掉，还是记录里只写配置项名和版本），另外决定。今天待提交的那批文件也卡在这个问题上。
- 不支持放在仓库目录里的配置文件（理由见"未确认"第 1 条）。
- 不管 Python 环境本身（openmm_dev 怎么装），只管外部工具在哪；env.sh 只做到"可覆盖"。
- NAGL 模型从 `openff.nagl_models` 包里取，跟着 Python 环境走，不进配置。
- ORCA、Gaussian、CP2K 只用来读它们的输出文件，本项目不调用它们的程序，不进配置。
