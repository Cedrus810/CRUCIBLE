"""GAFF2 typing 入口（设计 §3.1）：antechamber → parmchk2 → tleap → 现有 prmtop 导入器。

不自研 typing 引擎（决策 0024 裁定 1）：antechamber 是外部 typer，与 xTB/ORCA 作外部
Hessian provider 同一模式。开工裁定（2026-11-11 计划 §1）：本机 AmberTools 的
antechamber 没有 ``-c zero`` 电荷方法，故采用 §3.1 预留退路 ``-c bcc``——一次取得
typing + AM1-BCC，骨架电荷 provenance 记 am1bcc，charge provider 仍可覆盖。

配方与 ``scripts/make_fixtures.py``（12 固件的生成配方）逐行一致；``equalize_mol2_charges``
与 ``ambertools_version`` 自该脚本移入本库，make_fixtures 改为导入。
模块归属：§21 的 ``typing/gaff2.py``；SDF 检查合同在 ``molecule/io.py``。
"""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ffcompiler.importers.amber_prmtop import from_amber
from ffcompiler.molecule.io import MoleculeInput, read_sdf

if False:  # pragma: no cover - 仅类型标注
    from ffcompiler.ir.schema import ForceFieldIR

_LEAPRC = "leaprc.gaff2"
_CHARGE_SCALE = 18.2223  # prmtop 电荷单位（与 importers/amber_prmtop 一致）


class TypingError(RuntimeError):
    """入口配方（antechamber/parmchk2/tleap）失败；日志随消息给出。"""


def _tool(key: str) -> str:
    """工具位置走决策 0026 的统一解析（FFC_* → 配置 → PATH → 旧默认）。"""
    from ffcompiler.toolconfig import resolve

    return str(resolve(key).path)


def ambertools_version() -> str:
    """conda-meta 里的 ambertools 版本（原 scripts/make_fixtures.py 同名函数）。"""
    import os

    meta_dir = Path(os.environ.get("CONDA_PREFIX", "")) / "conda-meta"
    if meta_dir.is_dir():
        for meta in sorted(meta_dir.glob("ambertools-*.json")):
            data = json.loads(meta.read_text(encoding="utf-8"))
            return str(data.get("version", "unknown"))
    return "unknown"


def _run(cmd: list[str], cwd: Path, log_name: str) -> None:
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    except OSError as exc:  # 可执行文件不存在/不可执行
        raise TypingError(f"无法执行 {cmd[0]}：{exc}") from exc
    (cwd / log_name).write_text(
        f"$ {' '.join(cmd)}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n"
        f"--- exit {proc.returncode} ---\n",
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise TypingError(f"{cmd[0]} 退出码 {proc.returncode}；日志 {cwd / log_name}")


def equalize_mol2_charges(mol2: Path, charge: int, *, name: str = "molecule") -> str:
    """把 AM1-BCC 电荷的净残差均摊回每个原子，使 Σq 与 -nc 一致（1e-4 门）。

    antechamber -c bcc 的电荷经 4 位小数舍入后总和会残留 ~1e-3 的漂移，
    且 -eq 等价化不强制归零总电荷。这里只改 @<TRIPOS>ATOM 的电荷列，
    不影响原子类型，parmchk2/tleap 结果的其余部分不变。
    返回日志文本；失败抛 :class:`TypingError`（原 make_fixtures 版本 SystemExit）。
    """
    del name  # 与 make_fixtures 原签名兼容；错误消息自带 mol2 路径与行号
    lines = Path(mol2).read_text(encoding="utf-8").splitlines()
    atom_rows: list[int] = []
    in_atoms = False
    for n, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            continue
        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) < 9:
                raise TypingError(f"{mol2}:{n + 1}: mol2 ATOM 行列数不足（需 ≥9 列，含电荷列）")
            try:
                float(parts[8])
            except ValueError:
                raise TypingError(
                    f"{mol2}:{n + 1}: mol2 ATOM 行第 9 列（电荷）不是数字：{parts[8]!r}") from None
            atom_rows.append(n)
    if not atom_rows:
        raise TypingError(f"{mol2}：没有 @<TRIPOS>ATOM 记录")
    charges = [float(lines[n].split()[8]) for n in atom_rows]
    residual = charge - sum(charges)
    log = f"net charge target {charge}, raw sum {sum(charges):.6f}, residual {residual:.6f}\n"
    if abs(residual) > 1e-12:
        corr = residual / len(atom_rows)
        values = [q + corr for q in charges]
        # 以写出的精度为准：末原子吸收残余，使文件里的电荷和精确等于 -nc
        values[-1] += charge - sum(values)
        for n, value in zip(atom_rows, values):
            parts = lines[n].split()
            parts[8] = f"{value:.10f}"
            lines[n] = " ".join(parts)
        new_sum = sum(float(lines[n].split()[8]) for n in atom_rows)
        if abs(new_sum - charge) > 1e-8:
            raise TypingError(f"{mol2}：均摊后电荷和 {new_sum:.10f} 仍不等于 {charge}")
        log += f"adjusted by {corr:.6e} per atom; new sum {new_sum:.10f}\n"
    Path(mol2).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return log


@dataclass(frozen=True)
class EntryResult:
    """入口产物：骨架 IR（provenance 已增强）+ 坐标 + 配方记录。"""

    ir: "ForceFieldIR"
    positions_nm: np.ndarray
    workdir: Path
    artifacts: dict[str, Path]
    commands: list[list[str]]
    ambertools_version: str


def _file_stem(name: str) -> str:
    """文件名安全化（SDF _Name 可能带空格/斜杠）；IR 名仍用原名。"""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name) or "molecule"


def _enrich_entry_provenance(ir: "ForceFieldIR", mi: MoleculeInput,
                             commands: list[list[str]]) -> None:
    """§3.1 provenance：配方、版本、命令行、SDF hash、几何检查门（各项 source 不动）。"""
    ir.provenance["entry"] = {
        "recipe": ("antechamber -c bcc -at gaff2 → equalize → parmchk2 -s gaff2 → "
                   "tleap(saveamberparm) → importers.amber_prmtop.from_amber"),
        "recipe_note": ("§3.1 开工裁定（2026-11-11 计划 §1）：本 AmberTools 无 -c zero，"
                        "采用预留退路 -c bcc，一次取得 typing + AM1-BCC（决策 0024 裁定 1）"),
        "ambertools_version": ambertools_version(),
        "commands": [list(c) for c in commands],
        "leaprc": _LEAPRC,
        "sdf_sha256": mi.sha256,
        "geometry_checks": {"requires_3d": True, "clash_min_A": 0.5,
                            "no_geometry_generation": True},
        "charge_note": "骨架电荷为 antechamber AM1-BCC；charge provider 可覆盖（不强制）",
    }


def run_entry(sdf: Path, workdir: Path, *, name: str,
              net_charge: int | None = None) -> EntryResult:
    """SDF → 骨架 ForceFieldIR + 坐标（nm）。任何一步失败抛 :class:`TypingError`。"""
    sdf = Path(sdf)
    workdir = Path(workdir)
    mi = read_sdf(sdf, name=name, net_charge=net_charge)
    stem = _file_stem(mi.name)
    workdir.mkdir(parents=True, exist_ok=True)
    mol2 = workdir / f"{stem}.mol2"
    frcmod = workdir / f"{stem}.frcmod"
    prmtop = workdir / f"{stem}.prmtop"
    inpcrd = workdir / f"{stem}.inpcrd"

    commands: list[list[str]] = []
    cmd_ac = [_tool("antechamber"), "-i", str(sdf.resolve()), "-fi", "sdf",
              "-o", mol2.name, "-fo", "mol2", "-c", "bcc", "-at", "gaff2",
              "-nc", str(mi.net_charge), "-pf", "y"]
    _run(cmd_ac, workdir, "antechamber.log")
    commands.append(cmd_ac)
    if not mol2.is_file():
        raise TypingError(f"antechamber 未产出 {mol2}；见 {workdir / 'antechamber.log'}")
    log_eq = equalize_mol2_charges(mol2, mi.net_charge)
    (workdir / "charge_equalize.log").write_text(log_eq, encoding="utf-8")

    cmd_pc = [_tool("parmchk2"), "-i", mol2.name, "-f", "mol2", "-o", frcmod.name, "-s", "gaff2"]
    _run(cmd_pc, workdir, "parmchk2.log")
    commands.append(cmd_pc)
    if not frcmod.is_file():
        raise TypingError(f"parmchk2 未产出 {frcmod}；见 {workdir / 'parmchk2.log'}")

    (workdir / "tleap.in").write_text(
        f"source {_LEAPRC}\n"
        f"MOL = loadmol2 {mol2.name}\n"
        f"loadamberparams {frcmod.name}\n"
        f"saveamberparm MOL {prmtop.name} {inpcrd.name}\n"
        "quit\n",
        encoding="utf-8",
    )
    cmd_tl = [_tool("tleap"), "-f", "tleap.in"]
    _run(cmd_tl, workdir, "tleap.log")
    commands.append(cmd_tl)
    if not (prmtop.is_file() and inpcrd.is_file()):
        raise TypingError(f"tleap 未产出 prmtop/inpcrd；见 {workdir / 'tleap.log'}")
    (workdir / "leap.log").unlink(missing_ok=True)

    from ffcompiler.formats.amber import read_prmtop

    prm = read_prmtop(prmtop)
    if prm.pointer("NATOM") != mi.n_atoms:
        raise TypingError(f"prmtop 有 {prm.pointer('NATOM')} 个原子，SDF 有 {mi.n_atoms}")
    q_sum = sum(q / _CHARGE_SCALE for q in prm["CHARGE"])
    if abs(q_sum - mi.net_charge) > 1e-4:
        raise TypingError(f"prmtop 电荷和 {q_sum:.6f} ≠ 净电荷 {mi.net_charge}（1e-4 门）")

    ir, positions = from_amber(prmtop, inpcrd, sdf, name=mi.name)
    _enrich_entry_provenance(ir, mi, commands)
    artifacts = {"mol2": mol2, "frcmod": frcmod, "prmtop": prmtop, "inpcrd": inpcrd,
                 "tleap_in": workdir / "tleap.in"}
    return EntryResult(ir=ir, positions_nm=positions, workdir=workdir,
                       artifacts=artifacts, commands=commands,
                       ambertools_version=ambertools_version())
