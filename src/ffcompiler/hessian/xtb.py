"""xtb Hessian 运行器与读取器（MVP-1 计划 T13）。

命令：``xtb input.xyz --gfnff|--gfn 2 --ohess vtight --chrg <Q> --uhf 0``。
``--ohess`` 先优化到 vtight 再算 Hessian，视为驻点，gradient=None，
meta 记 ``"stationary_by": "xtb ohess vtight"``。
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from pathlib import Path
from typing import Literal

import numpy as np
from rdkit import Chem

from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.units import HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2

XTB_BINARY = Path("~/software/xtb-dist/bin/xtb").expanduser()


class XtbError(RuntimeError):
    """xtb 运行失败或其输出文件不符合格式。"""


def _element_symbol(atomic_number: int) -> str:
    return Chem.GetPeriodicTable().GetElementSymbol(int(atomic_number))


def _atomic_number(symbol: str) -> int:
    return int(Chem.GetPeriodicTable().GetAtomicNumber(symbol.capitalize()))


def read_xtb_hessian(
    hessian_file: Path,
    xyz_file: Path,
    *,
    method: str,
    version: str | None = None,
) -> CartesianHessian:
    """读取 xtb 的 ``hessian``（Ha/Bohr²，行主序 (3N)²）与 ``xtbopt.xyz``（Å）。

    hessian 文件取 ``$hessian`` 之后到 ``$end``（或文件尾）的全部浮点数，
    个数不对抛 XtbError。
    """
    hessian_file = Path(hessian_file)
    xyz_file = Path(xyz_file)

    lines = hessian_file.read_text().splitlines()
    values: list[float] = []
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("$hessian"):
            in_block = True
            continue
        if in_block and stripped.startswith("$"):
            break  # $end 或下一个块
        if in_block:
            for token in stripped.split():
                values.append(float(token.replace("D", "E").replace("d", "e")))
    if not in_block:
        raise XtbError(f"{hessian_file}: 找不到 $hessian 块")

    xyz_lines = xyz_file.read_text().splitlines()
    n_atoms = int(xyz_lines[0].split()[0])
    symbols: list[str] = []
    positions_nm: list[list[float]] = []
    for line in xyz_lines[2 : 2 + n_atoms]:
        parts = line.split()
        symbols.append(parts[0])
        positions_nm.append([float(p) for p in parts[1:4]])
    if len(symbols) != n_atoms:
        raise XtbError(f"{xyz_file}: 原子行数与声明的 {n_atoms} 不符")
    atomic_numbers = np.array([_atomic_number(s) for s in symbols], dtype=np.int64)

    n3 = 3 * n_atoms
    if len(values) != n3 * n3:
        raise XtbError(
            f"{hessian_file}: $hessian 数值个数 {len(values)} != (3N)² = {n3 * n3}"
        )
    hessian = (
        np.array(values, dtype=np.float64).reshape(n3, n3)
        * HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2
    )

    return CartesianHessian(
        atomic_numbers=atomic_numbers,
        positions_nm=np.array(positions_nm, dtype=np.float64) * 0.1,
        hessian=hessian,
        gradient=None,
        energy_kj=None,
        provider="xtb",
        method=method,
        source_file=str(hessian_file),
        source_sha256=hashlib.sha256(hessian_file.read_bytes()).hexdigest(),
        program_version=version,
        meta={"stationary_by": "xtb ohess vtight"},
    )


def _is_float(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def read_xtb_vibspectrum(path: Path, *, linear: bool = False) -> np.ndarray:
    """读取 xtb ``vibspectrum`` 的波数（cm⁻¹），去掉平移/转动的零模。

    行格式：``mode  [symmetry]  wave number  IR intensity  selection``；
    对称标签可有可无（零模行通常没有），波数 = mode 序号之后第一个数值字段。
    零模按个数去掉：|ν| 最小的 6 个（``linear=True`` 时 5 个），不用阈值，
    这样未投影干净的零模残值也不会混进来；虚频（负值）照常保留。
    返回值保持文件中的模式顺序。
    """
    wavenumbers: list[float] = []
    for lineno, line in enumerate(Path(path).read_text().splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("$end"):
            break
        if not stripped or stripped.startswith(("#", "$")):
            continue
        parts = stripped.split()
        if len(parts) < 2 or not parts[0].isdigit():
            raise ValueError(f"{path}:{lineno}: 无法解析的 vibspectrum 行 {line!r}")
        field = parts[1] if _is_float(parts[1]) else (parts[2] if len(parts) >= 3 else "")
        if not _is_float(field):
            raise ValueError(f"{path}:{lineno}: 找不到波数字段 {line!r}")
        wavenumbers.append(float(field))

    n_zero = 5 if linear else 6
    n_modes = len(wavenumbers)
    if n_modes % 3 != 0 or n_modes <= n_zero:
        raise ValueError(
            f"{path}: 模式数 {n_modes} 不是 3N（或不多于 {n_zero} 个零模）"
        )
    values = np.array(wavenumbers, dtype=np.float64)
    zero_idx = np.argsort(np.abs(values), kind="stable")[:n_zero]
    keep = np.ones(n_modes, dtype=bool)
    keep[zero_idx] = False
    return values[keep]


def _xtb_version(stdout: str) -> str | None:
    match = re.search(r"xtb\s+version\s+(\S+)", stdout)
    return match.group(1) if match else None


# 稀有气体核外电子数，valence(Z) = Z − 不超过 Z 的最大核
_NOBLE_CORES = (2, 10, 18, 36, 54, 86)


def _valence_electron_count(atomic_numbers: np.ndarray) -> int:
    """xtb 意义上的价电子总数（Z 减去稀有气体核）。

    主族与有机常用元素精确；对超出本范围体系 MVP-1 不适用。
    """
    total = 0
    for z in atomic_numbers:
        z = int(z)
        core = 0
        for c in _NOBLE_CORES:
            if z > c:
                core = c
        total += z - core
    return total


def _check_charge_and_electrons(atomic_numbers: np.ndarray, charge: int) -> None:
    """MVP-1 不支持开壳层（uhf 固定 0）：奇数价电子或负电子数在运行前拒绝。

    xtb 6.7.1 自身在 Fermi 展宽下对奇电子体系不报错，因此该校验必须由本
    运行器执行（计划 T13 测试用例 5、Review Focus 5）。
    """
    n_valence = _valence_electron_count(atomic_numbers)
    n_electrons = n_valence - int(charge)
    if n_electrons < 0:
        raise XtbError(
            f"电子数为负（价电子 {n_valence}，chrg {charge}）：--chrg 取值非法"
        )
    if n_electrons % 2 != 0:
        raise XtbError(
            f"价电子数为奇数（{n_electrons} 个，chrg {charge}）但 --uhf 0："
            "MVP-1 不支持开壳层，请检查电荷（uhf mismatch / odd electron count）"
        )


def _check_success(workdir: Path, returncode: int) -> None:
    """run_xtb_hessian 的成功条件（计划 T13）：退出码 0、CONVERGED、文件齐全。"""
    stdout_path = workdir / "xtb.out"
    stdout = stdout_path.read_text() if stdout_path.is_file() else ""

    def _fail(reason: str) -> XtbError:
        tail = "\n".join(stdout.splitlines()[-30:])
        raise XtbError(f"xtb 失败（{reason}）：\n{tail}")

    if returncode != 0:
        _fail(f"退出码 {returncode}")
    if "GEOMETRY OPTIMIZATION CONVERGED" not in stdout:
        _fail("xtb.out 没有 GEOMETRY OPTIMIZATION CONVERGED")
    if not (workdir / "hessian").is_file():
        _fail("没有生成 hessian 文件")
    if not (workdir / "xtbopt.xyz").is_file():
        _fail("没有生成 xtbopt.xyz")


def run_xtb_hessian(
    atomic_numbers: np.ndarray,
    positions_nm: np.ndarray,
    *,
    method: Literal["gfnff", "gfn2"],
    charge: int,
    workdir: Path,
    threads: int = 4,
    binary: Path = XTB_BINARY,
) -> CartesianHessian:
    """在 workdir 里跑 xtb 优化 + Hessian，返回 CartesianHessian（MVP-1 不支持开壳层）。

    成功条件：退出码 0、xtb.out 含 GEOMETRY OPTIMIZATION CONVERGED、
    hessian 与 xtbopt.xyz 存在；任何一条不满足抛 XtbError（附 xtb.out 末 30 行）。
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    _check_charge_and_electrons(
        np.asarray(atomic_numbers, dtype=np.int64), int(charge)
    )
    xyz_path = workdir / "input.xyz"
    n = len(atomic_numbers)
    xyz_lines = [str(n), "ffcompiler run_xtb_hessian"]
    for z, p in zip(atomic_numbers, positions_nm):
        xyz_lines.append(
            f"{_element_symbol(int(z))} {p[0] * 10.0:.10f} {p[1] * 10.0:.10f} {p[2] * 10.0:.10f}"
        )
    xyz_path.write_text("\n".join(xyz_lines) + "\n")

    command = [str(binary), "input.xyz"]
    if method == "gfnff":
        command.append("--gfnff")
    elif method == "gfn2":
        command += ["--gfn", "2"]
    else:
        raise ValueError(f"未知 xtb method {method!r}")
    command += ["--ohess", "vtight", "--chrg", str(int(charge)), "--uhf", "0"]

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_STACKSIZE"] = "500M"
    stdout_path = workdir / "xtb.out"
    with stdout_path.open("w") as stdout_file:
        try:
            result = subprocess.run(
                command, cwd=workdir, env=env, stdout=stdout_file,
                stderr=subprocess.STDOUT, check=False,
            )
        except OSError as exc:
            raise XtbError(f"无法启动 xtb（{binary}）：{exc}") from exc

    _check_success(workdir, result.returncode)

    return read_xtb_hessian(
        workdir / "hessian",
        workdir / "xtbopt.xyz",
        method=method,
        version=_xtb_version(stdout_path.read_text()),
    )
