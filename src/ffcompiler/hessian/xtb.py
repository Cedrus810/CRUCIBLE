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
from dataclasses import replace
from pathlib import Path
from typing import Literal

import numpy as np
from rdkit import Chem

from ffcompiler.hessian.base import CartesianHessian
from ffcompiler.units import HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2

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


def _ohess_once(
    binary: Path,
    atomic_numbers: np.ndarray,
    positions_nm: np.ndarray,
    *,
    method: Literal["gfnff", "gfn2"],
    charge: int,
    workdir: Path,
    threads: int,
) -> CartesianHessian:
    """一次 ``xtb --ohess vtight``（2026-11-24 之前 run_xtb_hessian 的函数体，逐行未改）。"""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
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


def _xyz_energy_hartree(path: Path) -> float | None:
    """xtbopt.xyz 注释行 ``energy: E gnorm: g ...`` 的 E；格式不符返回 None（只用于记录）。"""
    try:
        match = re.search(r"energy:\s*(\S+)", Path(path).read_text().splitlines()[1])
        return float(match.group(1)) if match else None
    except (OSError, IndexError, ValueError):
        return None


def _restart_trigger(workdir: Path) -> float | None:
    """该次 ohess 是否需要沿虚频重启：返回 xtb 报告的最低频率（cm⁻¹），不需要则 None。

    只在结果本来就会被 ``check_hessian`` 拒绝时触发——最低频率低于同一报错线
    ``_IMAGINARY_ERROR_CM``（−20 cm⁻¹）且 xtb 写出了 ``xtbhess.xyz``。xtb 自己的虚频
    阈值（5 cm⁻¹）更宽；按它触发会改动现在只带警告、照常通过的分子，所以不用。
    vibspectrum 解析失败 → 不重启（保持原行为，由 check_hessian 照常判定）。
    """
    from ffcompiler.hessian.base import _IMAGINARY_ERROR_CM

    if not (workdir / "xtbhess.xyz").is_file():
        return None
    try:
        freqs = read_xtb_vibspectrum(workdir / "vibspectrum")
    except (OSError, ValueError):
        return None
    lowest = float(freqs.min()) if freqs.size else 0.0
    return lowest if lowest < _IMAGINARY_ERROR_CM else None


def run_xtb_hessian(
    atomic_numbers: np.ndarray,
    positions_nm: np.ndarray,
    *,
    method: Literal["gfnff", "gfn2"],
    charge: int,
    workdir: Path,
    threads: int = 4,
    binary: Path | None = None,
    max_imaginary_restarts: int = 2,
) -> CartesianHessian:
    """在 workdir 里跑 xtb 优化 + Hessian，返回 CartesianHessian（MVP-1 不支持开壳层）。

    成功条件：退出码 0、xtb.out 含 GEOMETRY OPTIMIZATION CONVERGED、
    hessian 与 xtbopt.xyz 存在；任何一条不满足抛 XtbError（附 xtb.out 末 30 行）。
    ``binary`` 为 None 时按 toolconfig 解析（FFC_XTB → 配置 → 旧默认，决策 0026）。

    鞍点重启（2026-11-25，用户批准；诊断见 experiments/20261125_saddle_diagnosis/）：
    优化停在鞍点时（最低频率 < −20 cm⁻¹，即 check_hessian 会拒绝的情形），从 xtb 写出的
    ``xtbhess.xyz``（沿虚频扭曲的几何）在 ``workdir/imag_restart_<k>/`` 重跑同一 ohess，
    最多 ``max_imaginary_restarts`` 次；返回最后一次的 Hessian，``meta["imaginary_restarts"]``
    逐次记录。仍有虚频则原样返回，由 check_hessian 照常拒绝（不静默放行）。
    不触发时返回值与改动前逐位相同（meta 不加键）。
    """
    if binary is None:
        from ffcompiler.toolconfig import ToolNotFoundError, resolve

        try:
            binary = resolve("xtb").path
        except ToolNotFoundError as exc:
            raise XtbError(str(exc)) from exc
    binary = Path(binary)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    _check_charge_and_electrons(
        np.asarray(atomic_numbers, dtype=np.int64), int(charge)
    )
    h = _ohess_once(binary, atomic_numbers, positions_nm, method=method, charge=charge,
                    workdir=workdir, threads=threads)
    restarts: list[dict] = []
    current = workdir
    for k in range(1, int(max_imaginary_restarts) + 1):
        lowest = _restart_trigger(current)
        if lowest is None:
            break
        nxt = workdir / f"imag_restart_{k}"
        _, start_pos, _ = _read_xyz_angstrom(current / "xtbhess.xyz")
        h = _ohess_once(binary, atomic_numbers, start_pos * 0.1, method=method,
                        charge=charge, workdir=nxt, threads=threads)
        restarts.append({
            "round": k,
            "from": str(current / "xtbhess.xyz"),
            "lowest_freq_before_cm": lowest,
            "energy_before_hartree": _xyz_energy_hartree(current / "xtbopt.xyz"),
            "energy_after_hartree": _xyz_energy_hartree(nxt / "xtbopt.xyz"),
        })
        current = nxt
    if restarts:
        h = replace(h, meta={**dict(h.meta), "imaginary_restarts": restarts})
    return h


def _read_xyz_angstrom(path: Path) -> tuple[list[str], np.ndarray, str]:
    """xyz（Å）→ (元素, 坐标 Å, 注释行)。"""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    if len(symbols) != n:
        raise XtbError(f"{path}: 原子行数与声明的 {n} 不符")
    return symbols, np.asarray(coords, dtype=np.float64), lines[1]


def read_xtb_gradient(path: Path, n_atoms: int) -> tuple[float, np.ndarray]:
    """读 xtb ``--grad`` 写出的 Turbomole 格式 ``gradient`` 文件，取最后一个 cycle。

    结构：``$grad`` / 每 cycle 一行 ``cycle = k  SCF energy = E  |dE/dxyz| = g`` +
    n 行坐标（Bohr，带元素）+ n 行梯度（Eh/Bohr）/ ``$end``。同一目录重复运行时 xtb
    会追加 cycle，故只取最后一个。返回 ``(energy_hartree, gradient[n, 3])``，单位
    Eh 与 Eh/Bohr（不换算，换算留给调用方，免得单位在读取层被悄悄改掉）。
    """
    lines = Path(path).read_text().splitlines()
    cycle_idx = [i for i, line in enumerate(lines) if "cycle =" in line]
    if not cycle_idx:
        raise XtbError(f"{path}: 没有 'cycle =' 行，不是 xtb gradient 文件")
    start = cycle_idx[-1]
    match = re.search(r"SCF energy\s*=\s*(\S+)", lines[start])
    if match is None:
        raise XtbError(f"{path}: cycle 行没有 SCF energy：{lines[start]!r}")
    energy = float(match.group(1).replace("D", "E"))
    body = []
    for line in lines[start + 1:]:
        if line.strip().startswith("$"):
            break
        if line.strip():
            body.append(line)
    if len(body) != 2 * n_atoms:
        raise XtbError(
            f"{path}: 最后一个 cycle 有 {len(body)} 行数据，应为 2×{n_atoms}（坐标 + 梯度）")
    grad_rows = []
    for line in body[n_atoms:]:
        tokens = line.replace("D", "E").split()
        if len(tokens) != 3:
            raise XtbError(f"{path}: 梯度行不是 3 列：{line!r}")
        grad_rows.append([float(t) for t in tokens])
    return energy, np.asarray(grad_rows, dtype=np.float64)


def run_xtb_singlepoint(
    atomic_numbers: np.ndarray,
    positions_nm: np.ndarray,
    *,
    method: Literal["gfnff", "gfn2"],
    charge: int,
    workdir: Path,
    threads: int = 4,
    binary: Path | None = None,
) -> dict:
    """在给定几何上跑 xtb 单点 + 梯度（``--grad``，不优化），闭壳层（uhf 0）。

    用于覆盖集诊断的参考势适用性检查（2026-11-20 计划 §2.2 第 4 项）。成功条件：
    退出码 0、输出含 ``normal termination of xtb``、``gradient`` 文件可解析且全部有限；
    否则抛 :class:`XtbError`。开壳层/非法电荷在运行前由同一电子数校验拒绝。
    返回 ``{"energy_hartree", "gradient_hartree_per_bohr", "xtb_version", "command"}``。
    """
    if binary is None:
        from ffcompiler.toolconfig import ToolNotFoundError, resolve

        try:
            binary = resolve("xtb").path
        except ToolNotFoundError as exc:
            raise XtbError(str(exc)) from exc
    binary = Path(binary)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    atomic_numbers = np.asarray(atomic_numbers, dtype=np.int64)
    _check_charge_and_electrons(atomic_numbers, int(charge))
    for stale in ("gradient", "energy", "xtbrestart"):
        (workdir / stale).unlink(missing_ok=True)   # 防止读到上一次运行追加的 cycle
    n = len(atomic_numbers)
    xyz_lines = [str(n), "ffcompiler run_xtb_singlepoint"]
    for z, p in zip(atomic_numbers, positions_nm):
        xyz_lines.append(
            f"{_element_symbol(int(z))} {p[0] * 10.0:.10f} {p[1] * 10.0:.10f} {p[2] * 10.0:.10f}"
        )
    (workdir / "input.xyz").write_text("\n".join(xyz_lines) + "\n")

    command = [str(binary), "input.xyz"]
    if method == "gfnff":
        command.append("--gfnff")
    elif method == "gfn2":
        command += ["--gfn", "2"]
    else:
        raise ValueError(f"未知 xtb method {method!r}")
    command += ["--grad", "--chrg", str(int(charge)), "--uhf", "0"]

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
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
    stdout = stdout_path.read_text()
    tail = "\n".join(stdout.splitlines()[-30:])
    if result.returncode != 0:
        raise XtbError(f"xtb 单点失败（退出码 {result.returncode}）：\n{tail}")
    if "normal termination of xtb" not in stdout:
        raise XtbError(f"xtb 单点失败（xtb.out 没有 normal termination）：\n{tail}")
    if not (workdir / "gradient").is_file():
        raise XtbError(f"xtb 单点失败（没有生成 gradient 文件）：\n{tail}")
    energy, gradient = read_xtb_gradient(workdir / "gradient", n)
    if not (np.isfinite(energy) and np.all(np.isfinite(gradient))):
        raise XtbError(f"xtb 单点能量/梯度含非有限值（E={energy}）")
    return {
        "energy_hartree": energy,
        "gradient_hartree_per_bohr": gradient,
        "xtb_version": _xtb_version(stdout),
        "command": command,
    }
