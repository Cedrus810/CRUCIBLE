"""生成 Sobtop 参考数据（MVP-1 计划 T16）。

用法：``python scripts/make_sobtop_refs.py``（在仓库根目录，openmm_dev 环境）。

对每个 case（2 个 fchk 来源 + MVP-0 的 12 个 GAFF2 固件中除乙腈外的 11 个，
xtb GFN2 Hessian）× Sobtop 方法 2/3/4 各跑一次，param.txt 存为
``tests/fixtures/hessian/sobtop_refs/<case>/param_m{2,3,4}.txt``，
并写 ``manifest.json``（Sobtop/xtb 版本、每个 case 输入的 sha256）。

Sobtop 从仓库内 ``tools/sobtop/`` 运行（脚本自动复制并 chmod +x，
不改原始目录）；运行必须设 ``OMP_STACKSIZE=64M OMP_NUM_THREADS=4``
（2026-10-04 实测，默认线程数会报 OMP Error #34）。
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from ffcompiler.hessian.xtb import run_xtb_hessian  # noqa: E402
from ffcompiler.importers import from_amber  # noqa: E402
from ffcompiler.ir.topology import bond_graph  # noqa: E402

SOBTOP_SOURCE = Path("/home/user/software/sobtop_2026.1.16")
SOBTOP_DIR = REPO / "tools" / "sobtop"
SOBTOP_BINARY = SOBTOP_DIR / "sobtop"
SUPPORT_FILES = (
    "sobtop.ini", "assign_AT.dat", "bondcrit.dat", "bonded_param.dat",
    "LJ_param.dat", "ATOMTYPE_AMBER.DEF", "ATOMTYPE_GFF.DEF",
    "CORR_NAME_TYPE.DAT",
)
GAFF2 = REPO / "tests" / "fixtures" / "gaff2"
SOBTOP_EXAMPLES = REPO / "tests" / "fixtures" / "hessian" / "sobtop_examples"
OUT_ROOT = REPO / "tests" / "fixtures" / "hessian" / "sobtop_refs"
# 乙腈的 C-C≡N 是线性键角，Sobtop 与推导方法都跳过（计划 T16）
FCHK_CASES = {
    "sf6_fchk": (SOBTOP_EXAMPLES / "SF6.mol2", SOBTOP_EXAMPLES / "SF6.fchk"),
    "methyl_benzoate_fchk": (
        SOBTOP_EXAMPLES / "Methyl_benzoate.mol2",
        SOBTOP_EXAMPLES / "Methyl_benzoate.fchk",
    ),
}
GAFF2_NAMES = sorted(p.name for p in GAFF2.iterdir() if p.is_dir())
# 乙腈的 C-C≡N 是线性键角，Sobtop 与推导方法都跳过（计划 T16）
GAFF2_NAMES = [n for n in GAFF2_NAMES if n != "acetonitrile"]
METHODS = {2: "param_m2.txt", 3: "param_m3.txt", 4: "param_m4.txt"}


def log(message: str) -> None:
    elapsed = time.monotonic() - _START
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')} +{elapsed:7.1f}s] {message}", flush=True)


_START = time.monotonic()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_sobtop() -> None:
    """复制 Sobtop 到 tools/sobtop/（已存在则跳过），加执行位。"""
    if not SOBTOP_BINARY.is_file():
        log(f"复制 {SOBTOP_SOURCE} -> {SOBTOP_DIR}")
        shutil.copytree(SOBTOP_SOURCE, SOBTOP_DIR)
    if not (SOBTOP_SOURCE / "sobtop").is_file():
        raise SystemExit(f"原始目录里没有 sobtop 二进制：{SOBTOP_SOURCE}")
    SOBTOP_BINARY.chmod(0o755)


def run_sobtop(workdir: Path, mol2_path: Path, method: int, hessian_path: Path) -> Path:
    """在 workdir 里按 §2.8 的 stdin 序列跑 Sobtop，返回 param.txt 路径。"""
    workdir.mkdir(parents=True, exist_ok=True)
    for name in SUPPORT_FILES:
        shutil.copy2(SOBTOP_DIR / name, workdir / name)
    stdin_text = "\n".join(
        [str(mol2_path), "-1", str(method), "6", str(hessian_path), "2", "0", ""]
    )
    env = {"OMP_STACKSIZE": "64M", "OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"}
    result = subprocess.run(
        [str(SOBTOP_BINARY)], cwd=workdir, input=stdin_text, env=env,
        capture_output=True, text=True, check=False,
    )
    (workdir / "sobtop_stdout.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        raise SystemExit(
            f"Sobtop 方法 {method} 退出码 {result.returncode}：见 {workdir}/sobtop_stdout.log"
        )
    param = workdir / "param.txt"
    if not param.is_file():
        raise SystemExit(
            f"Sobtop 方法 {method} 没有生成 param.txt：见 {workdir}/sobtop_stdout.log"
        )
    return param


def counts_from_mol2(mol2_path: Path) -> tuple[int, int]:
    """mol2 键连推出的 (bond 数, angle 数)。angle：每个顶点取无序邻居对。"""
    lines = Path(mol2_path).read_text().splitlines()
    bonds: list[tuple[int, int]] = []
    in_bond_block = False
    for line in lines:
        if line.startswith("@<TRIPOS>BOND"):
            in_bond_block = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_bond_block = False
            continue
        if in_bond_block and line.strip():
            parts = line.split()
            bonds.append((int(parts[1]) - 1, int(parts[2]) - 1))
    neighbors: dict[int, set[int]] = {}
    for a, b in bonds:
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    n_angles = 0
    for vertex, nbrs in neighbors.items():
        n = len(nbrs)
        n_angles += n * (n - 1) // 2
    return len(bonds), n_angles


def write_ir_mol2(ir, positions_nm: np.ndarray, path: Path) -> None:
    """按 IR 的原子顺序和键连写 mol2，坐标来自 xtbopt（Å）。"""
    lines = [
        "@<TRIPOS>MOLECULE",
        ir.name,
        f"{len(ir.atoms)} {len(ir.bonds)}",
        "SMALL",
        "NO_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]
    for i, atom in enumerate(ir.atoms):
        x, y, z = positions_nm[i] * 10.0
        lines.append(
            f"{i + 1} {atom.name} {x:.6f} {y:.6f} {z:.6f} {atom.element} 1 MOL 0.000"
        )
    lines.append("@<TRIPOS>BOND")
    for idx, bond in enumerate(ir.bonds, start=1):
        lines.append(f"{idx} {bond.i + 1} {bond.j + 1} 1")
    path.write_text("\n".join(lines) + "\n")


def parse_sobtop_version() -> str:
    text = (REPO / "jobs" / "sobtop_version_probe" / "sobtop_stdout.log").read_text()
    match = re.search(r"Version\s*(\S+)", text, re.IGNORECASE)
    return match.group(1) if match else "unknown"


def probe_sobtop_version() -> None:
    """空跑一次 Sobtop 抓版本号。"""
    workdir = REPO / "jobs" / "sobtop_version_probe"
    workdir.mkdir(parents=True, exist_ok=True)
    for name in SUPPORT_FILES:
        shutil.copy2(SOBTOP_DIR / name, workdir / name)
    result = subprocess.run(
        [str(SOBTOP_BINARY)], cwd=workdir, input="\n", env={
            "OMP_STACKSIZE": "64M", "OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin",
        },
        capture_output=True, text=True, check=False,
    )
    (workdir / "sobtop_stdout.log").write_text(result.stdout + result.stderr)


def check_xtb_case(ir, input_positions_nm: np.ndarray, h) -> None:
    """参考数据只能来自合格的 Hessian：原子顺序与 IR 相同、优化后拓扑不变
    （Review Focus 1），且通过 §2.3 质量检查（用 IR 质量）。不合格就抛错。"""
    from ffcompiler.bonded.assign import check_hessian_topology
    from ffcompiler.hessian.base import check_hessian

    check_hessian_topology(ir, h, input_positions_nm)
    check_hessian(h, np.array([a.mass for a in ir.atoms]))


def main() -> None:
    prepare_sobtop()
    probe_sobtop_version()
    sobtop_version = parse_sobtop_version()
    log(f"Sobtop version: {sobtop_version}")

    xtb_version = "unknown"
    manifest: dict[str, object] = {"sobtop_version": sobtop_version, "cases": {}}
    if OUT_ROOT.exists():
        shutil.rmtree(OUT_ROOT)
    OUT_ROOT.mkdir(parents=True)

    cases: dict[str, tuple[Path, Path, str]] = {}
    for case, (mol2, fchk) in FCHK_CASES.items():
        cases[case] = (mol2, fchk, "fchk")
    for name in GAFF2_NAMES:
        mdir = GAFF2 / name
        ir, pos = from_amber(
            mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf"
        )
        z = np.array([a.atomic_number for a in ir.atoms])
        workdir = REPO / "jobs" / "sobtop_refs_xtb" / name
        log(f"xtb gfn2 Hessian：{name}（{len(z)} 原子，chrg {ir.total_charge}）")
        h = run_xtb_hessian(
            z, pos, method="gfn2", charge=ir.total_charge, workdir=workdir, threads=4,
        )
        check_xtb_case(ir, pos, h)
        if h.program_version:
            xtb_version = h.program_version
        case_dir = OUT_ROOT / f"{name}_gfn2"
        case_dir.mkdir(parents=True, exist_ok=True)
        mol2_path = case_dir / "input.mol2"
        # Sobtop 从 mol2 取几何：必须写 xtbopt.xyz 的坐标（计划 §2.8）
        write_ir_mol2(ir, h.positions_nm, mol2_path)
        hessian_path = case_dir / "hessian"
        shutil.copy2(workdir / "hessian", hessian_path)
        cases[f"{name}_gfn2"] = (mol2_path, hessian_path, "xtb-gfn2")

    for case, (mol2_path, hessian_path, source) in cases.items():
        case_dir = OUT_ROOT / case
        case_dir.mkdir(parents=True, exist_ok=True)
        if case in FCHK_CASES:
            shutil.copy2(mol2_path, case_dir / "input.mol2")
            shutil.copy2(hessian_path, case_dir / "input.fchk")
            hessian_path = case_dir / "input.fchk"
        n_bonds, n_angles = counts_from_mol2(case_dir / "input.mol2")
        log(f"case {case}（{source}）：mol2 推出 bond {n_bonds}、angle {n_angles}")
        for method, out_name in METHODS.items():
            workdir = REPO / "jobs" / "sobtop_refs_runs" / case / f"m{method}"
            param = run_sobtop(workdir, case_dir / "input.mol2", method, hessian_path)
            text = param.read_text()
            declared_bond = int(re.search(r"Number of bonds:\s*(\d+)", text).group(1))
            declared_angle = int(re.search(r"Number of angles:\s*(\d+)", text).group(1))
            if (declared_bond, declared_angle) != (n_bonds, n_angles):
                raise SystemExit(
                    f"{case} 方法 {method}：Sobtop 给出 bond {declared_bond}/angle "
                    f"{declared_angle}，mol2 推出 {n_bonds}/{n_angles}"
                )
            shutil.copy2(param, case_dir / out_name)
            log(f"  方法 {method} -> {out_name}")
        manifest["cases"][case] = {
            "mol2_sha256": sha256_file(case_dir / "input.mol2"),
            "hessian_sha256": sha256_file(hessian_path),
            "hessian_source": source,
        }

    manifest["xtb_version"] = xtb_version
    (OUT_ROOT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    )
    log(f"完成：{len(cases)} 个 case，manifest 写入 {OUT_ROOT / 'manifest.json'}")


if __name__ == "__main__":
    main()
