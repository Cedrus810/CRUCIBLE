"""ffcompiler 命令行接口（MVP-0 计划 T10）。

    ffcompiler convert --prmtop P --inpcrd C --sdf S --output gromacs amber openmm --outdir D
    ffcompiler check-equivalence --prmtop P --inpcrd C --sdf S --workdir D

退出码：0 成功 / 等价通过；1 失败或出错；2 参数缺失（argparse）。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ffcompiler.importers import from_amber


def _do_convert(args: argparse.Namespace) -> int:
    ir, positions = from_amber(args.prmtop, args.inpcrd, args.sdf)
    outdir = Path(args.outdir)
    written: dict[str, Path] = {}
    if "openmm" in args.output:
        from ffcompiler.exporters.openmm import write_openmm
        written.update(write_openmm(ir, positions, outdir / "openmm"))
    if "amber" in args.output:
        from ffcompiler.exporters.amber import write_amber
        written.update(write_amber(ir, positions, outdir / "amber"))
    if "gromacs" in args.output:
        from ffcompiler.exporters.gromacs import write_gromacs
        written.update(write_gromacs(ir, positions, outdir / "gromacs"))
    if args.ir_json:
        from ffcompiler.ir import io as ir_io
        ir_io.to_json(ir, Path(args.ir_json))
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def _do_check_equivalence(args: argparse.Namespace) -> int:
    from ffcompiler.validation.engines import compare

    ir, positions = from_amber(args.prmtop, args.inpcrd, args.sdf)
    report = compare(ir, positions, Path(args.workdir),
                     tleap_prmtop=args.prmtop, tleap_inpcrd=args.inpcrd)
    print(report.to_markdown())
    return 0 if report.passed else 1


def _do_derive_bonded(args: argparse.Namespace) -> int:
    import os

    from ffcompiler.bonded.assign import derive_bonded
    from ffcompiler.hessian.gaussian import read_fchk
    from ffcompiler.hessian.xtb import XTB_BINARY, run_xtb_hessian

    ir, positions = from_amber(args.prmtop, args.inpcrd, args.sdf)
    outdir = Path(args.outdir)
    if args.hessian_provider == "fchk":
        if not args.hessian_file:
            print("error: fchk provider 需要 --hessian-file", file=sys.stderr)
            return 2
        h = read_fchk(Path(args.hessian_file))
    else:
        # 计划 §1 固定 xtb 6.7.1；FFC_XTB 只用于显式替换（同 FFC_GMX），不查 PATH
        h = run_xtb_hessian(
            np.array([a.atomic_number for a in ir.atoms]), positions,
            method=args.hessian_provider, charge=ir.total_charge,
            workdir=outdir / "xtb",
            binary=Path(os.environ.get("FFC_XTB", str(XTB_BINARY))),
        )
    new_ir, hess_positions, report = derive_bonded(
        ir, h, method=args.method, symmetrize=not args.no_symmetrize,
        linear_angle_policy=args.linear_angles, input_positions_nm=positions,
    )
    written: dict[str, Path] = {}
    if "openmm" in args.output:
        from ffcompiler.exporters.openmm import write_openmm
        written.update(write_openmm(new_ir, hess_positions, outdir / "openmm"))
    if "amber" in args.output:
        from ffcompiler.exporters.amber import write_amber
        written.update(write_amber(new_ir, hess_positions, outdir / "amber"))
    if "gromacs" in args.output:
        from ffcompiler.exporters.gromacs import write_gromacs
        written.update(write_gromacs(new_ir, hess_positions, outdir / "gromacs"))
    print(
        f"derive-bonded: method={report.method} provider={report.hessian_provider}"
        f" 线性角保留={len(report.linear_angles_kept)}"
        f" 最大类内离散度={max(report.class_spread.values(), default=0.0):.3e}"
    )
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    from ffcompiler import __version__

    parser = argparse.ArgumentParser(
        prog="ffcompiler",
        description=(
            "引擎无关的小分子力场编译器：ForceFieldIR、GROMACS/AMBER/OpenMM 原生导出、"
            "由 Hessian 推导 bond/angle 参数"
        ),
    )
    parser.add_argument("--version", action="version", version=f"ffcompiler {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_convert = sub.add_parser("convert", help="导入 prmtop 并按 --output 写出")
    p_convert.add_argument("--prmtop", required=True, type=Path)
    p_convert.add_argument("--inpcrd", required=True, type=Path)
    p_convert.add_argument("--sdf", required=True, type=Path)
    p_convert.add_argument("--output", nargs="+", required=True,
                           choices=["gromacs", "amber", "openmm"])
    p_convert.add_argument("--outdir", required=True, type=Path)
    p_convert.add_argument("--ir-json", type=Path, default=None)

    p_check = sub.add_parser("check-equivalence", help="跨引擎能量/力等价检查")
    p_check.add_argument("--prmtop", required=True, type=Path)
    p_check.add_argument("--inpcrd", required=True, type=Path)
    p_check.add_argument("--sdf", required=True, type=Path)
    p_check.add_argument("--workdir", required=True, type=Path)

    p_derive = sub.add_parser(
        "derive-bonded", help="从 Hessian 推导 bond/angle 参数并导出（MVP-1）"
    )
    p_derive.add_argument("--prmtop", required=True, type=Path)
    p_derive.add_argument("--inpcrd", required=True, type=Path)
    p_derive.add_argument("--sdf", required=True, type=Path)
    p_derive.add_argument("--hessian-provider", required=True,
                          choices=["gfnff", "gfn2", "fchk"])
    p_derive.add_argument("--hessian-file", type=Path, default=None)
    p_derive.add_argument("--method", required=True,
                          choices=["mseminario", "m2seminario", "drih"])
    p_derive.add_argument("--no-symmetrize", action="store_true")
    p_derive.add_argument("--linear-angles", choices=["error", "library"],
                          default="error")
    p_derive.add_argument("--output", nargs="+", required=True,
                          choices=["gromacs", "amber", "openmm"])
    p_derive.add_argument("--outdir", required=True, type=Path)

    args = parser.parse_args(argv)
    try:
        if args.command == "convert":
            return _do_convert(args)
        if args.command == "check-equivalence":
            return _do_check_equivalence(args)
        if args.command == "derive-bonded":
            return _do_derive_bonded(args)
    except Exception as exc:  # noqa: BLE001 - CLI 边界：任何错误都以 1 退出
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 2  # pragma: no cover - argparse 已处理


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
