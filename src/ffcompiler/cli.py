"""ffcompiler 命令行接口（MVP-0 计划 T10；C1 坐标合同 2026-10-17 扩展）。

    ffcompiler convert --prmtop P --inpcrd C --sdf S --output gromacs amber openmm --outdir D
                       [--ir-in F] [--coords-in G]     # IR+坐标 sidecar 可替代 Amber 输入
    ffcompiler check-equivalence --prmtop P --inpcrd C --sdf S --workdir D [--ir-in F] [--coords-in G]
    ffcompiler derive-bonded --prmtop P --inpcrd C --sdf S [--ir-in F] [--coords-in G]
                             --hessian-provider gfnff|gfn2|fchk|orca-hess [--hessian-file F]
                             --method mseminario|m2seminario|drih --output ... --outdir D [--ir-json J]
    ffcompiler derive-torsions --prmtop P --inpcrd C --sdf S [--ir-in F] [--coords-in G]
                               [--reference mace|mm] [--scan-cache D] --output ... --outdir D
    ffcompiler assign-charges --prmtop P --inpcrd C --sdf S [--ir-in F] [--coords-in G]
                              --method loaded|am1bcc|nagl [--charges-file F] [--scan-cache D]
                              --output ... --outdir D
    ffcompiler global-fit --prmtop P --inpcrd C --sdf S [--ir-in F] [--coords-in G]
                          --hessian-provider gfnff|gfn2|fchk|orca-hess [--hessian-file F]
                          [--prior mseminario] [--alpha 1e-3] --output ... --outdir D
    ffcompiler check-staleness --ir-in F

坐标合同（2026-10-15 计划 C1，2026-10-17 施工）：``--coords-in`` 只与
``--ir-in`` 同用，读取 ffcompiler/ir/coords.py 的 sidecar JSON，与 IR 显式核对
原子序数/元素序列（不按同名文件自动配对）；不给时坐标取 ``--inpcrd``。
derive-bonded / derive-torsions / assign-charges 成功时输出 ``ir.json``、
``coords.json``（实际使用/返回的几何与来源）与各自报告；中途出错不写这些
产物（输出目录里不会出现冒充完整结果的残留文件）。

退出码：0 成功 / 等价通过；1 失败或出错；2 参数缺失（argparse）。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ffcompiler.importers import from_amber


def _do_convert(args: argparse.Namespace) -> int:
    ir, positions, _ = _load_ir_positions(args, need_mol=False)
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

    ir, positions, _ = _load_ir_positions(args, need_mol=False)
    kwargs = {}
    if args.ir_in:
        # sander-tleap 腿比的是 tleap prmtop 的 Hamiltonian；--ir-in 的 IR 可能已换过
        # bonded/电荷/二面角，原 prmtop 不代表它，比了必然虚报失败（2026-10-17 验收）
        print("note: --ir-in 模式不跑 sander-tleap 腿（原 prmtop 不代表该 IR 的 Hamiltonian）",
              file=sys.stderr)
    elif args.prmtop and args.inpcrd:
        kwargs["tleap_prmtop"] = args.prmtop
        kwargs["tleap_inpcrd"] = args.inpcrd
    report = compare(ir, positions, Path(args.workdir), **kwargs)
    print(report.to_markdown())
    return 0 if report.passed else 1


def _do_derive_bonded(args: argparse.Namespace) -> int:
    import json
    import os

    from ffcompiler.bonded.assign import derive_bonded
    from ffcompiler.hessian.gaussian import read_fchk
    from ffcompiler.hessian.xtb import XTB_BINARY, run_xtb_hessian
    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    ir, positions, _ = _load_ir_positions(args, need_mol=False)
    outdir = Path(args.outdir)
    if args.hessian_provider == "fchk":
        if not args.hessian_file:
            print("error: fchk provider 需要 --hessian-file", file=sys.stderr)
            return 2
        h = read_fchk(Path(args.hessian_file))
    elif args.hessian_provider == "orca-hess":
        if not args.hessian_file:
            print("error: orca-hess provider 需要 --hessian-file", file=sys.stderr)
            return 2
        from ffcompiler.hessian.orca import read_orca_hess
        h = read_orca_hess(Path(args.hessian_file))
    elif args.hessian_provider == "cp2k-out":
        if not args.hessian_file:
            print("error: cp2k-out provider 需要 --hessian-file", file=sys.stderr)
            return 2
        from ffcompiler.hessian.cp2k import read_cp2k_hess
        h = read_cp2k_hess(Path(args.hessian_file), positions_nm=positions)
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
    ir_json = Path(args.ir_json) if args.ir_json else outdir / "ir.json"
    coords_json = Path(args.coords_json) if args.coords_json else outdir / "coords.json"
    ir_io.to_json(new_ir, ir_json)
    ir_coords.to_json(
        ir_coords.from_positions(
            hess_positions, new_ir, source=f"hessian:{report.hessian_provider}",
            hessian_method=report.hessian_method,
        ),
        coords_json,
    )
    (outdir / "bonded_report.json").write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"derive-bonded: method={report.method} provider={report.hessian_provider}"
        f" 线性角保留={len(report.linear_angles_kept)}"
        f" 最大类内离散度={max(report.class_spread.values(), default=0.0):.3e}"
    )
    print(f"ir: {ir_json}")
    print(f"coords: {coords_json}")
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def _do_derive_torsions(args: argparse.Namespace) -> int:
    import json

    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io
    from ffcompiler.torsion.mlp import MLPTorsionProvider

    ir, positions, mol, positions_source = _load_ir_positions(
        args, need_mol=True, with_source=True)
    outdir = Path(args.outdir)
    scan_cache = Path(args.scan_cache) if args.scan_cache else outdir / "scans"
    if args.reference == "mm":
        # 调试/测试用：参考势取输入 IR 自身的 MM，在扫描时冻结（键里记电荷）；之后换电荷重拟合仍用这次的参考能
        from ffcompiler.torsion.mm import OpenMMCalculator

        provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_),
                                      supported=frozenset(range(1, 119)),
                                      scan_kwargs={"step_deg": args.step_deg}, model_name="mm-self",
                                      scan_cache=scan_cache,
                                      refinement_rounds=args.refinement_rounds,
                                      reference={"kind": "mm-self",
                                                 "charges": [a.partial_charge for a in ir.atoms]})
    else:
        from ffcompiler.torsion.mace import default_model_path, polar_calculator

        model_path = default_model_path()
        provider = MLPTorsionProvider(lambda ir_: polar_calculator(model_path=model_path, device=args.device),
                                      scan_kwargs={"step_deg": args.step_deg}, model_path=model_path,
                                      scan_cache=scan_cache, refinement_rounds=args.refinement_rounds)
    new_ir, report = provider.parameterize(ir, positions, mol)
    outdir.mkdir(parents=True, exist_ok=True)
    written = _write_engines(new_ir, positions, outdir, args.output)
    ir_json = Path(args.ir_json) if args.ir_json else outdir / "ir.json"
    coords_json = Path(args.coords_json) if args.coords_json else outdir / "coords.json"
    ir_io.to_json(new_ir, ir_json)
    ir_coords.to_json(
        ir_coords.from_positions(positions, new_ir, source=positions_source), coords_json)
    (outdir / "torsion_report.json").write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    for row in report.rotors:
        print(f"rotor {row['bond']}: {row['status']}"
              + (f" ({row['reason']})" if row["reason"] else f" rmse={row['fit_rmse']:.3f} kJ/mol"))
    print(f"ir: {ir_json}")
    print(f"coords: {coords_json}")
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def _load_ir_positions(args: argparse.Namespace, *, need_mol: bool, with_source: bool = False):
    """统一加载（C1 坐标合同）。返回 (ir, positions, mol[, positions_source])。

    - 无 ``--ir-in``：从 Amber 三件套导入（现状行为，prmtop/inpcrd/sdf 必需，
      缺了报错并列出缺的参数）。
    - 有 ``--ir-in``：IR 取 JSON；prmtop+sdf 仍提供时做拓扑交叉检查（名字、
      元素、键集合一致）；``--coords-in`` 提供坐标 sidecar 并与 IR 显式核对，
      不给则仍取 ``--inpcrd``。``--coords-in`` 不得脱离 ``--ir-in`` 单独使用。
    """
    from rdkit import Chem

    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    if args.coords_in and not args.ir_in:
        raise ValueError("--coords-in 只能与 --ir-in 同用（坐标 sidecar 与 IR 显式配对）")

    positions_source = None
    ir = None
    positions = None

    if args.ir_in:
        ir = ir_io.from_json(args.ir_in)
        if args.prmtop and args.sdf and args.inpcrd:
            _cross_check_topology(ir, args.ir_in, args.prmtop, args.inpcrd, args.sdf)
        elif args.prmtop:
            print("note: 未做 --ir-in 与 --prmtop 的拓扑交叉检查（需要 --prmtop --inpcrd --sdf 三件套）",
                  file=sys.stderr)
        if args.coords_in:
            bundle = ir_coords.from_json(args.coords_in)
            ir_coords.check_against_ir(bundle, ir)
            positions = bundle.positions_nm
            positions_source = f"coords:{args.coords_in}"
        else:
            if not args.inpcrd:
                raise ValueError("--ir-in 模式没有 --coords-in 时需要 --inpcrd 提供坐标")
            positions = _positions_from_inpcrd(args.inpcrd)
            positions_source = f"inpcrd:{args.inpcrd}"
    else:
        missing = [
            name for name, value in (
                ("--prmtop", args.prmtop), ("--inpcrd", args.inpcrd), ("--sdf", args.sdf),
            ) if not value
        ]
        if missing:
            raise ValueError(f"缺少参数：{' '.join(missing)}（或改用 --ir-in [--coords-in]）")
        ir, positions = from_amber(args.prmtop, args.inpcrd, args.sdf)
        positions_source = f"inpcrd:{args.inpcrd}"

    mol = None
    if need_mol:
        if not args.sdf:
            raise ValueError("该命令需要 --sdf（RDKit 分子用于可转键识别）")
        mol = Chem.MolFromMolFile(str(args.sdf), removeHs=False)
    if with_source:
        return ir, positions, mol, positions_source
    return ir, positions, mol


def _positions_from_inpcrd(inpcrd: Path) -> np.ndarray:
    """--ir-in 模式下只读坐标，不导入整套 prmtop。"""
    from ffcompiler.formats.amber import read_rst7

    return read_rst7(inpcrd)


def _cross_check_topology(ir, ir_in_path: Path, prmtop: Path, inpcrd: Path, sdf: Path) -> None:
    """--ir-in 的 IR 与 prmtop/sdf 的拓扑交叉检查（沿用既有规则，不自动重排）。"""
    from ffcompiler.bonded.assign import TopologyMismatchError

    imported, _ = from_amber(prmtop, inpcrd, sdf)
    atoms = [(a.name, a.element) for a in imported.atoms]
    atoms_in = [(a.name, a.element) for a in ir.atoms]
    bonds = {frozenset((b.i, b.j)) for b in imported.bonds}
    bonds_in = {frozenset((b.i, b.j)) for b in ir.bonds}
    if atoms != atoms_in or bonds != bonds_in:
        raise TopologyMismatchError(
            f"--ir-in {ir_in_path} 与 --prmtop/--sdf 的原子（名字、元素）或键集合不一致")


def _write_engines(ir, positions, outdir: Path, outputs) -> dict[str, Path]:
    written: dict[str, Path] = {}
    if "openmm" in outputs:
        from ffcompiler.exporters.openmm import write_openmm
        written.update(write_openmm(ir, positions, outdir / "openmm"))
    if "amber" in outputs:
        from ffcompiler.exporters.amber import write_amber
        written.update(write_amber(ir, positions, outdir / "amber"))
    if "gromacs" in outputs:
        from ffcompiler.exporters.gromacs import write_gromacs
        written.update(write_gromacs(ir, positions, outdir / "gromacs"))
    return written


def _do_assign_charges(args: argparse.Namespace) -> int:
    import json

    from ffcompiler.charge.apply import assign_charges
    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    ir, positions, mol, positions_source = _load_ir_positions(
        args, need_mol=True, with_source=True)
    outdir = Path(args.outdir)
    if args.method == "loaded":
        from ffcompiler.charge.loaded import LoadedCharges
        provider = LoadedCharges(args.charges_file)
    elif args.method == "am1bcc":
        from ffcompiler.charge.am1bcc import AM1BCCCharges
        provider = AM1BCCCharges(outdir / "am1bcc")
    else:
        from ffcompiler.charge.nagl import NAGLCharges
        provider = NAGLCharges()
    new_ir, report = assign_charges(ir, provider, positions, mol, scan_cache=args.scan_cache)
    outdir.mkdir(parents=True, exist_ok=True)
    written = _write_engines(new_ir, positions, outdir, args.output)
    ir_json = Path(args.ir_json) if args.ir_json else outdir / "ir.json"
    coords_json = Path(args.coords_json) if args.coords_json else outdir / "coords.json"
    ir_io.to_json(new_ir, ir_json)
    ir_coords.to_json(
        ir_coords.from_positions(positions, new_ir, source=positions_source), coords_json)
    (outdir / "charge_report.json").write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"charges: {report.method} raw_sum={report.raw_sum:.6f} residual={report.residual:.3e}"
          + (" (MLP torsions refit)" if report.torsion_refit else ""))
    print(f"ir: {ir_json}")
    print(f"coords: {coords_json}")
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def _do_global_fit(args: argparse.Namespace) -> int:
    """global_fit 的标准导出路径（C4：先定更新规则，再暴露 CLI）。

    下界命中（n_at_bound > 0，决策 0014 的 k=0 现象）→ 报告与诊断 IR 落盘
    （``ir_failed_lower_bound.json`` / ``coords_failed_lower_bound.json``，不用
    标准名）但**不导出三后端**、退出码 1（产物状态明确：report 的 status 字段）；
    不得静默截成正值或换回 GAFF2。成功 → status=ok + 三后端 + IR + 坐标。
    """
    import hashlib
    import json

    from ffcompiler.bonded.global_fit import global_fit
    from ffcompiler.hessian.gaussian import read_fchk
    from ffcompiler.hessian.xtb import XTB_BINARY, run_xtb_hessian
    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    ir, positions, _ = _load_ir_positions(args, need_mol=False)
    outdir = Path(args.outdir)
    if args.hessian_provider == "fchk":
        if not args.hessian_file:
            print("error: fchk provider 需要 --hessian-file", file=sys.stderr)
            return 2
        h = read_fchk(Path(args.hessian_file))
    elif args.hessian_provider == "orca-hess":
        if not args.hessian_file:
            print("error: orca-hess provider 需要 --hessian-file", file=sys.stderr)
            return 2
        from ffcompiler.hessian.orca import read_orca_hess
        h = read_orca_hess(Path(args.hessian_file))
    elif args.hessian_provider == "cp2k-out":
        if not args.hessian_file:
            print("error: cp2k-out provider 需要 --hessian-file", file=sys.stderr)
            return 2
        from ffcompiler.hessian.cp2k import read_cp2k_hess
        h = read_cp2k_hess(Path(args.hessian_file), positions_nm=positions)
    else:
        import os

        h = run_xtb_hessian(
            np.array([a.atomic_number for a in ir.atoms]), positions,
            method=args.hessian_provider, charge=ir.total_charge,
            workdir=outdir / "xtb",
            binary=Path(os.environ.get("FFC_XTB", str(XTB_BINARY))),
        )
    new_ir, hess_positions, report = global_fit(
        ir, h,
        prior_method=args.prior, alpha=args.alpha,
        symmetrize=not args.no_symmetrize,
        linear_angle_policy=args.linear_angles,
        subtract_fixed=not args.no_subtract_fixed,
        input_positions_nm=positions,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "ok" if report.n_at_bound == 0 else "failed_lower_bound",
        "n_at_bound": report.n_at_bound,
        "alpha": report.alpha,
        "lam": report.lam,
        "cond": report.cond,
        "rel_residual_fit": report.rel_residual_fit,
        "rel_residual_prior": report.rel_residual_prior,
        "scale_by_class": report.scale_by_class,
        "max_force_at_x0": report.max_force_at_x0,
        "subtract_fixed": report.subtract_fixed,
        "hessian_positions_sha256": hashlib.sha256(
            np.ascontiguousarray(hess_positions).tobytes()).hexdigest(),
    }
    (outdir / "global_fit_report.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    if report.n_at_bound > 0:
        # 诊断产物不用标准名，避免下游把 k=0 的 IR 当成本次完整结果（2026-10-17 验收）
        ir_json = outdir / "ir_failed_lower_bound.json"
        coords_json = outdir / "coords_failed_lower_bound.json"
    else:
        ir_json = Path(args.ir_json) if args.ir_json else outdir / "ir.json"
        coords_json = Path(args.coords_json) if args.coords_json else outdir / "coords.json"
    ir_io.to_json(new_ir, ir_json)
    ir_coords.to_json(
        ir_coords.from_positions(
            hess_positions, new_ir, source=f"hessian:{h.provider}",
            hessian_method=h.method),
        coords_json)
    print(f"ir: {ir_json}")
    print(f"coords: {coords_json}")
    if report.n_at_bound > 0:
        print(
            f"error: global_fit 有 {report.n_at_bound} 个等价类落在下界（k=0，"
            "决策 0014 方法局限）；按产物状态约定不导出三后端。"
            "诊断见 global_fit_report.json 与 ir_failed_lower_bound.json", file=sys.stderr)
        return 1
    written = _write_engines(new_ir, hess_positions, outdir, args.output)
    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


def _do_check_staleness(args: argparse.Namespace) -> int:
    """C4：检查 IR 的拟合依赖是否过期；有过期项 → 退出码 1。"""
    from ffcompiler.bonded.staleness import staleness
    from ffcompiler.ir import io as ir_io

    ir = ir_io.from_json(args.ir_in)
    items = staleness(ir)
    if items:
        for line in items:
            print(f"stale: {line}")
        print(f"共 {len(items)} 项过期：需按依赖重拟合后再声称组合结果", file=sys.stderr)
        return 1
    print("ok: 无过期依赖")
    return 0


def _add_common_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--prmtop", type=Path, default=None)
    parser.add_argument("--inpcrd", type=Path, default=None)
    parser.add_argument("--sdf", type=Path, default=None)
    parser.add_argument("--ir-in", type=Path, default=None,
                        help="从 IR JSON 读 IR（配合 --prmtop/--inpcrd/--sdf 时做拓扑交叉检查）")
    parser.add_argument("--coords-in", type=Path, default=None,
                        help="坐标 sidecar JSON（ffcoords/0.2）；只与 --ir-in 同用")


def _validate_common_inputs(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """输入组合校验（argparse 报错 → 退出码 2，保持既有语义）。

    - ``--coords-in`` 必须与 ``--ir-in`` 同用；
    - 无 ``--ir-in``（Amber 导入模式）时 prmtop/inpcrd/sdf 三件套必需。

    只对挂了 :func:`_add_common_inputs` 的子命令做（check-staleness 只有 --ir-in）。
    """
    if not hasattr(args, "coords_in"):
        return
    if args.coords_in and not args.ir_in:
        parser.error("--coords-in 只能与 --ir-in 同用（坐标 sidecar 与 IR 显式配对）")
    if not args.ir_in:
        missing = [
            f"--{name}" for name in ("prmtop", "inpcrd", "sdf")
            if getattr(args, name) is None
        ]
        if missing:
            parser.error(f"缺少参数：{' '.join(missing)}（或改用 --ir-in [--coords-in]）")


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

    p_convert = sub.add_parser("convert", help="导入 prmtop（或 --ir-in）并按 --output 写出")
    _add_common_inputs(p_convert)
    p_convert.add_argument("--output", nargs="+", required=True,
                           choices=["gromacs", "amber", "openmm"])
    p_convert.add_argument("--outdir", required=True, type=Path)
    p_convert.add_argument("--ir-json", type=Path, default=None)

    p_check = sub.add_parser("check-equivalence", help="跨引擎能量/力等价检查")
    _add_common_inputs(p_check)
    p_check.add_argument("--workdir", required=True, type=Path)

    p_derive = sub.add_parser(
        "derive-bonded", help="从 Hessian 推导 bond/angle 参数并导出（MVP-1）"
    )
    _add_common_inputs(p_derive)
    p_derive.add_argument("--hessian-provider", required=True,
                          choices=["gfnff", "gfn2", "fchk", "orca-hess", "cp2k-out"])
    p_derive.add_argument("--hessian-file", type=Path, default=None)
    p_derive.add_argument("--method", required=True,
                          choices=["mseminario", "m2seminario", "drih", "localhess"],  # localhess experimental（决策 0023）
    )
    p_derive.add_argument("--no-symmetrize", action="store_true")
    p_derive.add_argument("--linear-angles", choices=["error", "library"],
                          default="error")
    p_derive.add_argument("--output", nargs="+", required=True,
                          choices=["gromacs", "amber", "openmm"])
    p_derive.add_argument("--outdir", required=True, type=Path)
    p_derive.add_argument("--ir-json", type=Path, default=None)
    p_derive.add_argument("--coords-json", type=Path, default=None)

    p_tors = sub.add_parser("derive-torsions", help="MLP（MACE-POLAR）二面角扫描与拟合（MVP-2）")
    _add_common_inputs(p_tors)
    p_tors.add_argument("--reference", choices=["mace", "mm"], default="mace")
    p_tors.add_argument("--device", default="cpu")
    p_tors.add_argument("--step-deg", type=float, default=10.0)
    p_tors.add_argument("--refinement-rounds", type=int, choices=range(4), default=0,
                        help="单遍之后额外做 0–3 轮顺序重拟合（默认 0，不重扫）；换电荷时沿用记录值")
    p_tors.add_argument("--output", nargs="+", required=True,
                        choices=["gromacs", "amber", "openmm"])
    p_tors.add_argument("--outdir", required=True, type=Path)
    p_tors.add_argument("--ir-json", type=Path, default=None)
    p_tors.add_argument("--coords-json", type=Path, default=None)
    p_tors.add_argument("--scan-cache", type=Path, default=None, help="扫描缓存目录（默认 <outdir>/scans）")

    p_chg = sub.add_parser("assign-charges", help="电荷 provider：载入 / AM1-BCC / NAGL（MVP-3）")
    _add_common_inputs(p_chg)
    p_chg.add_argument("--method", required=True, choices=["loaded", "am1bcc", "nagl"])
    p_chg.add_argument("--charges-file", type=Path, default=None, help="loaded：mol2 或每行一个电荷的文本")
    p_chg.add_argument("--scan-cache", type=Path, default=None, help="IR 有 MLP 二面角时重拟合用的扫描缓存")
    p_chg.add_argument("--output", nargs="+", required=True, choices=["gromacs", "amber", "openmm"])
    p_chg.add_argument("--outdir", required=True, type=Path)
    p_chg.add_argument("--ir-json", type=Path, default=None)
    p_chg.add_argument("--coords-json", type=Path, default=None)

    p_gfit = sub.add_parser(
        "global-fit", help="QMDFF 式全局 Hessian 拟合（MVP-1b；C4 起暴露 CLI）"
    )
    _add_common_inputs(p_gfit)
    p_gfit.add_argument("--hessian-provider", required=True,
                        choices=["gfnff", "gfn2", "fchk", "orca-hess", "cp2k-out"])
    p_gfit.add_argument("--hessian-file", type=Path, default=None)
    p_gfit.add_argument("--prior", choices=["mseminario", "m2seminario", "drih"],
                        default="mseminario")
    p_gfit.add_argument("--alpha", type=float, default=1e-3)
    p_gfit.add_argument("--no-symmetrize", action="store_true")
    p_gfit.add_argument("--no-subtract-fixed", action="store_true",
                        help="不扣 H_fixed（仅供诊断；标准导出必扣）")
    p_gfit.add_argument("--linear-angles", choices=["error", "library"], default="error")
    p_gfit.add_argument("--output", nargs="+", required=True,
                        choices=["gromacs", "amber", "openmm"])
    p_gfit.add_argument("--outdir", required=True, type=Path)
    p_gfit.add_argument("--ir-json", type=Path, default=None)
    p_gfit.add_argument("--coords-json", type=Path, default=None)

    p_stale = sub.add_parser(
        "check-staleness", help="检查拟合依赖（H_fixed / MM 扣除）是否过期（C4）"
    )
    p_stale.add_argument("--ir-in", required=True, type=Path)

    args = parser.parse_args(argv)
    _validate_common_inputs(parser, args)
    try:
        if args.command == "convert":
            return _do_convert(args)
        if args.command == "check-equivalence":
            return _do_check_equivalence(args)
        if args.command == "derive-bonded":
            return _do_derive_bonded(args)
        if args.command == "derive-torsions":
            return _do_derive_torsions(args)
        if args.command == "assign-charges":
            return _do_assign_charges(args)
        if args.command == "global-fit":
            return _do_global_fit(args)
        if args.command == "check-staleness":
            return _do_check_staleness(args)
    except Exception as exc:  # noqa: BLE001 - CLI 边界：任何错误都以 1 退出
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 2  # pragma: no cover - argparse 已处理


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
