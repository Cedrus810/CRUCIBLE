"""覆盖集诊断驱动（2026-11-20 计划 §2）：逐分子走入口配方并记录诊断，产出支持范围报告。

默认只做轻阶段（antechamber/parmchk2/tleap + 纯 Python + OpenMM Reference 单点），
不跑 MACE 扫描、不做 MD（扫描/液体属计划 §B/§C，另行报资源）；xtb 单点 opt-in。
逐分子阶段（--stages 逗号分隔）：

    entry      SDF 合同 + run_entry（antechamber -c bcc → parmchk2 → tleap，MVP-6 配方）；
               其余阶段都依赖它，必须排第一
    rotors     find_rotors 计数与等价类（纯 Python）
    frcmod     解析 parmchk2 frcmod 的 ATTN 缺项（纯 Python）
    export     validation.engines.compare 在入口 IR 上跑三后端等价（skip_missing_sander，
               gmx 缺失自动 SKIPPED；sander 存在时会做两腿单点，秒级）
    refcheck   参考势适用性的纯检查（MACE 支持元素/闭壳层：molecule_fallback_reason）
    mace_smoke MACE 单点冒烟（加载模型，分钟级；--stages 显式点名才跑）
    xtb_sp     xtb GFN2 单点 + 梯度（入口几何、不优化，``--grad --chrg Q --uhf 0``；
               秒级，--stages 显式点名才跑；线程数 --xtb-threads，默认 4）。
               开壳层/非法电荷在运行前拒绝（记 fail → 缺项，不阻塞）
    fallback   计划 §2.2 第 3 项：fast preset 的 bonded 拟合（xtb GFN-FF --ohess →
               mSeminario，线性角策略取 library 以便计数）；记兜底到 GAFF2 库值的项数
               （线性角）。fast preset 本身用 linear_angles=error，有线性角时 run 会拒绝——
               这条也记进缺项。分钟级以内，显式点名才跑（线程 --xtb-threads）
    scan_probe 计划 §2.2 第 6 项：取 find_rotors 的第一根可转键，用参考势（决策 0028：离子/两性离子
               xtb GFN2，其余 MACE）做一次真实
               relaxed scan（与 G2 同参数，10°）并计时；扫描缓存落在 <name>/scan_probe/scans，
               重跑同一分子直接命中缓存。每分子分钟到几十分钟，显式点名才跑
               （torch 线程 --mace-threads，设备 --device）

verdict 规则与字段定义在 ``ffcompiler.validation.coverage``（施工前锁定，不改结果再改规则）。

用法：
    python scripts/coverage_diagnose.py --table <molecules.json> \
        --outdir jobs/coverage_20261122 [--molecules ethanol,acetone] \
        [--stages entry,rotors,frcmod,export,refcheck] [--xtb-threads 4] \
        [--mace-threads 4] [--device cpu] [--force]

输出：``<outdir>/<name>/``（entry 配方产物 + record.json）、``summary.csv``、
``report.md``、``manifest.json``。``<name>/record.json`` 已存在时跳过该分子
（续跑按文件存在判断，不用 hash；同 scripts/torsion_relaxed_check.py 惯例）。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ffcompiler.validation.coverage import (  # noqa: E402
    MoleculeRecord,
    load_table,
    parse_frcmod_gaps,
    summarize,
    write_summary_csv,
)

LIGHT_STAGES = ("entry", "rotors", "frcmod", "export", "refcheck")
ALL_STAGES = LIGHT_STAGES + ("mace_smoke", "xtb_sp", "fallback", "scan_probe")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_record(record: MoleculeRecord, path: Path) -> None:
    """落 record.json（续跑按文件存在判断；diagnose_one 读回同一组字段）。先写临时文件再改名，
    中断时不会留下半截 record 被当成已完成。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record.to_dict(), ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    tmp.replace(path)


def _fail(error: str, **extra) -> dict:
    return {"status": "fail", "error": error, **extra}


def stage_entry(sdf: Path, workdir: Path, *, name: str, net_charge: int | None) -> dict:
    """SDF 合同 + 入口配方。失败返回 fail（MoleculeInputError/TypingError 消息自带日志路径）。"""
    from ffcompiler.molecule.io import MoleculeInputError
    from ffcompiler.typing.gaff2 import TypingError, run_entry

    try:
        entry = run_entry(sdf, workdir, name=name, net_charge=net_charge)
    except (MoleculeInputError, TypingError) as exc:
        return _fail(str(exc))
    ir = entry.ir
    # 下游阶段（diagnose_one）从这两个文件读回入口 IR；必须在这里落盘
    from ffcompiler.ir import coords as ir_coords
    from ffcompiler.ir import io as ir_io

    ir_io.to_json(ir, workdir / "ir.json")
    ir_coords.to_json(ir_coords.from_positions(entry.positions_nm, ir, source=f"entry:{sdf}"),
                      workdir / "coords.json")
    frcmod_artifact = entry.artifacts.get("frcmod")
    return {
        "status": "ok",
        "n_atoms": len(ir.atoms),
        "total_charge": int(ir.total_charge),
        "ambertools_version": entry.ambertools_version,
        "ir_json": str(workdir / "ir.json"),
        "coords_json": str(workdir / "coords.json"),
        "prmtop": str(entry.artifacts.get("prmtop", "")),
        "inpcrd": str(entry.artifacts.get("inpcrd", "")),
        "frcmod": str(frcmod_artifact) if frcmod_artifact else "",
        "sdf_sha256": ir.provenance["entry"]["sdf_sha256"],
    }


def stage_rotors(sdf: Path, ir, positions_nm: np.ndarray) -> dict:
    """可转键计数（rotors.py 规则）；entry 失败时不调用。"""
    from rdkit import Chem

    from ffcompiler.torsion.rotors import find_rotors

    mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
    if mol is None:
        return _fail(f"SDF 无法被 RDKit 解析：{sdf}")
    try:
        rotors = find_rotors(ir, mol, positions_nm)
    except Exception as exc:  # noqa: BLE001 - 诊断口径：拓扑不一致等按 fail 记录
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    return {
        "status": "ok",
        "n_rotors": len(rotors),
        "bonds": [list(r.bond) for r in rotors],
        "classes": [int(r.rotor_class) for r in rotors],
    }


def stage_frcmod(frcmod_path: str) -> dict:
    """frcmod ATTN 缺项；ATTN > 0 记为 gap（GAFF2 库缺参数，兜底来源）。"""
    if not frcmod_path or not Path(frcmod_path).is_file():
        return _fail(f"frcmod 不存在：{frcmod_path!r}")
    gaps = parse_frcmod_gaps(frcmod_path)
    out = {"status": "ok", **gaps}
    if gaps["n_attn"]:
        out["gaps"] = [f"frcmod_missing_params({gaps['n_attn']})"]
    return out


def stage_export(ir, positions_nm: np.ndarray, workdir: Path, *,
                 prmtop: str, inpcrd: str) -> dict:
    """入口 IR 的三后端等价（compare 口径不新增判据；sander/gmx 缺失走 SKIPPED）。"""
    from ffcompiler.validation.engines import compare

    kwargs = {}
    if prmtop and inpcrd and Path(prmtop).is_file() and Path(inpcrd).is_file():
        kwargs = {"tleap_prmtop": Path(prmtop), "tleap_inpcrd": Path(inpcrd)}
    try:
        report = compare(ir, positions_nm, workdir / "export", skip_missing_sander=True,
                         **kwargs)
    except Exception as exc:  # noqa: BLE001 - 诊断口径：引擎异常按 fail 记录
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    failed = [row for row in report.rows if not row["passed"]]
    failed_f = [row for row in report.force_rows if not row["passed"]]
    max_de = max((row["delta"] for row in report.rows), default=0.0)
    max_df = max((row["max_df"] for row in report.force_rows), default=0.0)
    gmx_precision = report.metadata.get("gmx", {}).get("precision")
    if (failed and not failed_f and not report.errors and gmx_precision == "mixed"
            and all(row["leg"] == "gmx" for row in failed)):
        # 用户决定 2026-11-24：混合精度 GROMACS 腿的超差只记精度缺项，不判 blocked
        # （同决策 0013"错在精度，那不应该成为门的问题"）；其余腿超差/出错仍是 blocker
        worst = max(failed, key=lambda r: r["delta"])
        rel = worst["delta"] / abs(worst["E_ref"]) if worst["E_ref"] else float("inf")
        cats = ",".join(sorted({r["category"] for r in failed}))
        return {"status": "ok", "n_rows": len(report.rows),
                "n_force_rows": len(report.force_rows), "max_dE": max_de, "max_dF": max_df,
                "skipped": dict(report.skipped), "gmx_precision": gmx_precision,
                "gmx_mixed_precision_failed_rows": failed,
                "gaps": [f"gmx_mixed_precision({cats}; max dE {worst['delta']:.3g} kJ/mol, "
                         f"rel {rel:.2g})"]}
    if report.errors or failed or failed_f:
        detail = json.dumps({"errors": report.errors,
                             "failed_energy": failed[:3], "failed_force": failed_f[:3]},
                            ensure_ascii=False)
        return _fail(detail[:2000], max_dE=max_de, max_dF=max_df)
    return {"status": "ok", "n_rows": len(report.rows),
            "n_force_rows": len(report.force_rows), "max_dE": max_de, "max_dF": max_df,
            "skipped": dict(report.skipped)}


def stage_refcheck(ir) -> dict:
    """参考势适用性纯检查（不加载 MACE/torch）。

    参考势按决策 0028 选择（与 derive-torsions --reference auto 同一判定）：离子/两性离子用
    xtb GFN2（元素范围 H–Rn），其余用 MACE。GFN2 能否在该分子上给出有限结果由 xtb_sp 阶段实测。
    """
    from ffcompiler.torsion import xtb_ref
    from ffcompiler.torsion.mlp import molecule_fallback_reason

    if xtb_ref.is_ion_or_zwitterion(ir):
        reference, supported = xtb_ref.MODEL_NAME, xtb_ref.SUPPORTED
    else:
        reference, supported = "mace", frozenset(range(1, 119))
    reason = molecule_fallback_reason(ir, supported)
    if reason == "odd_electrons":
        return _fail("odd_electrons（参考势闭壳层假设不成立）")
    gaps = []
    if reason == "unsupported_element":
        gaps.append(f"reference_{reference}_unsupported_element")
    return {"status": "ok", "reference": reference, "molecule_fallback_reason": reason,
            "gaps": gaps}


def stage_mace_smoke(ir, positions_nm: np.ndarray) -> dict:
    """MACE 单点冒烟：能量/力有限即通过。重 import（torch/mace），opt-in。"""
    from ase import Atoms

    from ffcompiler.torsion.mace import molecule_info, polar_calculator

    try:
        calc = polar_calculator()
        atoms = Atoms(numbers=[a.atomic_number for a in ir.atoms],
                      positions=np.asarray(positions_nm) * 10.0,  # nm → Å
                      info=molecule_info(ir))
        atoms.calc = calc
        e = float(atoms.get_potential_energy())
        f = atoms.get_forces()
    except Exception as exc:  # noqa: BLE001 - 诊断口径
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    if not np.all(np.isfinite(f)):
        return _fail(f"forces non-finite (E={e})")
    return {"status": "ok", "energy_eV": e}


def stage_xtb_sp(ir, positions_nm: np.ndarray, workdir: Path, *, threads: int = 4) -> dict:
    """xtb GFN2 单点 + 梯度（计划 §2.2 第 4 项）：能量/梯度有限即通过。

    只判"参考势能否在此分子上给出有限结果"，不设梯度阈值（入口几何不是 GFN2 极小，
    梯度大小只记录供看分布，不判定）。失败（含开壳层预检拒绝）记 fail → 缺项。
    """
    from ffcompiler.hessian.xtb import run_xtb_singlepoint
    from ffcompiler.units import BOHR_NM, HARTREE_KJ_MOL

    try:
        out = run_xtb_singlepoint(
            np.array([a.atomic_number for a in ir.atoms]), np.asarray(positions_nm),
            method="gfn2", charge=int(ir.total_charge), workdir=workdir, threads=threads)
    except Exception as exc:  # noqa: BLE001 - 诊断口径
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    grad = out["gradient_hartree_per_bohr"]
    max_grad = float(np.max(np.linalg.norm(grad, axis=1)))
    return {
        "status": "ok",
        "method": "gfn2",
        "energy_hartree": out["energy_hartree"],
        "max_atom_gradient_hartree_per_bohr": max_grad,
        "max_atom_gradient_kj_mol_nm": max_grad * HARTREE_KJ_MOL / BOHR_NM,
        "xtb_version": out["xtb_version"],
        "workdir": str(workdir),
    }


def stage_fallback(ir, positions_nm: np.ndarray, workdir: Path, *, threads: int = 4) -> dict:
    """计划 §2.2 第 3 项：fast preset 的 bonded 拟合后兜底计数（沿用现行 provenance 口径）。

    bonded 里唯一的兜底路径是线性角保留 GAFF2 库值（``linear_angle_policy="library"``，
    provenance ``bonded.linear_angles_kept``）；非正 k 按决策 0014 直接报错、不回退，
    记 fail。fast preset 的 run 固定 ``linear_angles=error``，因此有线性角的分子在单命令
    入口会被拒绝——同时记为缺项。Hessian 质量警告只记录，不判定。
    """
    from ffcompiler.bonded.assign import derive_bonded
    from ffcompiler.hessian.xtb import run_xtb_hessian

    try:
        h = run_xtb_hessian(
            np.array([a.atomic_number for a in ir.atoms]), np.asarray(positions_nm),
            method="gfnff", charge=int(ir.total_charge), workdir=workdir / "xtb",
            threads=threads)
        new_ir, _, report = derive_bonded(
            ir, h, method="mseminario", linear_angle_policy="library",
            input_positions_nm=np.asarray(positions_nm))
    except Exception as exc:  # noqa: BLE001 - 诊断口径
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    kept = [list(a) for a in report.linear_angles_kept]
    gaps = []
    if kept:
        gaps.append(f"bonded_linear_angles_library(n={len(kept)}; "
                    "fast preset linear_angles=error 会拒绝)")
    return {
        "status": "ok",
        "hessian": "gfnff",
        "method": "mseminario",
        "n_bonds": len(new_ir.bonds),
        "n_angles": len(new_ir.angles),
        "n_linear_angles_library": len(kept),
        "linear_angles_library": kept,
        "n_imaginary": int(report.diagnostics.n_imaginary),
        "hessian_warnings": list(report.diagnostics.warnings),
        "max_bond_ratio_dev": float(report.max_bond_ratio),
        "max_class_spread": max(report.class_spread.values(), default=0.0),
        "gaps": gaps,
    }


def stage_scan_probe(sdf: Path, ir, positions_nm: np.ndarray, workdir: Path, *,
                     threads: int = 4, device: str = "cpu", step_deg: float = 10.0) -> dict:
    """计划 §2.2 第 6 项：MACE 实扫一根可转键并计时（取 find_rotors 的第一根，确定性）。

    走现行 MLPTorsionProvider（fit_bonds 只选这一根，其余键 not_selected 不扫），所以
    记录的就是管线里的真实耗时与结果；该键兜底（scan_not_converged、mlp_error、整分子
    兜底等）记缺项。没有可转键 → skip（扫描就绪度 = 无需扫描）。
    """
    import time

    from rdkit import Chem

    from ffcompiler.torsion import xtb_ref
    from ffcompiler.torsion.mace import default_model_path, polar_calculator
    from ffcompiler.torsion.mlp import MLPTorsionProvider
    from ffcompiler.torsion.rotors import find_rotors

    try:
        mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
        if mol is None:
            return _fail(f"SDF 无法被 RDKit 解析：{sdf}")
        rotors = find_rotors(ir, mol, np.asarray(positions_nm))
        if not rotors:
            return {"status": "skip", "note": "no rotatable bonds"}
        bond = tuple(int(a) for a in rotors[0].bond)
        if xtb_ref.is_ion_or_zwitterion(ir):  # 决策 0028：离子/两性离子用 xtb GFN2
            reference = xtb_ref.MODEL_NAME
            binary = xtb_ref.xtb_binary()
            charge = int(ir.total_charge)
            provider = MLPTorsionProvider(
                lambda ir_: xtb_ref.XtbGFN2Calculator(charge, threads=threads, binary=binary),
                supported=xtb_ref.SUPPORTED, scan_kwargs={"step_deg": float(step_deg)},
                model_name=xtb_ref.MODEL_NAME, model_path=binary,
                scan_cache=workdir / "scans", fit_bonds=[bond])
        else:
            reference = "mace"
            model_path = default_model_path()
            provider = MLPTorsionProvider(
                lambda ir_: polar_calculator(model_path=model_path, device=device, threads=threads),
                scan_kwargs={"step_deg": float(step_deg)}, model_path=model_path,
                scan_cache=workdir / "scans", fit_bonds=[bond])
        t0 = time.perf_counter()
        _, report = provider.parameterize(ir, np.asarray(positions_nm), mol)
        wall = time.perf_counter() - t0
    except Exception as exc:  # noqa: BLE001 - 诊断口径
        return _fail(f"{type(exc).__module__}.{type(exc).__name__}: {exc}")
    row = next(r for r in report.rotors if tuple(r["bond"]) == bond)
    out = {
        "status": "ok",
        "bond": list(bond),
        "reference": reference,
        "n_rotors": len(rotors),
        "step_deg": float(step_deg),
        "threads": int(threads),
        "device": device,
        "wall_s": float(wall),
        "row_status": row["status"],
        "reason": row.get("reason"),
        "fit_rmse": row.get("fit_rmse"),
        "max_hysteresis": row.get("max_hysteresis"),
        "max_bond_change_pm": row.get("max_bond_change_pm"),
        "molecule_fallback": report.molecule_fallback,
        "gaps": [],
    }
    if row["status"] != "mlp":
        out["gaps"].append(f"scan_probe_fallback({row.get('reason')})")
    return out


def diagnose_one(row: dict, outdir: Path, stages: tuple[str, ...], *,
                 force: bool, xtb_threads: int = 4, mace_threads: int = 4,
                 device: str = "cpu") -> MoleculeRecord:
    """跑一个分子的诊断阶段并落 record.json；已存在且未 --force 则只读回。"""
    name = row["name"]
    workdir = outdir / name
    record_path = workdir / "record.json"
    if record_path.is_file() and not force:
        data = json.loads(record_path.read_text(encoding="utf-8"))
        return MoleculeRecord(
            name=data["name"], category=data["category"], net_charge=data["net_charge"],
            n_atoms=data["n_atoms"], sdf_sha256=data["sdf_sha256"],
            stages=data["stages"], verdict=data["verdict"], blockers=data["blockers"])

    if stages[0] != "entry":
        raise ValueError(f"{name}: stages 必须以 entry 开头（其余阶段依赖入口 IR）")
    sdf: Path = row["sdf"]
    stages_dict: dict[str, dict] = {"entry": stage_entry(
        sdf, workdir / "entry", name=name, net_charge=row.get("net_charge"))}
    entry_ir = None
    entry_pos = None
    if stages_dict["entry"]["status"] == "ok":
        entry_stage = stages_dict["entry"]
        from ffcompiler.ir import coords as ir_coords
        from ffcompiler.ir import io as ir_io

        entry_ir = ir_io.from_json(Path(entry_stage["ir_json"]))
        # 入口阶段落盘的 coords sidecar（nm，CLI 同一合同）；不要直接读 inpcrd——read_rst7 是 Å
        entry_pos = ir_coords.from_json(Path(entry_stage["coords_json"])).positions_nm
        for stage in stages[1:]:
            if stage == "rotors":
                stages_dict["rotors"] = stage_rotors(sdf, entry_ir, entry_pos)
            elif stage == "frcmod":
                stages_dict["frcmod"] = stage_frcmod(entry_stage.get("frcmod", ""))
            elif stage == "export":
                stages_dict["export"] = stage_export(
                    entry_ir, entry_pos, workdir,
                    prmtop=entry_stage.get("prmtop", ""),
                    inpcrd=entry_stage.get("inpcrd", ""))
            elif stage == "refcheck":
                stages_dict["refcheck"] = stage_refcheck(entry_ir)
            elif stage == "mace_smoke":
                stages_dict["mace_smoke"] = stage_mace_smoke(entry_ir, entry_pos)
            elif stage == "xtb_sp":
                stages_dict["xtb_sp"] = stage_xtb_sp(
                    entry_ir, entry_pos, workdir / "xtb_sp", threads=xtb_threads)
            elif stage == "fallback":
                stages_dict["fallback"] = stage_fallback(
                    entry_ir, entry_pos, workdir / "fallback", threads=xtb_threads)
            elif stage == "scan_probe":
                stages_dict["scan_probe"] = stage_scan_probe(
                    sdf, entry_ir, entry_pos, workdir / "scan_probe",
                    threads=mace_threads, device=device)
            else:
                raise ValueError(f"未知 stage：{stage}")
    record = MoleculeRecord(
        name=name, category=row["category"], net_charge=int(row.get("net_charge", 0)),
        n_atoms=len(entry_ir.atoms) if entry_ir is not None else 0,
        sdf_sha256=_sha256(sdf))
    record.stages = stages_dict
    record.finalize()
    _write_record(record, record_path)
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--table", required=True, type=Path,
                        help="覆盖集清单 JSON（name/category/sdf/net_charge）")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--molecules", default=None,
                        help="逗号分隔的分子名子集；缺省跑全表")
    parser.add_argument("--stages", default=",".join(LIGHT_STAGES),
                        help=f"逗号分隔，默认轻阶段 {','.join(LIGHT_STAGES)}；"
                             f"可选 {','.join(ALL_STAGES)}")
    parser.add_argument("--xtb-threads", type=int, default=4,
                        help="xtb_sp / fallback 阶段的 OMP/MKL 线程数（默认 4）")
    parser.add_argument("--mace-threads", type=int, default=4,
                        help="scan_probe 阶段的 torch 线程数（默认 4）")
    parser.add_argument("--device", default="cpu", help="scan_probe 的 MACE 设备（默认 cpu）")
    parser.add_argument("--force", action="store_true", help="已有 record.json 也重跑")
    args = parser.parse_args(argv)

    stages = tuple(s.strip() for s in args.stages.split(",") if s.strip())
    unknown = [s for s in stages if s not in ALL_STAGES]
    if unknown:
        parser.error(f"未知 stage：{unknown}；可选 {','.join(ALL_STAGES)}")
    if not stages or stages[0] != "entry":
        parser.error("--stages 必须以 entry 开头")

    table = load_table(args.table)
    if args.molecules:
        wanted = {s.strip() for s in args.molecules.split(",") if s.strip()}
        missing = wanted - {row["name"] for row in table}
        if missing:
            parser.error(f"--molecules 不在清单里：{sorted(missing)}")
        table = [row for row in table if row["name"] in wanted]

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    records: list[MoleculeRecord] = []
    for row in table:
        print(f"[coverage] {row['name']} ({row['category']}) ...", flush=True)
        try:
            record = diagnose_one(row, outdir, stages, force=args.force,
                                  xtb_threads=args.xtb_threads,
                                  mace_threads=args.mace_threads, device=args.device)
        except Exception as exc:  # noqa: BLE001 - 诊断驱动不得因单分子中断整批
            print(f"[coverage] {row['name']}: UNEXPECTED {exc}", flush=True)
            traceback.print_exc()
            record = MoleculeRecord(name=row["name"], category=row["category"],
                                    net_charge=int(row.get("net_charge", 0)),
                                    n_atoms=0, sdf_sha256="")
            record.stages = {"entry": _fail(f"unexpected: {exc}")}
            record.finalize()
        records.append(record)
        print(f"[coverage] {row['name']}: {record.verdict}"
              + (f" ({'; '.join(record.blockers)})" if record.blockers else ""), flush=True)

    write_summary_csv(records, outdir / "summary.csv")
    (outdir / "report.md").write_text(summarize(records), encoding="utf-8")
    manifest = {
        "task": "coverage diagnostic (plan 2026-11-20 §2)",
        "argv": sys.argv,
        "table": str(args.table),
        "table_sha256": _sha256(args.table),
        "stages": list(stages),
        "n_molecules": len(records),
        "verdicts": {v: sum(1 for r in records if r.verdict == v)
                     for v in ("ok", "ok_with_gaps", "blocked")},
    }
    (outdir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[coverage] report: {outdir / 'report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
