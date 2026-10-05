"""引擎适配器（MVP-0 计划 T9）：各引擎求能量/力、sander 单点、比较与报告。"""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import openmm
import openmm.app
from openmm import unit

from ffcompiler.formats.amber import read_prmtop, write_rst7
from ffcompiler.exporters.openmm import FORCE_GROUPS, to_openmm_system
from ffcompiler.ir import ForceFieldIR
from ffcompiler.ir.energy import COULOMB_KJ_NM, _bond_angle, evaluate

KCAL_KJ = 4.184


def reference_context(system: openmm.System) -> openmm.Context:
    """Reference 平台 Context（计划 §1.3：全部 float64、只用 Reference）。

    不能裸建 Context：默认平台在本机会落到 CUDA（混合精度），能量只在
    ~1e-7 相对误差内一致。Reference 平台本身是双精度、无 Precision 属性。
    """
    platform = openmm.Platform.getPlatformByName("Reference")
    return openmm.Context(system, openmm.VerletIntegrator(0.001), platform)


# §2.4 容差（决策 0004.1：内存内容差相对项放宽到 1e-8，覆盖双精度求和地板）
LINEAR_ANGLE_DEG = 175.0
_IN_MEM_E = lambda e: 1e-6 + 1e-8 * abs(e)          # 内存内，kJ/mol
_IN_MEM_F = 1e-6                                     # 内存内，kJ/mol/nm
_VIA_FILE_E = lambda e: 1e-4 + 1e-7 * abs(e)         # 经文件，kJ/mol
_VIA_FILE_F = 1e-3                                   # 经文件，kJ/mol/nm
_SANDER_E = lambda e: 2e-4 + 1e-7 * abs(e)           # sander 非静电，kcal/mol
_GMX_MIXED_REL = 1e-5                                # mixed 版 gmx 相对项（决策 0013）
_SANDER_EEL = lambda e: 2e-4 + 5e-5 * abs(e)         # sander 静电（库仑常数差异），kcal/mol


def openmm_energy_forces(
    system: openmm.System,
    positions_nm: np.ndarray,
    *,
    groups: dict[str, int] | None = None,
) -> tuple[dict[str, float], np.ndarray]:
    """按力组（给定映射）或按力的类别汇总能量，并返回力 (N,3) kJ/mol/nm。

    groups=None 时按力的类别：HarmonicBondForce→bond、HarmonicAngleForce→angle、
    PeriodicTorsionForce→torsion_all、NonbondedForce/CustomBondForce/
    CustomNonbondedForce→nonbonded_all。
    """
    platform = openmm.Platform.getPlatformByName("Reference")
    if groups is None:
        # 按力的类别聚合：把每个力的组号设为索引（必须在建 Context 之前）
        for idx in range(system.getNumForces()):
            system.getForce(idx).setForceGroup(idx)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    context.setPositions([openmm.Vec3(*p) for p in positions_nm] * unit.nanometer)

    energies: dict[str, float] = {}
    if groups is not None:
        for name, group in groups.items():
            state = context.getState(getEnergy=True, groups=1 << group)
            energies[name] = state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )
        energies["total"] = sum(energies.values())
    else:
        categories = {
            "HarmonicBondForce": "bond",
            "HarmonicAngleForce": "angle",
            "PeriodicTorsionForce": "torsion_all",
            "NonbondedForce": "nonbonded_all",
            "CustomBondForce": "nonbonded_all",
            "CustomNonbondedForce": "nonbonded_all",
        }
        energies = {v: 0.0 for v in set(categories.values())}
        for idx in range(system.getNumForces()):
            cls = system.getForce(idx).__class__.__name__
            key = categories.get(cls)
            if key is None:
                continue
            state = context.getState(getEnergy=True, groups=1 << idx)
            energies[key] += state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )
        energies["total"] = sum(energies.values())

    forces = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    del context
    return energies, np.asarray(forces, dtype=np.float64)


_SANDER_KEYS = ("bond", "angle", "dihed", "vdw", "eel", "vdw14", "eel14")

_MDIN = """single point energy at input coordinates
 &cntrl
  imin=1, maxcyc=1, ntb=0, igb=0, cut=999.0, ntpr=1,
 /
"""


_SANDER_ENERGY_PATTERNS = {
    "bond": r"BOND\s*=\s*([-\d.EeD+]+)",
    "angle": r"ANGLE\s*=\s*([-\d.EeD+]+)",
    "dihed": r"DIHED\s*=\s*([-\d.EeD+]+)",
    "vdw": r"\bVDWAALS\s*=\s*([-\d.EeD+]+)",
    "eel": r"\bEEL\s*=\s*([-\d.EeD+]+)",
    "vdw14": r"1-4\s+VDW\w*\s*=\s*([-\d.EeD+]+)",
    "eel14": r"1-4\s+EEL\s*=\s*([-\d.EeD+]+)",
}


def parse_mdout_energies(text: str) -> dict[str, float]:
    """解析 sander mdout 的第一个能量块（kcal/mol）；缺任何一个键就报错。

    Amber 24 的 imin=1 输出：能量表以 ``NSTEP ENERGY RMS ...`` 列头行开始、
    分量行以 ``RESTRAINT`` 结尾；中间可能夹杂空行，所以不能用空行截断。
    """
    header = re.search(r"^\s*NSTEP\b", text, re.M)
    if header is None:
        raise ValueError("no energy block found in sander output")
    rest = text[header.start():]
    end = rest.find("RESTRAINT")
    if end >= 0:
        newline = rest.find("\n", end)
        block = rest[: newline + 1 if newline >= 0 else len(rest)]
    else:
        block = "\n".join(rest.splitlines()[:40])
    out: dict[str, float] = {}
    missing = []
    for key, pat in _SANDER_ENERGY_PATTERNS.items():
        mm = re.search(pat, block)
        if mm is None:
            missing.append(key)
            continue
        out[key] = float(mm.group(1).replace("D", "E"))
    if missing:
        raise ValueError(f"missing energies {missing} in sander output block:\n{block[:400]}")
    return out


def sander_energies(prmtop: Path, inpcrd: Path, workdir: Path) -> dict[str, float]:
    """sander 单点能量（kcal/mol），取输入坐标处的第一个能量块。"""
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    mdin = workdir / "sander.mdin"
    mdout = workdir / "sander.mdout"
    mdin.write_text(_MDIN, encoding="utf-8")
    proc = subprocess.run(
        ["sander", "-O", "-i", str(mdin), "-o", str(mdout),
         "-p", str(Path(prmtop).resolve()), "-c", str(Path(inpcrd).resolve())],
        capture_output=True, text=True, timeout=600, cwd=workdir,
    )
    if proc.returncode != 0:
        tail = "\n".join((mdout.read_text(errors="replace").splitlines()[-30:])
                         if mdout.exists() else proc.stderr.splitlines()[-30:])
        raise RuntimeError(f"sander failed (exit {proc.returncode}):\n{tail}")
    return parse_mdout_energies(mdout.read_text(errors="replace"))


def _numpy_kcal_categories(ir: ForceFieldIR, positions_nm: np.ndarray) -> dict[str, float]:
    """numpy-IR 能量换算到 sander 的类别与单位（kcal/mol）。"""
    e = evaluate(ir, positions_nm)
    return {
        "bond": e["bond"] / KCAL_KJ,
        "angle": e["angle"] / KCAL_KJ,
        "dihed": (e["proper"] + e["improper"]) / KCAL_KJ,
        "vdw": e["lj"] / KCAL_KJ,
        "eel": e["coulomb"] / KCAL_KJ,
        "vdw14": e["lj14"] / KCAL_KJ,
        "eel14": e["coulomb14"] / KCAL_KJ,
        "total": e["total"] / KCAL_KJ,
    }


@dataclass
class CrossEngineReport:
    rows: list[dict] = field(default_factory=list)
    force_rows: list[dict] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return bool(self.rows) and not self.errors and all(row["passed"] for row in self.rows) and \
            all(row["passed"] for row in self.force_rows)

    def to_markdown(self) -> str:
        lines = [
            "| leg | reference | category | E_leg | E_ref | delta | tol | passed |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for row in self.rows:
            lines.append(
                f"| {row['leg']} | {row['reference']} | {row['category']} "
                f"| {row['E_leg']:.10g} | {row['E_ref']:.10g} "
                f"| {row['delta']:.3g} | {row['tol']:.3g} | {row['passed']} |"
            )
        lines.append("")
        lines.append("| leg | reference | max|dF| | tol | passed |")
        lines.append("|---|---|---|---|---|")
        for row in self.force_rows:
            lines.append(
                f"| {row['leg']} | {row['reference']} | {row['max_df']:.3g} "
                f"| {row['tol']:.3g} | {row['passed']} |"
            )
        for leg, reason in self.skipped.items():
            lines.append(f"SKIPPED {leg}: {reason}")
        for leg, reason in self.errors.items():
            lines.append(f"ERROR {leg}: {reason}")
        if "gmx" in self.metadata:
            lines.append(f"GROMACS precision={self.metadata['gmx']['precision']}; strict file tolerance")
        return "\n".join(lines)

    def to_json(self) -> dict:
        return {
            "rows": self.rows,
            "force_rows": self.force_rows,
            "skipped": self.skipped,
            "passed": self.passed,
            "errors": self.errors,
            "metadata": self.metadata,
        }


def gmx_energy_tolerance(
    category: str, e_ref: float, total_mag: float, precision: str | None
) -> float:
    """gmx 腿能量容差（kJ/mol，决策 0013）。

    mixed 版：坐标、参数、运算都是 fp32，实测偏差 ~1e-6 相对（≈4–23 ε₃₂），
    容差取 1e-4 + 1e-5|E|；total 行的相对项按最大分量量级（同 0004.1）。
    真实转换错误 ≥1e-3 相对，仍有约 100 倍判别余量。
    double 版或精度未知：原经文件容差 1e-4 + 1e-7|E|。
    """
    if precision != "mixed":
        return _VIA_FILE_E(e_ref)
    scale = total_mag if category == "total" else abs(e_ref)
    return 1e-4 + _GMX_MIXED_REL * scale


def _add_row(report: CrossEngineReport, leg: str, reference: str, category: str,
             e_leg: float, e_ref: float, tol: float) -> None:
    delta = abs(e_leg - e_ref)
    report.rows.append({
        "leg": leg, "reference": reference, "category": category,
        "E_leg": e_leg, "E_ref": e_ref, "delta": delta, "tol": tol,
        "passed": bool(delta <= tol),
    })


def compare(
    ir: ForceFieldIR,
    positions_nm: np.ndarray,
    workdir: Path,
    *,
    tleap_prmtop: Path | None = None,
    tleap_inpcrd: Path | None = None,
    amber_ir: ForceFieldIR | None = None,
    gromacs_ir: ForceFieldIR | None = None,
    include_gromacs: bool = True,
) -> CrossEngineReport:
    """跑 §5 T9 的各条腿并套 §2.4 容差。

    amber_ir / gromacs_ir 默认与 ir 相同（测试注入篡改导出，只进对应格式的文件）。
    """
    from ffcompiler.exporters.amber import write_amber
    from ffcompiler.exporters.gromacs import write_gromacs

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    report = CrossEngineReport()

    ref = evaluate(ir, positions_nm)
    for key in ("bond", "angle", "proper", "improper", "lj", "coulomb", "lj14",
                "coulomb14", "total"):
        _add_row(report, "numpy-ir", "reference-itself", key, ref[key], ref[key],
                 _IN_MEM_E(ref[key]))

    torsion_all = ref["proper"] + ref["improper"]
    nonbonded_all = ref["lj"] + ref["coulomb"] + ref["lj14"] + ref["coulomb14"]
    ref5 = {
        "bond": ref["bond"], "angle": ref["angle"],
        "torsion_all": torsion_all, "nonbonded_all": nonbonded_all,
        "total": ref["total"],
    }
    # total 行受非键/成键抵消影响，噪声地板由最大分量决定（决策 0004.1）
    total_mag = max(abs(ref[k]) for k in
                    ("bond", "angle", "proper", "improper", "lj", "coulomb",
                     "lj14", "coulomb14"))

    # ---- openmm-ir（内存内，8 项 + total）------------------------------
    e_groups, f_openmm_ir = openmm_energy_forces(
        to_openmm_system(ir, decompose=True), positions_nm, groups=FORCE_GROUPS,
    )
    for key in FORCE_GROUPS:
        _add_row(report, "openmm-ir", "numpy-ir", key, e_groups[key], ref[key],
                 _IN_MEM_E(ref[key]))
    _add_row(report, "openmm-ir", "numpy-ir", "total", e_groups["total"], ref["total"],
             _IN_MEM_E(total_mag))

    # ---- openmm-xml（生产模式写盘再反序列化，内存内容差）----------------
    system_prod = to_openmm_system(ir, decompose=False)
    xml_path = workdir / "system.xml"
    xml_path.write_text(openmm.XmlSerializer.serialize(system_prod), encoding="utf-8")
    system_xml = openmm.XmlSerializer.deserialize(xml_path.read_text())
    e_xml, f_xml = openmm_energy_forces(system_xml, positions_nm)
    for key in ("bond", "angle", "torsion_all", "nonbonded_all", "total"):
        e_ref = ref5["total"] if key == "total" else ref5[key]
        tol = _IN_MEM_E(total_mag) if key == "total" else _IN_MEM_E(e_ref)
        _add_row(report, "openmm-xml", "numpy-ir", key, e_xml[key], e_ref, tol)

    # ---- prmtop / itp 导出（经文件容差）--------------------------------
    amber_target = amber_ir if amber_ir is not None else ir
    amber_files = write_amber(amber_target, positions_nm, workdir / "amber")
    gromacs_target = gromacs_ir if gromacs_ir is not None else ir
    gromacs_files = write_gromacs(gromacs_target, positions_nm, workdir / "gromacs")

    amber = openmm.app.AmberPrmtopFile(str(amber_files["prmtop"]))
    system_prm = amber.createSystem(
        nonbondedMethod=openmm.app.NoCutoff, constraints=None, rigidWater=False,
        removeCMMotion=False,
    )
    e_prm, f_prm = openmm_energy_forces(system_prm, positions_nm)
    for key in ("bond", "angle", "torsion_all", "nonbonded_all", "total"):
        e_ref = ref5["total"] if key == "total" else ref5[key]
        tol = _VIA_FILE_E(total_mag) if key == "total" else _VIA_FILE_E(e_ref)
        _add_row(report, "openmm-prmtop", "numpy-ir", key, e_prm[key], e_ref, tol)

    gmx_top = openmm.app.GromacsTopFile(str(gromacs_files["top"]))
    system_top = gmx_top.createSystem(
        nonbondedMethod=openmm.app.NoCutoff, constraints=None, rigidWater=False,
        removeCMMotion=False,
    )
    e_top, f_top = openmm_energy_forces(system_top, positions_nm)
    for key in ("bond", "angle", "torsion_all", "nonbonded_all", "total"):
        e_ref = ref5["total"] if key == "total" else ref5[key]
        tol = _VIA_FILE_E(total_mag) if key == "total" else _VIA_FILE_E(e_ref)
        _add_row(report, "openmm-top", "numpy-ir", key, e_top[key], e_ref, tol)

    # ---- sander 两条腿 --------------------------------------------------
    numpy_kcal = _numpy_kcal_categories(ir, positions_nm)
    sander_ours = sander_energies(
        amber_files["prmtop"], amber_files["inpcrd"], workdir / "sander_ours",
    )
    # 近线性角（θ > 175°）：sander 与 OpenMM/numpy 的角度求值路径存在引擎间
    # 系统差（决策 0004.2），sander-vs-numpy 的 angle/dihed 行跳过；同引擎的
    # sander-tleap 对比保留。
    near_linear: list[tuple[int, int, int]] = []
    for angle in ir.angles:
        theta_deg = math.degrees(_bond_angle(
            positions_nm[angle.i], positions_nm[angle.j], positions_nm[angle.k_atom],
        ))
        if theta_deg > LINEAR_ANGLE_DEG:
            near_linear.append((angle.i, angle.j, angle.k_atom))
    if near_linear:
        report.skipped["sander-ours:angle/dihed"] = (
            f"near-linear angle(s) {near_linear[:3]} (> {LINEAR_ANGLE_DEG} deg); "
            "engine-specific sander-vs-numpy difference (decision 0004.2)"
        )
    for key in _SANDER_KEYS:
        if near_linear and key in ("angle", "dihed"):
            continue
        _add_row(report, "sander-ours", "numpy-ir", key, sander_ours[key],
                 numpy_kcal[key], _SANDER_EEL(numpy_kcal[key])
                 if key in ("eel", "eel14") else _SANDER_E(numpy_kcal[key]))

    if tleap_prmtop is not None and tleap_inpcrd is not None:
        # sander-tleap 腿必须在同一几何上比 Hamiltonian：给 tleap prmtop 配
        # 当前 positions 的 rst7（tleap 自带的 inpcrd 是原始几何，扰动构型下
        # 几何不一致会虚报失败）。
        tleap_dir = workdir / "sander_tleap"
        tleap_dir.mkdir(parents=True, exist_ok=True)
        natom_tleap = read_prmtop(tleap_prmtop).pointer("NATOM")
        if natom_tleap != len(positions_nm):
            raise ValueError(f"tleap prmtop 有 {natom_tleap} 个原子，坐标有 {len(positions_nm)} 个")
        tleap_crd = tleap_dir / "tleap.rst7"
        write_rst7(np.asarray(positions_nm, dtype=np.float64) * 10.0, tleap_crd, title="tleap")
        sander_tleap = sander_energies(tleap_prmtop, tleap_crd, tleap_dir)
        for key in _SANDER_KEYS:
            _add_row(report, "sander-tleap", "sander-ours", key, sander_tleap[key],
                     sander_ours[key], _SANDER_E(sander_ours[key]))

    # ---- gmx 腿（没有就 SKIPPED，决策 D1）-------------------------------
    if include_gromacs:
        from .gromacs import find_gromacs, energies as gromacs_energies
        try:
            executable = find_gromacs()
            if executable is None:
                report.skipped["gmx"] = "gmx: not installed (roadmap decision D1)"
            else:
                e_gmx, metadata = gromacs_energies(ir, positions_nm, gromacs_files,
                                                 workdir / "gmx", executable=executable)
                report.metadata["gmx"] = metadata
                for category in ref:
                    _add_row(report, "gmx", "numpy-ir", category, e_gmx[category],
                             ref[category],
                             gmx_energy_tolerance(category, ref[category], total_mag,
                                                  metadata.get("precision")))
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            report.errors["gmx"] = str(exc)

    # ---- 力（OpenMM 系四条腿，参照 openmm-ir）---------------------------
    # 参照取 production 内存的系统：与 XML/prmtop/top 同为 NonbondedForce 内核。
    # decompose 系统用 CustomNonbondedForce 表达同一 Hamiltonian，但算术路径
    # 不同，大力（~700 kJ/mol/nm）下与 NonbondedForce 的差异 ~1e-9 相对
    # （实测 1.1e-6 绝对），会虚报失败。
    _, f_openmm_ir = openmm_energy_forces(system_prod, positions_nm)
    for leg, forces in (
        ("openmm-xml", f_xml), ("openmm-prmtop", f_prm), ("openmm-top", f_top),
    ):
        max_df = float(np.max(np.abs(forces - f_openmm_ir)))
        tol = _IN_MEM_F if leg == "openmm-xml" else _VIA_FILE_F
        report.force_rows.append({
            "leg": leg, "reference": "openmm-ir", "max_df": max_df, "tol": tol,
            "passed": bool(max_df <= tol),
        })

    return report
