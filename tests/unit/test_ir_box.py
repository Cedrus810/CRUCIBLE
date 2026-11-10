"""IR 周期盒（决策 0023，语义按决策 0025 修订）单元测试：schema/JSON 往返、
盒子不改变参考能量且与 OpenMM 导出一致、三导出器盒写出、AMBER 周期 prmtop 可被
sander 读入、validation。"""

from __future__ import annotations

from pathlib import Path

import dataclasses
import pytest

from ffcompiler.exporters.amber import build_prmtop, write_amber
from ffcompiler.exporters.gromacs import write_gromacs
from ffcompiler.exporters.openmm import to_openmm_system
from ffcompiler.formats.amber import read_prmtop, read_rst7
from ffcompiler.importers import from_amber
from ffcompiler.ir import io as ir_io
from ffcompiler.ir.energy import evaluate
from ffcompiler.ir.schema import ForceFieldIR
from ffcompiler.ir.validation import IRValidationError, validate
from ffcompiler.validation.engines import openmm_energy_forces, sander_energies

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
BOX = ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))


@pytest.fixture()
def ethanol():
    ir, pos = from_amber(
        FIXTURES / "ethanol/ethanol.prmtop",
        FIXTURES / "ethanol/ethanol.inpcrd",
        FIXTURES / "ethanol/ethanol.sdf",
    )
    return ir, pos


def _with_box(ir: ForceFieldIR, box) -> ForceFieldIR:
    return dataclasses.replace(ir, box_nm=box)


def test_schema_default_and_validation(ethanol):
    ir, _ = ethanol
    assert ir.box_nm is None
    assert validate(_with_box(ir, BOX)) is None
    with pytest.raises(IRValidationError, match="行列式"):
        validate(_with_box(ir, ((3.0, 0, 0), (0, 3.0, 0), (0, 0, -3.0))))
    with pytest.raises(IRValidationError, match="box_nm"):
        validate(_with_box(ir, ((3.0, 0.0), (0.0, 3.0))))


def test_json_roundtrip(ethanol, tmp_path):
    ir, _ = ethanol
    boxed = _with_box(ir, ((2.5, 0.1, 0.0), (0.0, 2.5, 0.0), (0.2, 0.0, 2.6)))
    p = tmp_path / "ir.json"
    ir_io.to_json(boxed, p)
    back = ir_io.from_json(p)
    assert back.box_nm == boxed.box_nm
    ir_io.to_json(ir, tmp_path / "ir0.json")
    assert ir_io.from_json(tmp_path / "ir0.json").box_nm is None


def test_box_is_metadata_for_reference_energy(ethanol):
    # box_nm 只写进输出文件，参考能量仍是非周期哈密顿量（决策 0025）
    ir, pos = ethanol
    assert evaluate(_with_box(ir, BOX), pos) == evaluate(ir, pos)


@pytest.mark.parametrize("shift", [0.0, 2.9])
def test_boxed_openmm_export_matches_reference(ethanol, shift):
    # 带盒导出的 OpenMM System 与 NumPy 参考是同一个量；整体平移到盒边附近
    # （原子坐标跨过 L=3 nm）也一致——两边都不取镜像
    ir, pos = ethanol
    boxed = _with_box(ir, BOX)
    x = pos + shift
    energies, _ = openmm_energy_forces(to_openmm_system(boxed), x)
    assert energies["total"] == pytest.approx(evaluate(boxed, x)["total"], abs=1e-6, rel=1e-8)


def test_gromacs_gro_box_line(ethanol, tmp_path):
    ir, pos = ethanol
    out = write_gromacs(_with_box(ir, BOX), pos, tmp_path)
    last = out["gro"].read_text().splitlines()[-1].split()
    assert [float(v) for v in last] == pytest.approx([3.0, 3.0, 3.0])
    # 三斜：9 值，含非对角
    tri = ((3.0, 0.0, 0.0), (0.5, 3.0, 0.0), (0.2, 0.1, 3.0))
    out2 = write_gromacs(_with_box(ir, tri), pos, tmp_path / "tri")
    vals = [float(v) for v in out2["gro"].read_text().splitlines()[-1].split()]
    assert len(vals) == 9
    assert vals == pytest.approx([3.0, 3.0, 3.0, 0.0, 0.0, 0.5, 0.0, 0.2, 0.1])
    # 无盒：保持 10 Å 立方旧行为
    out3 = write_gromacs(ir, pos, tmp_path / "vac")
    assert out3["gro"].read_text().splitlines()[-1].split()[:3] == ["10.00000"] * 3
    # 非下三角 → 报错（gro 盒行表达不了，与 OpenMM 同一检查）
    with pytest.raises(ValueError, match="下三角"):
        write_gromacs(_with_box(ir, ((3.0, 0.5, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))),
                      pos, tmp_path / "bad")


def test_amber_prmtop_ifbox_and_rst7(ethanol, tmp_path):
    ir, pos = ethanol
    out = write_amber(_with_box(ir, BOX), pos, tmp_path)
    prm = read_prmtop(out["prmtop"])
    assert prm.pointer("IFBOX") == 1
    box_dims = prm["BOX_DIMENSIONS"]
    assert box_dims[0] == pytest.approx(90.0)
    assert box_dims[1:] == pytest.approx([30.0, 30.0, 30.0])
    assert list(prm["SOLVENT_POINTERS"]) == [1, 1, 2]
    assert list(prm["ATOMS_PER_MOLECULE"]) == [len(ir.atoms)]
    lines = out["inpcrd"].read_text().splitlines()
    assert [float(v) for v in lines[-1].split()] == pytest.approx(
        [30.0, 30.0, 30.0, 90.0, 90.0, 90.0])
    # 无盒：IFBOX=0，无盒行
    out2 = write_amber(ir, pos, tmp_path / "vac")
    assert read_prmtop(out2["prmtop"]).pointer("IFBOX") == 0
    # 坐标读回不受盒行影响
    assert read_rst7(out["inpcrd"]).shape == (len(ir.atoms), 3)


@pytest.mark.external
def test_sander_reads_periodic_prmtop(ethanol, tmp_path):
    # IFBOX=1 的 prmtop 必须带 SOLVENT_POINTERS / ATOMS_PER_MOLECULE，否则 sander
    # 读入即报错；单点（ntb=0）能量与无盒导出相同
    ir, pos = ethanol
    boxed = write_amber(_with_box(ir, BOX), pos, tmp_path / "box")
    plain = write_amber(ir, pos, tmp_path / "vac")
    e_box = sander_energies(boxed["prmtop"], boxed["inpcrd"], tmp_path / "run_box")
    e_vac = sander_energies(plain["prmtop"], plain["inpcrd"], tmp_path / "run_vac")
    assert e_box == pytest.approx(e_vac, abs=1e-4)


def test_amber_rejects_triclinic(ethanol):
    ir, _ = ethanol
    tri = ((3.0, 0.0, 0.0), (0.5, 3.0, 0.0), (0.2, 0.1, 3.0))
    with pytest.raises(ValueError, match="正交"):
        build_prmtop(_with_box(ir, tri))


def test_openmm_periodic_vectors(ethanol):
    ir, pos = ethanol
    system = to_openmm_system(_with_box(ir, BOX))
    (a, b, c) = system.getDefaultPeriodicBoxVectors()
    assert (a.x, a.y, a.z) == pytest.approx((3.0, 0.0, 0.0))
    assert (b.x, b.y, b.z) == pytest.approx((0.0, 3.0, 0.0))
    assert (c.x, c.y, c.z) == pytest.approx((0.0, 0.0, 3.0))
    # 下三角三斜
    tri = ((3.0, 0.0, 0.0), (0.5, 3.0, 0.0), (0.2, 0.1, 3.0))
    system2 = to_openmm_system(_with_box(ir, tri))
    (a2, b2, c2) = system2.getDefaultPeriodicBoxVectors()
    assert (b2.x, b2.y) == pytest.approx((0.5, 3.0))
    # 非下三角 → 报错
    with pytest.raises(ValueError, match="下三角"):
        to_openmm_system(_with_box(ir, ((3.0, 0.5, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))))
