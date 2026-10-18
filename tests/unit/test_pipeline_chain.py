"""C2/C5：完整局部路线串联与组合验收（2026-10-15 计划 §5；2026-10-17 写，未运行）。

链条：基础 AMBER IR → assign-charges(loaded) → derive-bonded(--ir-in/--coords-in,
合成谐振 fchk) → convert(--ir-in/--coords-in, 三后端)。外加 C4 的
global-fit / check-staleness 组合行为与 C1 合同的拒绝路径。

2026-10-17 验收后在 181 上首次运行并修正（全文件约 1 s）。

合成 fchk 的 Hessian 是 IR 自身 bond/angle 谐振项的精确 Hessian
H = Σ_t k_t b_tᵀb_t（§5.1），PSD、平移不变，不需要 xtb。它没有二面角刚度：
乙醇的两个扭转模（C–C、C–O）频率为 0（数值噪声 ~1e-5 cm⁻¹，可正可负）。
global_fit 默认扣 H_fixed（二面角 + 非键），所以拟合解不是 s=1，只是有限正值。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from ffcompiler.cli import main
from ffcompiler.ir import io as ir_io

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
GAFF2_ETHANOL = FIXTURES / "ethanol"


def _ethanol_args() -> list[str]:
    return [
        "--prmtop", str(GAFF2_ETHANOL / "ethanol.prmtop"),
        "--inpcrd", str(GAFF2_ETHANOL / "ethanol.inpcrd"),
        "--sdf", str(GAFF2_ETHANOL / "ethanol.sdf"),
    ]


def _write_synthetic_fchk(path: Path) -> None:
    """ethanol IR 谐振 Hessian 写成 fchk（Ha/Bohr²、下三角行主序）。"""
    from ffcompiler.bonded.internal import wilson_b
    from ffcompiler.importers import from_amber
    from ffcompiler.units import BOHR_NM, HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2

    ir, positions = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    coords = (
        [("bond", (b.i, b.j)) for b in ir.bonds]
        + [("angle", (a.i, a.j, a.k_atom)) for a in ir.angles]
    )
    ks = np.array(
        [b.k for b in ir.bonds] + [a.force_constant for a in ir.angles],
        dtype=np.float64,
    )
    b_matrix = wilson_b(positions, coords)
    h = (b_matrix.T * ks) @ b_matrix  # kJ/mol/nm²
    n = len(ir.atoms)
    n3 = 3 * n
    tri = h[np.tril_indices(n3)] / HARTREE_KJ_MOL_PER_BOHR2_TO_KJ_MOL_PER_NM2
    coords_bohr = (positions / BOHR_NM).ravel()

    def _block(name: str, kind: str, count: int, values: list[float]) -> list[str]:
        lines = [f"{name:<42}{kind:>3}   N=  {count}"]
        for start in range(0, len(values), 5):
            lines.append("  " + " ".join(f"{v:.16e}" for v in values[start:start + 5]))
        return lines

    lines = [
        "synthetic harmonic fchk",
        "Freq      RB3LYP                                                      test",
        f"Number of atoms                            I                {n}",
    ]
    lines += _block("Atomic numbers", "I", n, [int(a.atomic_number) for a in ir.atoms])
    lines += _block("Current cartesian coordinates", "R", 3 * n, coords_bohr.tolist())
    lines += _block("Cartesian Force Constants", "R", len(tri), tri.tolist())
    lines.append("Total Energy                               R     -1.550000000000000E+02")
    path.write_text("\n".join(lines) + "\n")


def _charges_file(tmp_path: Path, ir, shift: float = 0.0) -> Path:
    q = np.array([a.partial_charge for a in ir.atoms])
    if shift:
        q[0] += shift
        q[1] -= shift
    path = tmp_path / f"q_{abs(shift):g}.txt"
    path.write_text("\n".join(repr(float(x)) for x in q))
    return path


# ---- C1 合同 ----

def test_coords_in_requires_ir_in(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        main([
            "convert", "--coords-in", str(tmp_path / "c.json"),
            "--output", "openmm", "--outdir", str(tmp_path / "out"),
        ])
    assert exc.value.code == 2


def test_missing_amber_inputs_exit_2(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["convert", "--prmtop", "x", "--output", "openmm",
              "--outdir", str(tmp_path / "out")])
    assert exc.value.code == 2


def test_derive_bonded_rejects_mismatched_sdf(tmp_path: Path, capsys) -> None:
    from ffcompiler.importers import from_amber

    ir, _ = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    ir_io.to_json(ir, tmp_path / "ir.json")
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)
    rc = main([
        "derive-bonded", "--ir-in", str(tmp_path / "ir.json"),
        "--prmtop", str(GAFF2_ETHANOL / "ethanol.prmtop"),
        "--inpcrd", str(GAFF2_ETHANOL / "ethanol.inpcrd"),
        "--sdf", str(FIXTURES / "acetate" / "acetate.sdf"),
        "--hessian-provider", "fchk", "--hessian-file", str(fchk),
        "--method", "mseminario", "--output", "openmm",
        "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 1
    assert "SDF" in capsys.readouterr().err  # 交叉检查导入三件套时即拒绝（原子数不符）


# ---- C2 串联 ----

def test_local_route_chain_preserves_params_and_coords(tmp_path: Path) -> None:
    """基础 IR → 电荷 → 局部 Hessian bond/angle → 同一最终 IR 三后端。"""
    from ffcompiler.importers import from_amber

    base_ir, _ = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    qfile = _charges_file(tmp_path, base_ir)

    # 1) 电荷
    q0 = tmp_path / "q0"
    assert main(["assign-charges", *_ethanol_args(), "--method", "loaded",
                 "--charges-file", str(qfile), "--output", "openmm",
                 "--outdir", str(q0)]) == 0
    assert (q0 / "ir.json").is_file() and (q0 / "coords.json").is_file()
    assert json.loads((q0 / "coords.json").read_text())["source"].startswith("inpcrd:")

    # 2) 局部 Hessian bond/angle（坐标取 q0 的 sidecar）
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)
    b1 = tmp_path / "b1"
    assert main([
        "derive-bonded", "--ir-in", str(q0 / "ir.json"),
        "--coords-in", str(q0 / "coords.json"),
        "--hessian-provider", "fchk", "--hessian-file", str(fchk),
        "--method", "mseminario", "--output", "openmm", "--outdir", str(b1),
    ]) == 0
    bonded_ir = ir_io.from_json(b1 / "ir.json")
    assert all(b.source == "mseminario" for b in bonded_ir.bonds)
    assert [a.partial_charge for a in bonded_ir.atoms] == pytest.approx(
        [a.partial_charge for a in ir_io.from_json(q0 / "ir.json").atoms])
    assert json.loads((b1 / "coords.json").read_text())["source"] == "hessian:gaussian-fchk"
    report = json.loads((b1 / "bonded_report.json").read_text())
    assert report["method"] == "mseminario"
    freqs = np.array(report["diagnostics"]["frequencies_cm"])
    assert freqs.size == 3 * 9 - 6
    # 合成 Hessian 无二面角刚度：恰两个零频扭转模，其余为正
    assert int(np.sum(np.abs(freqs) < 1.0)) == 2
    assert np.all(freqs[np.abs(freqs) >= 1.0] > 0.0)

    # 3) 同一最终 IR + 坐标 → 三后端
    final = tmp_path / "final"
    assert main([
        "convert", "--ir-in", str(b1 / "ir.json"), "--coords-in", str(b1 / "coords.json"),
        "--output", "amber", "gromacs", "openmm",
        "--outdir", str(final), "--ir-json", str(final / "ir.json"),
    ]) == 0
    final_ir = ir_io.from_json(final / "ir.json")
    assert final_ir.bonds == bonded_ir.bonds          # 非默认参数完整保留
    assert final_ir.atoms == bonded_ir.atoms
    for engine, fname in (("openmm", "system.xml"), ("amber", "molecule.prmtop"),
                          ("gromacs", "molecule.itp")):
        assert (final / engine / fname).is_file(), f"缺少 {engine}/{fname}"


@pytest.mark.parametrize("refinement_rounds", [0, 3])
def test_local_route_chain_with_torsions(tmp_path: Path, refinement_rounds: int) -> None:
    """C2 完整局部路线（10-15 计划 §5 C2 的顺序）：基础 IR → 局部 Hessian bond/angle → 最终电荷
    → 二面角 → 同一最终 IR + 坐标 → 三后端。每步都用 --ir-in/--coords-in 接上一步。

    二面角参考势用 ``--reference mm``（输入 IR 自身的 MM，扫描时冻结，单元测试的快速替身，
    不加载 MACE）；扫描步长 20°。"""
    from ffcompiler.importers import from_amber
    from ffcompiler.ir import coords as ir_coords

    sdf = str(GAFF2_ETHANOL / "ethanol.sdf")
    base_ir, _ = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)

    # 1) 局部 Hessian bond/angle（从 Amber 三件套起步）
    b1 = tmp_path / "b1"
    assert main(["derive-bonded", *_ethanol_args(), "--hessian-provider", "fchk",
                 "--hessian-file", str(fchk), "--method", "mseminario",
                 "--output", "openmm", "--outdir", str(b1)]) == 0
    bonded_ir = ir_io.from_json(b1 / "ir.json")
    hess_coords = ir_coords.from_json(b1 / "coords.json")

    # 2) 最终电荷：非默认值（两个原子 ±0.05 e），接 b1 的 IR 与 Hessian 几何
    qfile = _charges_file(tmp_path, bonded_ir, shift=0.05)
    q1 = tmp_path / "q1"
    assert main(["assign-charges", "--ir-in", str(b1 / "ir.json"), "--coords-in", str(b1 / "coords.json"),
                 "--sdf", sdf, "--method", "loaded", "--charges-file", str(qfile),
                 "--output", "openmm", "--outdir", str(q1)]) == 0
    charged_ir = ir_io.from_json(q1 / "ir.json")
    q_final = [a.partial_charge for a in charged_ir.atoms]
    assert q_final != [a.partial_charge for a in base_ir.atoms]
    assert charged_ir.bonds == bonded_ir.bonds and charged_ir.angles == bonded_ir.angles

    # 3) 二面角：在最终电荷 + Hessian bond/angle 上扫描、扣 MM、拟合
    t1 = tmp_path / "t1"
    cache = tmp_path / "scans"
    assert main(["derive-torsions", "--ir-in", str(q1 / "ir.json"), "--coords-in", str(q1 / "coords.json"),
                 "--sdf", sdf, "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
                 "--refinement-rounds", str(refinement_rounds),
                 "--output", "openmm", "--outdir", str(t1)]) == 0
    torsion_ir = ir_io.from_json(t1 / "ir.json")
    rotors = torsion_ir.provenance["torsion"]["rotors"]
    assert [r["status"] for r in rotors] == ["mlp", "mlp"]
    assert torsion_ir.provenance["torsion"]["charge_method"] == "loaded"
    assert torsion_ir.provenance["torsion"]["reference"]["charges"] == q_final
    assert torsion_ir.provenance["torsion"]["fit_protocol"]["refinement_rounds"] == refinement_rounds

    # 4) 同一最终 IR + 坐标 → 三后端
    final = tmp_path / "final"
    assert main(["convert", "--ir-in", str(t1 / "ir.json"), "--coords-in", str(t1 / "coords.json"),
                 "--output", "amber", "gromacs", "openmm",
                 "--outdir", str(final), "--ir-json", str(final / "ir.json")]) == 0
    final_ir = ir_io.from_json(final / "ir.json")
    # 三步的参数都完整保留
    assert final_ir.bonds == bonded_ir.bonds and final_ir.angles == bonded_ir.angles
    assert all(b.source == "mseminario" for b in final_ir.bonds)
    assert [a.partial_charge for a in final_ir.atoms] == q_final
    assert final_ir.proper_torsions == torsion_ir.proper_torsions
    assert any(t.source == "mlp" for t in final_ir.proper_torsions)
    for engine, fname in (("openmm", "system.xml"), ("amber", "molecule.prmtop"),
                          ("gromacs", "molecule.itp")):
        assert (final / engine / fname).is_file(), f"缺少 {engine}/{fname}"
    # 坐标：Hessian 几何一路沿 sidecar 传递，没有回退到 inpcrd
    for step, prev in ((q1, b1), (t1, q1)):
        c = ir_coords.from_json(step / "coords.json")
        assert c.source == f"coords:{prev / 'coords.json'}"
        assert c.positions_nm.tolist() == hess_coords.positions_nm.tolist()
    # 按设计顺序走完：依赖记录与最终 IR 一致，不过期
    assert main(["check-staleness", "--ir-in", str(final / "ir.json")]) == 0

    # 5) 再换一次电荷：用缓存重拟合二面角，依赖记录随之刷新，仍不过期
    q2file = _charges_file(tmp_path, final_ir, shift=-0.03)
    q2 = tmp_path / "q2"
    assert main(["assign-charges", "--ir-in", str(final / "ir.json"), "--coords-in", str(t1 / "coords.json"),
                 "--sdf", sdf, "--method", "loaded", "--charges-file", str(q2file),
                 "--scan-cache", str(cache), "--output", "openmm", "--outdir", str(q2)]) == 0
    refit_ir = ir_io.from_json(q2 / "ir.json")
    assert refit_ir.proper_torsions != final_ir.proper_torsions
    assert [r["status"] for r in refit_ir.provenance["torsion"]["rotors"]] == ["mlp", "mlp"]
    assert main(["check-staleness", "--ir-in", str(q2 / "ir.json")]) == 0
    assert refit_ir.provenance["torsion"]["fit_protocol"]["refinement_rounds"] == refinement_rounds


def test_derive_bonded_rejects_permuted_coords(tmp_path: Path, capsys) -> None:
    from ffcompiler.importers import from_amber

    ir, positions = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    from ffcompiler.ir import coords as ir_coords

    ir_io.to_json(ir, tmp_path / "ir.json")
    payload = json.loads((tmp_path / "ir.json").read_text())
    # 0、1 号原子都是 C：元素序列不变，只有原子名能查出换序
    payload["atoms"][0], payload["atoms"][1] = payload["atoms"][1], payload["atoms"][0]
    (tmp_path / "swapped.json").write_text(json.dumps(payload))
    # 坐标文件按未换序 IR 生成 → 与换序 IR 不配对
    ir_coords.to_json(ir_coords.from_positions(positions, ir, source="test"),
                      tmp_path / "coords.json")
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)
    rc = main([
        "derive-bonded", "--ir-in", str(tmp_path / "swapped.json"),
        "--coords-in", str(tmp_path / "coords.json"),
        "--hessian-provider", "fchk", "--hessian-file", str(fchk),
        "--method", "mseminario", "--output", "openmm",
        "--outdir", str(tmp_path / "out"),
    ])
    assert rc == 1
    assert "原子顺序不一致" in capsys.readouterr().err


# ---- C4 组合 ----

def test_global_fit_cli_ok_then_staleness_after_charge_change(tmp_path: Path) -> None:
    """global-fit 成功导出；换电荷后 check-staleness 报过期（拒绝静默沿用）。"""
    from ffcompiler.importers import from_amber

    base_ir, _ = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)
    gf = tmp_path / "gf"
    assert main([
        "global-fit", *_ethanol_args(),
        "--hessian-provider", "fchk", "--hessian-file", str(fchk),
        "--prior", "mseminario", "--output", "openmm", "--outdir", str(gf),
    ]) == 0
    report = json.loads((gf / "global_fit_report.json").read_text())
    assert report["status"] == "ok" and report["n_at_bound"] == 0
    fitted = ir_io.from_json(gf / "ir.json")
    assert fitted.provenance["bonded"]["method"] == "global_fit"
    assert (gf / "openmm" / "system.xml").is_file()

    # 刚拟合完：不过期
    assert main(["check-staleness", "--ir-in", str(gf / "ir.json")]) == 0

    # 换电荷（无 MLP 二面角，不需扫描缓存）→ H_fixed 依赖过期
    qfile = _charges_file(tmp_path, fitted, shift=0.05)
    q1 = tmp_path / "q1"
    assert main(["assign-charges", *_ethanol_args(), "--ir-in", str(gf / "ir.json"),
                 "--method", "loaded", "--charges-file", str(qfile),
                 "--output", "openmm", "--outdir", str(q1)]) == 0
    assert main(["check-staleness", "--ir-in", str(q1 / "ir.json")]) == 1


def test_global_fit_cli_lower_bound_writes_status_and_skips_engines(
    tmp_path: Path, monkeypatch, capsys,
) -> None:
    """下界命中（决策 0014）：report 记 failed_lower_bound、IR 落盘、不导出
    三后端、退出码 1。数学边界复用既有 test_lower_bound_hit，这里只测 CLI
    合同（monkeypatch global_fit 本体）。"""
    import ffcompiler.bonded.global_fit as gf_mod
    from ffcompiler.importers import from_amber

    base_ir, positions = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )

    class _Report:
        alpha = lam = cond = rel_residual_fit = rel_residual_prior = 0.0
        n_at_bound = 2
        scale_by_class = {}
        max_force_at_x0 = 0.0
        subtract_fixed = True

    monkeypatch.setattr(
        gf_mod, "global_fit",
        lambda *a, **k: (base_ir, positions, _Report()))
    out = tmp_path / "out"
    fchk = tmp_path / "synth.fchk"
    _write_synthetic_fchk(fchk)  # 读取器要有真文件；拟合本体被替换
    rc = main([
        "global-fit", *_ethanol_args(),
        "--hessian-provider", "fchk", "--hessian-file", str(fchk),
        "--output", "openmm", "--outdir", str(out),
    ])
    assert rc == 1
    report = json.loads((out / "global_fit_report.json").read_text())
    assert report["status"] == "failed_lower_bound" and report["n_at_bound"] == 2
    assert (out / "ir_failed_lower_bound.json").is_file()      # 诊断产物落盘，不用标准名
    assert not (out / "ir.json").exists()
    assert not (out / "openmm").exists()    # 不导出三后端
    assert "下界" in capsys.readouterr().err


def test_check_equivalence_ir_in_skips_tleap_leg(tmp_path: Path, monkeypatch, capsys) -> None:
    """--ir-in 的 IR 可能已改参数，原 prmtop 不代表它：不挂 sander-tleap 腿（2026-10-17 验收）。"""
    import ffcompiler.validation.engines as engines
    from ffcompiler.importers import from_amber

    seen = {}

    class _Report:
        passed = True

        def to_markdown(self) -> str:
            return "stub"

    def fake_compare(ir, positions, workdir, **kwargs):
        seen.update(kwargs)
        return _Report()

    monkeypatch.setattr(engines, "compare", fake_compare)
    ir, _ = from_amber(
        GAFF2_ETHANOL / "ethanol.prmtop", GAFF2_ETHANOL / "ethanol.inpcrd",
        GAFF2_ETHANOL / "ethanol.sdf",
    )
    ir_io.to_json(ir, tmp_path / "ir.json")
    assert main(["check-equivalence", "--ir-in", str(tmp_path / "ir.json"), *_ethanol_args(),
                 "--workdir", str(tmp_path / "w")]) == 0
    assert "tleap_prmtop" not in seen
    assert "sander-tleap" in capsys.readouterr().err
    # Amber 导入模式照常挂 tleap 腿
    assert main(["check-equivalence", *_ethanol_args(), "--workdir", str(tmp_path / "w2")]) == 0
    assert seen["tleap_prmtop"] == GAFF2_ETHANOL / "ethanol.prmtop"
