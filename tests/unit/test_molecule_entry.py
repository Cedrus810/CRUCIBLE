"""MVP-6 验收要点 2（设计 §3.1）：12 固件由原始 SDF 重建骨架 IR，与 prmtop 导入逐项等价。

只调 antechamber（-c bcc，2026-11-11 计划 §1 的开工裁定）/ parmchk2 / tleap，
不调 xtb/MACE/sander；单分子秒级，12 个合计约 1–3 min。
"""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
NAMES = sorted(p.name for p in FIXTURES.iterdir() if p.is_dir())


@pytest.mark.external
@pytest.mark.parametrize("name", NAMES)
def test_entry_rebuild_matches_fixture_prmtop(name, tmp_path):
    from ffcompiler.importers import from_amber
    from ffcompiler.typing.gaff2 import run_entry

    mdir = FIXTURES / name
    sdf = mdir / f"{name}.sdf"
    entry = run_entry(sdf, tmp_path / "entry", name=name)
    ref, _ = from_amber(mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", sdf)
    ir = entry.ir

    # 拓扑：名字、元素、GAFF2 类型
    assert [(a.name, a.atomic_number, a.atom_type) for a in ir.atoms] == \
           [(a.name, a.atomic_number, a.atom_type) for a in ref.atoms]
    # LJ 库值逐原子一致
    for a, b in zip(ir.atoms, ref.atoms):
        assert a.sigma == pytest.approx(b.sigma, abs=1e-12)
        assert a.epsilon == pytest.approx(b.epsilon, abs=1e-12)
    # 连接与 1-2/1-3/1-4 表
    assert {(min(b.i, b.j), max(b.i, b.j)) for b in ir.bonds} == \
           {(min(b.i, b.j), max(b.i, b.j)) for b in ref.bonds}
    assert sorted(ir.exclusions) == sorted(ref.exclusions)
    assert sorted(ir.pairs14) == sorted(ref.pairs14)
    assert ir.coulomb_14_scale == ref.coulomb_14_scale
    assert ir.vdw_14_scale == ref.vdw_14_scale
    # 库值参数逐项一致（bond/angle/proper/improper；顺序无关）
    bkey = lambda b: (min(b.i, b.j), max(b.i, b.j), round(b.r0, 9), round(b.k, 9))
    akey = lambda a: (a.i, a.j, a.k_atom, round(a.theta0, 9), round(a.force_constant, 9))
    tkey = lambda t: (t.i, t.j, t.k, t.l, t.periodicity,
                      round(t.phase, 6), round(t.force_constant, 9))
    assert sorted(map(bkey, ir.bonds)) == sorted(map(bkey, ref.bonds))
    assert sorted(map(akey, ir.angles)) == sorted(map(akey, ref.angles))
    assert sorted(map(tkey, ir.proper_torsions)) == sorted(map(tkey, ref.proper_torsions))
    assert sorted(map(tkey, ir.impropers)) == sorted(map(tkey, ref.impropers))
    # 电荷：同一配方同一 AmberTools，期望一致；按 1e-6 宽容（不作逐位硬判据）
    for a, b in zip(ir.atoms, ref.atoms):
        assert a.partial_charge == pytest.approx(b.partial_charge, abs=1e-6)
    assert ir.total_charge == ref.total_charge
    # provenance：source 不变 + entry 块（验收要点 3 前半）
    assert all(b.source == "gaff2-library" for b in ir.bonds)
    assert all(a.source == "gaff2-library" for a in ir.angles)
    assert all(t.source == "gaff2-library" for t in ir.proper_torsions)
    assert ir.provenance["entry"]["ambertools_version"]
    assert len(ir.provenance["entry"]["sdf_sha256"]) == 64
