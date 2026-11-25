"""MVP-7 多起点扫描与分支诊断（设计 §10.4；2026-11-20 计划 §B）的单元测试。

纯逻辑部分（聚类/合并/缓存身份）用合成 ScanResult；扫描层用 GAFF2 ethanol 的 MM
作参考势（同 test_torsion_scan.py 的既有做法，不用 MACE）。默认 {x0} 路径逐位不变
的回归检查在这里（缓存键无 start 块、主数组不受 collect_candidates 影响）。

provider 端到端用 24° 网格（15 点 ≥ Fourier 的 13 点下限）并 fit_bonds 限一根键，
控制 OpenMM Reference 上的松弛次数。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from ffcompiler.importers import from_amber
from ffcompiler.torsion.mlp import MLPTorsionProvider, scan_cache_key
from ffcompiler.torsion.mm import OpenMMCalculator
from ffcompiler.torsion.multistart import etkdg_starts, merge_start_scans
from ffcompiler.torsion.rotors import find_rotors
from ffcompiler.torsion.scan import ScanResult, relaxed_scan

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
ETHANOL = FIXTURES / "ethanol"


def _ethanol():
    ir, pos = from_amber(*(ETHANOL / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(ETHANOL / "ethanol.sdf"), removeHs=False)
    return ir, pos, mol


def _heavy_mask(ir) -> np.ndarray:
    return np.array([a.atomic_number > 1 for a in ir.atoms])


def _synthetic_scan(quartet, positions, energy, *, candidates=True):
    """合成 ScanResult：positions (n, N, 3) nm；候选块按 relaxed_scan 的形状约定。"""
    n, n_atoms = len(energy), positions.shape[1]
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    common = dict(
        quartet=quartet, phi_target_rad=phi, phi_rad=phi,
        energy=np.asarray(energy, dtype=np.float64),
        positions=np.asarray(positions, dtype=np.float64),
        hysteresis=np.zeros(n), converged=np.ones(n, dtype=bool),
        max_bond_change_pm=0.0, params={"step_deg": 60.0})
    if not candidates:
        return ScanResult(candidates=None, **common)
    cands = {
        "energy": np.stack([common["energy"], common["energy"]]),
        "converged": np.ones((2, n), dtype=bool),
        "positions": np.stack([common["positions"], common["positions"]]),
        "phi_target_rad": np.stack([phi, phi]),
        "phi_rad": np.stack([phi, phi]),
    }
    return ScanResult(candidates=cands, **common)


# ---- clustering --------------------------------------------------------------

def test_kabsch_align_recovers_rigid_motion():
    from ffcompiler.clustering import kabsch_align

    rng = np.random.default_rng(7)
    src = rng.normal(size=(12, 3))
    angle = 0.7
    rot = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                    [np.sin(angle), np.cos(angle), 0.0],
                    [0.0, 0.0, 1.0]])
    dst = src @ rot.T + np.array([1.0, -2.0, 0.5])
    assert np.max(np.abs(kabsch_align(src, dst) - dst)) < 1e-10


def test_kabsch_does_not_use_reflection():
    from ffcompiler.clustering import heavy_rmsd_nm

    # 非共面（手性）四点：镜像体只能靠反射对齐，纯旋转对不齐
    src = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    mirrored = src.copy()
    mirrored[3, 2] *= -1.0
    assert heavy_rmsd_nm(src, mirrored, np.ones(4, dtype=bool)) > 0.1
    # 平面分子不适用：平面图形的镜像在 3D 里可用 180° 旋转对齐（不用在这里测）


def test_greedy_cluster_two_well_separated_groups():
    from ffcompiler.clustering import greedy_cluster

    heavy = np.ones(8, dtype=bool)
    base = np.zeros((8, 3))
    base[:, 0] = np.arange(8) * 0.15  # 链式参考几何
    bent = base.copy()
    bent[4:, 1] += 0.3  # 内部位移（4/8 原子挪 0.3 nm），对齐后 RMSD ≈ 0.21 nm
    positions = np.stack([base, bent, base + 1e-3 * np.arange(8)[:, None], bent + 1e-3])
    energies = np.array([10.0, 0.0, 12.0, 2.0])  # bent 组整体更低
    ids = greedy_cluster(positions, energies, heavy, threshold_nm=0.05)
    # bent 组（低能量）成为簇 0，base 组成为簇 1
    assert ids.tolist() == [1, 0, 1, 0]


def test_greedy_cluster_threshold_merges_everything():
    from ffcompiler.clustering import greedy_cluster

    heavy = np.ones(3, dtype=bool)
    positions = np.stack([np.zeros((3, 3)), np.full((3, 3), 0.01), np.full((3, 3), 0.02)])
    ids = greedy_cluster(positions, np.array([5.0, 1.0, 3.0]), heavy, threshold_nm=0.05)
    assert np.all(ids == 0)


# ---- merge_start_scans -------------------------------------------------------

N = 6
BONDS = [(0, 1)]
HEAVY8 = np.ones(8, dtype=bool)


def _two_scans(bend_nm=0.0, extra_energy=0.0, second_unconverged_at=None):
    """x0 与一个内部构象不同的起点（原子 4–7 沿 y 挪 bend_nm；整体平移会被 Kabsch
    对齐消去，必须用内部位移才能形成不同分支）；extra_energy 加到第二起点的全部格点。"""
    rng = np.random.default_rng(3)
    x0_pos = rng.normal(size=(N, 8, 3)) * 0.1
    x0 = _synthetic_scan((0, 1, 2, 3), x0_pos, np.arange(N) * 1.0)
    b_pos = x0_pos.copy()
    b_pos[:, 4:, 1] += bend_nm
    b = _synthetic_scan((0, 1, 2, 3), b_pos, np.arange(N) * 1.0 + extra_energy)
    if second_unconverged_at is not None:
        b.candidates["converged"][0, second_unconverged_at] = False
        b.candidates["converged"][1, second_unconverged_at] = False
    return x0, b


def test_merge_identical_scans_reproduces_x0_exactly():
    x0, b = _two_scans()
    env, diag = merge_start_scans([x0, b], BONDS, HEAVY8, threshold_nm=0.05)
    assert np.array_equal(env.energy, x0.energy)
    assert np.array_equal(env.positions, x0.positions)
    assert diag["branch_count_max"] == 1
    assert diag["n_points_improved"] == 0
    assert diag["envelope_gain_max"] == pytest.approx(0.0, abs=1e-15)


def test_merge_envelope_picks_lower_branch_and_reports_it():
    x0, b = _two_scans(bend_nm=0.3, extra_energy=-5.0)  # 内部构象差 ≥ 阈值，能量低 5
    env, diag = merge_start_scans([x0, b], BONDS, HEAVY8, threshold_nm=0.05)
    assert diag["branch_count_max"] == 2
    assert diag["branch_gap_max"] == pytest.approx(5.0)
    assert diag["n_points_multibranch"] == N
    assert diag["envelope_gain_max"] == pytest.approx(5.0)
    assert diag["n_points_improved"] == N
    assert np.allclose(env.energy, np.arange(N) * 1.0 - 5.0)
    expected = x0.positions[0].copy()
    expected[4:, 1] += 0.3
    assert np.allclose(env.positions[0], expected)  # 包络取自第二起点的几何


def test_merge_prefers_converged_over_lower_unconverged():
    x0, b = _two_scans(bend_nm=0.3, extra_energy=-5.0, second_unconverged_at=2)
    env, diag = merge_start_scans([x0, b], BONDS, HEAVY8, threshold_nm=0.05)
    assert env.energy[2] == pytest.approx(x0.energy[2])  # 格点 2 退回 x0 的收敛候选
    assert env.converged[2]
    assert diag["n_points_multibranch"] == N - 1


def test_merge_old_cache_without_candidates_still_merges():
    x0 = _synthetic_scan((0, 1, 2, 3), np.zeros((N, 8, 3)), np.arange(N) * 1.0,
                         candidates=False)
    b_pos = np.zeros((N, 8, 3))
    b_pos[:, 4:, 1] = 0.3
    b = _synthetic_scan((0, 1, 2, 3), b_pos, np.arange(N) - 5.0, candidates=False)
    env, diag = merge_start_scans([x0, b], BONDS, HEAVY8, threshold_nm=0.05)
    assert np.allclose(env.energy, np.arange(N) - 5.0)
    assert diag["branch_count_max"] == 2


def test_merge_rejects_inconsistent_grids():
    x0 = _synthetic_scan((0, 1, 2, 3), np.zeros((N, 8, 3)), np.arange(N) * 1.0)
    short = _synthetic_scan((0, 1, 2, 3), np.zeros((N - 1, 8, 3)), np.arange(N - 1) * 1.0)
    with pytest.raises(ValueError, match="网格"):
        merge_start_scans([x0, short], BONDS, HEAVY8, threshold_nm=0.05)


# ---- 缓存身份与候选块 --------------------------------------------------------

def test_scan_cache_key_start_identity():
    base = dict(model="mm-self", quartet=(0, 1, 2, 3), info={"charge": 0, "spin": 1},
                scan_kwargs={"step_deg": 10.0}, model_file=None,
                numbers=[6, 8], bonds=[(0, 1)])
    legacy = scan_cache_key(**base)
    assert "start" not in json.loads(legacy)
    assert scan_cache_key(**base, start=None) == legacy  # x0 腿 = 现行键
    identity = {"kind": "etkdg", "label": "etkdg-0", "index": 0, "seed": 1}
    with_start = scan_cache_key(**base, start=identity)
    assert json.loads(with_start)["start"] == identity
    assert with_start != legacy


def test_scanresult_candidates_save_load_roundtrip(tmp_path):
    from ffcompiler.torsion.mlp import load_scan, save_scan

    n = 4
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    cands = {
        "energy": np.stack([np.arange(n) * 1.0, np.arange(n) * 2.0]),
        "converged": np.ones((2, n), dtype=bool),
        "positions": np.zeros((2, n, 5, 3)),
        "phi_target_rad": np.stack([phi, phi]),
        "phi_rad": np.stack([phi, phi]),
    }
    scan = ScanResult(quartet=(0, 1, 2, 3), phi_target_rad=phi, phi_rad=phi,
                      energy=np.arange(n) * 1.0, positions=np.zeros((n, 5, 3)),
                      hysteresis=np.zeros(n), converged=np.ones(n, dtype=bool),
                      max_bond_change_pm=0.0, params={"step_deg": 90.0},
                      candidates=cands)
    path = tmp_path / "1-2.npz"
    save_scan(path, "k1", scan)
    loaded = load_scan(path, "k1")
    assert loaded is not None and loaded.candidates is not None
    assert np.array_equal(loaded.candidates["energy"], cands["energy"])
    assert np.array_equal(loaded.candidates["phi_rad"], cands["phi_rad"])
    # 旧格式（无 cand_* 键）仍可读，candidates=None
    with np.load(path) as d:
        plain = {k: d[k] for k in d.files if not k.startswith("cand_")}
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, **plain)
    tmp.replace(path)
    old = load_scan(path, "k1")
    assert old is not None and old.candidates is None
    assert np.array_equal(old.energy, scan.energy)


# ---- relaxed_scan 扩展（MM 参考势，3 点小网格，不拟合） ------------------------

@pytest.fixture(scope="module")
def _ethanol_scan_pair():
    ir, pos, mol = _ethanol()
    numbers = [a.atomic_number for a in ir.atoms]
    bonds = [(b.i, b.j) for b in ir.bonds]
    rotor = {r.bond: r for r in find_rotors(ir, mol, pos)}[(1, 2)]  # C–O
    calc = OpenMMCalculator(ir)
    plain = relaxed_scan(calc, numbers, pos, rotor.scan_quartet, bonds, step_deg=120.0)
    with_c = relaxed_scan(calc, numbers, pos, rotor.scan_quartet, bonds,
                          step_deg=120.0, collect_candidates=True)
    origin = relaxed_scan(calc, numbers, pos, rotor.scan_quartet, bonds,
                          step_deg=120.0, grid_origin_rad=0.0)
    return ir, plain, with_c, origin


def test_collect_candidates_leaves_main_arrays_bit_identical(_ethanol_scan_pair):
    _, plain, with_c, _ = _ethanol_scan_pair
    assert np.array_equal(plain.energy, with_c.energy)
    assert np.array_equal(plain.positions, with_c.positions)
    assert np.array_equal(plain.phi_target_rad, with_c.phi_target_rad)
    assert np.array_equal(plain.hysteresis, with_c.hysteresis)
    assert plain.candidates is None and with_c.candidates is not None
    c = with_c.candidates
    assert c["energy"].shape == (2, 3) and c["positions"].shape == (2, 3, 9, 3)
    # 逐方向候选与主数组一致：主数组 = 收敛候选的最小值
    for s in range(3):
        pool = [c["energy"][d, s] for d in (0, 1) if c["converged"][d, s]] or \
               [c["energy"][d, s] for d in (0, 1)]
        assert with_c.energy[s] == pytest.approx(min(pool), abs=1e-12)


def test_grid_origin_aligns_grid_to_requested_angle(_ethanol_scan_pair):
    _, _, _, origin = _ethanol_scan_pair
    dev = np.angle(np.exp(1j * (origin.phi_target_rad - np.arange(3) * np.radians(120.0))))
    assert np.max(np.abs(dev)) < 1e-4  # 网格对齐 [0°, 120°, 240°]


# ---- ETKDG 起点（MM 参考势） --------------------------------------------------

def test_etkdg_starts_ethanol():
    from ffcompiler.clustering import heavy_rmsd_nm

    ir, pos, mol = _ethanol()
    heavy = _heavy_mask(ir)
    starts = etkdg_starts(ir, mol, pos, OpenMMCalculator(ir), n_confs=3, seed=11,
                          dedup_nm=0.05, opt_fmax=0.005, opt_max_steps=200)
    assert 0 < len(starts) <= 3
    assert len({s.label for s in starts}) == len(starts)
    for s in starts:
        assert s.identity["kind"] == "etkdg" and s.identity["label"] == s.label
        assert s.identity["seed"] == 11 and s.identity["n_confs"] == 3
        assert s.positions_nm.shape == (9, 3) and np.all(np.isfinite(s.positions_nm))
        assert heavy_rmsd_nm(s.positions_nm, pos, heavy) >= 0.05  # 对 x0 已去重
    for i, a in enumerate(starts):
        for b in starts[i + 1:]:
            assert heavy_rmsd_nm(a.positions_nm, b.positions_nm, heavy) >= 0.05


# ---- provider 端到端（MM 参考势，缓存真实落盘） --------------------------------

@pytest.fixture(scope="module")
def _ethanol_rotors():
    ir, pos, mol = _ethanol()
    rotors = find_rotors(ir, mol, pos)
    return rotors


def test_provider_multistart_end_to_end_and_refit(tmp_path):
    ir, pos, mol = _ethanol()
    cache = tmp_path / "scans"
    provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_),
                                  supported=frozenset(range(1, 119)),
                                  scan_kwargs={"step_deg": 24.0},  # 15 点 ≥ 13
                                  model_name="mm-self", scan_cache=cache,
                                  fit_bonds=[(1, 2)],  # 只扫 C–O，控制时长
                                  multistart={"n_confs": 2, "seed": 5})
    new_ir, report = provider.parameterize(ir, pos, mol)
    mlp_rows = [r for r in report.rotors if r["status"] == "mlp"]
    assert len(mlp_rows) == 1
    assert "multistart" in mlp_rows[0]
    assert mlp_rows[0]["multistart"]["n_points"] == 15
    assert mlp_rows[0]["multistart"]["start_labels"][0] == "x0"
    assert report.multistart is not None
    assert report.multistart["identities"][0] is None
    assert all(i["kind"] == "etkdg" for i in report.multistart["identities"][1:])
    prov = new_ir.provenance["torsion"]
    assert prov["multistart"]["spec"] == {"n_confs": 2, "seed": 5}
    # 缓存：x0 腿用现行文件名，ETKDG 腿带 label
    names = {p.name for p in cache.glob("*.npz")}
    assert any(".etkdg-" not in name for name in names)
    assert any(".etkdg-" in name for name in names)
    # 换电荷重拟合走缓存（禁止参考势），同一包络复现同一组拟合项
    from ffcompiler.charge.apply import assign_charges
    from ffcompiler.charge.loaded import LoadedCharges
    from ffcompiler.torsion.mlp import refit_torsions

    # IR 带 MLP 二面角：assign_charges 内部就会走 refit，必须给它缓存
    charged, _ = assign_charges(new_ir, LoadedCharges(), pos, mol, scan_cache=cache)
    refit_ir, refit_report = refit_torsions(charged, pos, mol, scan_cache=cache)
    refit_rows = [r for r in refit_report.rotors if r["status"] == "mlp"]
    assert len(refit_rows) == 1

    def by_bond(ir_):
        out: dict = {}
        for t in ir_.proper_torsions:
            out.setdefault(frozenset((t.j, t.k)), []).append(t)
        return out

    # 换电荷后 y 依据新电荷背景重算，force_constant 本来就允许不同；
    # 这里验证的是"同缓存重拟合确定"：显式 refit_torsions 与 assign_charges 内部
    # 那次 refit 的结果逐位一致（都只读缓存、不再调用参考势）。
    charged_map, refit_map = by_bond(charged), by_bond(refit_ir)
    for bond, terms in charged_map.items():
        for a, b in zip(terms, refit_map[bond]):
            assert (a.source == "mlp") == (b.source == "mlp")
            if a.source == "mlp":
                assert a.force_constant == b.force_constant
                assert a.phase == b.phase


def test_provider_default_path_has_no_multistart(tmp_path):
    ir, pos, mol = _ethanol()
    cache = tmp_path / "scans"
    provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_),
                                  supported=frozenset(range(1, 119)),
                                  scan_kwargs={"step_deg": 24.0},
                                  model_name="mm-self", scan_cache=cache,
                                  fit_bonds=[(1, 2)])
    new_ir, report = provider.parameterize(ir, pos, mol)
    assert report.multistart is None
    assert all("multistart" not in r for r in report.rotors)
    assert "multistart" not in new_ir.provenance["torsion"]
    # 缓存键与现行逐位一致：x0 文件的 key JSON 无 start 块、无候选块
    npz = next(cache.glob("*.npz"))
    with np.load(npz) as d:
        payload = json.loads(str(d["key"]))
        assert "start" not in payload
        assert "cand_energy" not in d.files
