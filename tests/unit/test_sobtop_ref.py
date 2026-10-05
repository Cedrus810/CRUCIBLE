"""T16：Sobtop 参考数据读取与一致性的必测用例（非 external）。"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ffcompiler.bonded.sobtop_ref import SobtopParams, read_sobtop_params

REFS = Path(__file__).resolve().parents[1] / "fixtures" / "hessian" / "sobtop_refs"


def _counts_from_mol2(mol2_path: Path) -> tuple[int, int]:
    lines = mol2_path.read_text().splitlines()
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
            bonds.append((int(parts[1]), int(parts[2])))
    neighbors: dict[int, set[int]] = {}
    for a, b in bonds:
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    n_angles = sum(
        len(nbrs) * (len(nbrs) - 1) // 2 for nbrs in neighbors.values()
    )
    return len(bonds), n_angles


def test_read_sf6_drih_params() -> None:
    params = read_sobtop_params(REFS / "sf6_fchk" / "param_m4.txt")
    assert isinstance(params, SobtopParams)
    assert len(params.bonds) == 6
    assert all(k == pytest.approx(2928.942) for _, k in params.bonds.values())
    assert len(params.angles) == 15
    # 原子号 0 起算：S = 0（顶点），F = 1..6；90° 角与 180° 角都在
    assert (1, 0, 2) in params.angles
    angle_values = sorted(angle for angle, _ in params.angles.values())
    assert angle_values[0] == pytest.approx(90.0)
    assert angle_values[-1] == pytest.approx(180.0)


def test_case_counts_match_mol2_connectivity() -> None:
    manifest = json.loads((REFS / "manifest.json").read_text())
    cases = manifest["cases"]
    assert set(cases) == {
        "sf6_fchk", "methyl_benzoate_fchk",
        *(f"{n}_gfn2" for n in (
            "acetate", "benzene", "chlorobenzene", "cyclobutane", "cyclopropane",
            "dmso", "ethanol", "methylammonium", "methyl_benzoate", "nma",
            "norbornane",
        )),
    }
    for case in cases:
        n_bonds, n_angles = _counts_from_mol2(REFS / case / "input.mol2")
        for param_name in ("param_m2.txt", "param_m3.txt", "param_m4.txt"):
            params = read_sobtop_params(REFS / case / param_name)
            assert len(params.bonds) == n_bonds, (case, param_name)
            assert len(params.angles) == n_angles, (case, param_name)


def test_manifest_records_versions_and_hashes() -> None:
    manifest = json.loads((REFS / "manifest.json").read_text())
    assert manifest["sobtop_version"] == "2026.1.16"
    assert manifest["xtb_version"].startswith("6.")
    for case in manifest["cases"]:
        assert len(manifest["cases"][case]["mol2_sha256"]) == 64
        assert len(manifest["cases"][case]["hessian_sha256"]) == 64


def _make_refs_module():
    import importlib.util

    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location(
        "make_sobtop_refs", root / "scripts" / "make_sobtop_refs.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, root


def test_xtb_case_checks_reject_bad_hessian() -> None:
    """生成参考数据前，xtb 结果必须通过原子顺序、拓扑与 Hessian 质量检查。"""
    import numpy as np
    import pytest

    from ffcompiler.bonded.assign import TopologyMismatchError
    from ffcompiler.hessian.base import CartesianHessian
    from ffcompiler.hessian.xtb import read_xtb_hessian
    from ffcompiler.importers import from_amber

    refs, root = _make_refs_module()
    fx = root / "tests/fixtures/hessian/xtb_ethanol_gfn2"
    d = root / "tests/fixtures/gaff2/ethanol/ethanol"
    ir, pos = from_amber(Path(f"{d}.prmtop"), Path(f"{d}.inpcrd"), Path(f"{d}.sdf"))
    h = read_xtb_hessian(fx / "hessian", fx / "xtbopt.xyz", method="gfn2")
    refs.check_xtb_case(ir, pos, h)  # 合格的不报错

    perm = np.arange(len(h.atomic_numbers))[::-1]
    order = (3 * perm[:, None] + np.arange(3)).ravel()
    shuffled = CartesianHessian(
        atomic_numbers=h.atomic_numbers[perm], positions_nm=h.positions_nm[perm],
        hessian=h.hessian[np.ix_(order, order)], gradient=None, energy_kj=None,
        provider="xtb", method="gfn2",
    )
    with pytest.raises(TopologyMismatchError):
        refs.check_xtb_case(ir, pos, shuffled)
