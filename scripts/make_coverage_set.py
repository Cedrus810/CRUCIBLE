#!/usr/bin/env python
"""生成覆盖集清单（2026-11-20 计划 §2.1；用户 2026-11-24 认可 34 个分子）。

34 = 12 个现有固件（直接引用 ``tests/fixtures/gaff2/<name>/<name>.sdf``，不复制、不重生成）
   + 22 个新增分子（本脚本生成 SDF）。

新增分子的 SDF 配方与固件完全相同，直接复用 ``scripts/make_fixtures.py::build_sdf``：
RDKit 加氢 → ETKDGv3（种子 20261004）→ MMFF94 优化 → 写 SDF。构象生成是本脚本这个显式
步骤（设计 §3.1/§26：入口本身不生成构象），种子固定，可重放。

输出（默认 ``benchmarks/coverage_set/``）：
    sdf/<name>.sdf        新增 22 个分子
    molecules.json        覆盖集清单（name/category/sdf/net_charge/smiles/source），
                          sdf 路径相对清单文件，供 ``scripts/coverage_diagnose.py --table``
已存在的 SDF 不覆盖（要重生成就先删掉；与 make_fixtures 一样是显式动作）。

    python scripts/make_coverage_set.py [--out benchmarks/coverage_set]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "src"))

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "gaff2"

# (name, category, SMILES, 净电荷)；类别与计划 §2.1 表逐行对应
FIXTURES: list[tuple[str, str]] = [
    ("ethanol", "alcohol_ether"),
    ("nma", "carbonyl_ester_amide"),
    ("methyl_benzoate", "carbonyl_ester_amide"),
    ("benzene", "aromatic_conjugated_fused"),
    ("chlorobenzene", "aromatic_conjugated_fused"),
    ("norbornane", "aromatic_conjugated_fused"),
    ("cyclopropane", "small_ring"),
    ("cyclobutane", "small_ring"),
    ("dmso", "heteroatom_halogen"),
    ("acetonitrile", "heteroatom_halogen"),
    ("acetate", "ion_zwitterion"),
    ("methylammonium", "ion_zwitterion"),
]

NEW: list[tuple[str, str, str, int]] = [
    ("methanol", "alcohol_ether", "CO", 0),
    ("diethyl_ether", "alcohol_ether", "CCOCC", 0),
    ("thf", "alcohol_ether", "C1CCOC1", 0),
    ("acetone", "carbonyl_ester_amide", "CC(C)=O", 0),
    ("methyl_acetate", "carbonyl_ester_amide", "COC(C)=O", 0),
    ("methylamine", "amine", "CN", 0),
    ("aniline", "amine", "Nc1ccccc1", 0),
    ("furan", "heteroaromatic", "c1ccoc1", 0),
    ("imidazole", "heteroaromatic", "c1c[nH]cn1", 0),
    ("pyridine", "heteroaromatic", "c1ccncc1", 0),
    ("styrene", "aromatic_conjugated_fused", "C=Cc1ccccc1", 0),
    ("naphthalene", "aromatic_conjugated_fused", "c1ccc2ccccc2c1", 0),
    ("anisole", "aromatic_conjugated_fused", "COc1ccccc1", 0),
    ("oxirane", "small_ring", "C1CO1", 0),
    ("methanethiol", "heteroatom_halogen", "CS", 0),
    ("dimethyl_sulfide", "heteroatom_halogen", "CSC", 0),
    ("thiophene", "heteroatom_halogen", "c1ccsc1", 0),
    ("dcm", "heteroatom_halogen", "ClCCl", 0),
    ("chloroform", "heteroatom_halogen", "ClC(Cl)Cl", 0),
    ("tfe", "heteroatom_halogen", "OCC(F)(F)F", 0),
    ("trimethyl_phosphate", "heteroatom_halogen", "COP(=O)(OC)OC", 0),
    ("glycine_zwitterion", "ion_zwitterion", "[NH3+]CC(=O)[O-]", 0),
]


def _formal_charge(sdf: Path) -> int:
    from rdkit import Chem

    mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
    if mol is None:
        raise SystemExit(f"RDKit 读不回 {sdf}")
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def build(out: Path) -> list[dict]:
    from make_fixtures import SEED, build_sdf

    manifest = json.loads((FIXTURE_DIR / "manifest.json").read_text(encoding="utf-8"))
    rows: list[dict] = []
    for name, category in FIXTURES:
        sdf = FIXTURE_DIR / name / f"{name}.sdf"
        if not sdf.is_file():
            raise SystemExit(f"固件缺 SDF：{sdf}")
        info = manifest["molecules"][name]
        rows.append({"name": name, "category": category,
                     "sdf": os.path.relpath(sdf, out), "net_charge": int(info["charge"]),
                     "smiles": info["smiles"], "source": "fixture (scripts/make_fixtures.py)"})
    sdf_dir = out / "sdf"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    for name, category, smiles, charge in NEW:
        sdf = sdf_dir / f"{name}.sdf"
        if not sdf.is_file():
            build_sdf(name, smiles, sdf)
        q = _formal_charge(sdf)
        if q != charge:
            raise SystemExit(f"{name}: SDF 形式电荷和 {q} ≠ 清单净电荷 {charge}")
        rows.append({"name": name, "category": category,
                     "sdf": os.path.relpath(sdf, out), "net_charge": charge, "smiles": smiles,
                     "source": f"RDKit ETKDGv3 seed {SEED} + MMFF94 (make_fixtures.build_sdf)"})
    names = [r["name"] for r in rows]
    if len(names) != len(set(names)) or len(rows) != 34:
        raise SystemExit(f"清单应为 34 个不重复分子，得到 {len(rows)}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "benchmarks" / "coverage_set")
    args = parser.parse_args()
    out = args.out.resolve()
    rows = build(out)
    (out / "molecules.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[coverage-set] {len(rows)} 个分子 → {out / 'molecules.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
