#!/usr/bin/env python
"""T4：生成 GAFF2 固件（12 个分子，MVP-0 计划）。

流程：RDKit 加氢 → ETKDGv3（种子 20261004）→ MMFF94 优化 → 写 SDF
→ antechamber (-c bcc -at gaff2) → parmchk2 → tleap (saveamberparm)。

任何一个分子任何一步失败：非 0 退出，并报告是哪个分子、哪一步。
固件提交进仓库；测试只读这些文件。重新生成要显式运行本脚本。

    python scripts/make_fixtures.py                 # 全部 12 个
    python scripts/make_fixtures.py --only ethanol  # 指定分子
    python scripts/make_fixtures.py --out <dir>     # 另一个输出目录
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
DEFAULT_OUT = REPO_ROOT / "tests" / "fixtures" / "gaff2"
SEED = 20261004

FIXTURES: list[tuple[str, str, int]] = [
    # (name, SMILES, 总电荷)
    ("ethanol", "CCO", 0),
    ("benzene", "c1ccccc1", 0),
    ("nma", "CC(=O)NC", 0),
    ("cyclopropane", "C1CC1", 0),
    ("cyclobutane", "C1CCC1", 0),
    ("norbornane", "C1CC2CCC1C2", 0),
    ("acetate", "CC(=O)[O-]", -1),
    ("methylammonium", "C[NH3+]", 1),
    ("chlorobenzene", "Clc1ccccc1", 0),
    ("methyl_benzoate", "COC(=O)c1ccccc1", 0),
    ("dmso", "CS(C)=O", 0),
    ("acetonitrile", "CC#N", 0),
]

REQUIRED_FILES = ("{name}.sdf", "{name}.mol2", "{name}.frcmod", "{name}.prmtop", "{name}.inpcrd")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ambertools_version() -> str:
    import os

    meta_dir = Path(os.environ.get("CONDA_PREFIX", "")) / "conda-meta"
    if meta_dir.is_dir():
        for meta in sorted(meta_dir.glob("ambertools-*.json")):
            data = json.loads(meta.read_text(encoding="utf-8"))
            return str(data.get("version", "unknown"))
    return "unknown"


def step_fail(name: str, step: str, detail: str) -> None:
    print(f"FAIL [{name}] {step}: {detail}", file=sys.stderr)
    raise SystemExit(2)


def _run(cmd: list[str], cwd: Path, name: str, step: str, log_name: str) -> None:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    (cwd / log_name).write_text(
        f"$ {' '.join(cmd)}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n"
        f"--- exit {proc.returncode} ---\n",
        encoding="utf-8",
    )
    if proc.returncode != 0:
        step_fail(name, step, f"{' '.join(cmd)} exited {proc.returncode}; see {log_name}")


def build_sdf(name: str, smiles: str, sdf: Path) -> int:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        step_fail(name, "rdkit", f"cannot parse SMILES {smiles!r}")
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = SEED
    if rdDistGeom.EmbedMolecule(mol, params) != 0:
        step_fail(name, "etkdg", "EmbedMolecule failed")
    if AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) != 0:
        print(f"WARN [{name}] MMFF94 optimization not fully converged", file=sys.stderr)
    mol.SetProp("_Name", name)
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    return mol.GetNumAtoms()


def _equalize_mol2_charges(mol2: Path, charge: int, name: str) -> None:
    """把 AM1-BCC 电荷的净残差均摊回每个原子，使 Σq 与 -nc 一致（计划要求 1e-4）。

    antechamber -c bcc 的电荷经 4 位小数舍入后总和会残留 ~1e-3 的漂移，
    且 -eq 等价化不强制归零总电荷。这里只改 @<TRIPOS>ATOM 的电荷列，
    不影响原子类型，parmchk2/tleap 结果的其余部分不变。
    """
    lines = mol2.read_text(encoding="utf-8").splitlines()
    atom_rows: list[int] = []
    in_atoms = False
    for n, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            continue
        if in_atoms and line.strip():
            atom_rows.append(n)
    if not atom_rows:
        step_fail(name, "charge-equalize", f"no ATOM records found in {mol2}")
    charges = [float(lines[n].split()[8]) for n in atom_rows]
    residual = charge - sum(charges)
    log = f"net charge target {charge}, raw sum {sum(charges):.6f}, residual {residual:.6f}\n"
    if abs(residual) > 1e-12:
        corr = residual / len(atom_rows)
        values = [q + corr for q in charges]
        # 以写出的精度为准：末原子吸收残余，使文件里的电荷和精确等于 -nc
        values[-1] += charge - sum(values)
        for n, value in zip(atom_rows, values):
            parts = lines[n].split()
            parts[8] = f"{value:.10f}"
            lines[n] = " ".join(parts)
        new_sum = sum(float(lines[n].split()[8]) for n in atom_rows)
        if abs(new_sum - charge) > 1e-8:
            step_fail(name, "charge-equalize", f"after adjustment sum is {new_sum}")
        log += f"adjusted by {corr:.6e} per atom; new sum {new_sum:.10f}\n"
    mol2.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (mol2.parent / "charge_equalize.log").write_text(log, encoding="utf-8")


def generate_one(name: str, smiles: str, charge: int, outdir: Path) -> dict:
    mdir = outdir / name
    mdir.mkdir(parents=True, exist_ok=True)
    sdf = mdir / f"{name}.sdf"
    mol2 = mdir / f"{name}.mol2"
    frcmod = mdir / f"{name}.frcmod"
    prmtop = mdir / f"{name}.prmtop"
    inpcrd = mdir / f"{name}.inpcrd"

    n_atoms = build_sdf(name, smiles, sdf)

    _run(
        ["antechamber", "-i", sdf.name, "-fi", "sdf", "-o", mol2.name, "-fo", "mol2",
         "-c", "bcc", "-at", "gaff2", "-nc", str(charge), "-pf", "y"],
        mdir, name, "antechamber", "antechamber.log",
    )
    if not mol2.is_file():
        step_fail(name, "antechamber", f"{mol2} not produced")
    _equalize_mol2_charges(mol2, charge, name)

    _run(
        ["parmchk2", "-i", mol2.name, "-f", "mol2", "-o", frcmod.name, "-s", "gaff2"],
        mdir, name, "parmchk2", "parmchk2.log",
    )
    if not frcmod.is_file():
        step_fail(name, "parmchk2", f"{frcmod} not produced")

    tleap_in = mdir / "tleap.in"
    tleap_in.write_text(
        "source leaprc.gaff2\n"
        f"MOL = loadmol2 {mol2.name}\n"
        f"loadamberparams {frcmod.name}\n"
        f"saveamberparm MOL {prmtop.name} {inpcrd.name}\n"
        "quit\n",
        encoding="utf-8",
    )
    _run(["tleap", "-f", "tleap.in"], mdir, name, "tleap", "tleap.log")
    if not (prmtop.is_file() and inpcrd.is_file()):
        step_fail(name, "tleap", "prmtop/inpcrd not produced")
    (mdir / "leap.log").unlink(missing_ok=True)

    from ffcompiler.formats.amber import read_prmtop

    prm = read_prmtop(prmtop)
    if prm.pointer("NATOM") != n_atoms:
        step_fail(name, "verify", f"prmtop has {prm.pointer('NATOM')} atoms, SDF has {n_atoms}")
    q_sum = sum(q / 18.2223 for q in prm["CHARGE"])
    if abs(q_sum - charge) > 1e-4:
        step_fail(name, "verify", f"charge sum {q_sum} != -nc {charge}")

    files = {tpl.format(name=name): sha256_file(mdir / tpl.format(name=name)) for tpl in REQUIRED_FILES}
    return {"smiles": smiles, "charge": charge, "n_atoms": n_atoms, "files": files}


def load_manifest(outdir: Path) -> dict:
    manifest_path = outdir / "manifest.json"
    if manifest_path.is_file():
        return json.loads(manifest_path.read_text(encoding="utf-8"))
    return {"molecules": {}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", action="append", default=[],
                        help="只生成指定分子（可重复）；默认全部")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)

    wanted = FIXTURES
    if args.only:
        known = {n for n, _, _ in FIXTURES}
        unknown = set(args.only) - known
        if unknown:
            raise SystemExit(f"unknown molecules: {sorted(unknown)}")
        wanted = [f for f in FIXTURES if f[0] in set(args.only)]

    manifest = load_manifest(outdir)
    manifest.setdefault("molecules", {})
    for name, smiles, charge in wanted:
        print(f"[make_fixtures] generating {name} ...", flush=True)
        manifest["molecules"][name] = generate_one(name, smiles, charge, outdir)
        print(f"[make_fixtures] {name} done", flush=True)

    manifest["generated"] = datetime.now().isoformat(timespec="seconds")
    manifest["seed"] = SEED
    manifest["ambertools_version"] = ambertools_version()
    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )
    print(f"[make_fixtures] manifest written to {outdir / 'manifest.json'}")


if __name__ == "__main__":
    main()
