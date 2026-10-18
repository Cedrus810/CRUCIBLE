"""AM1-BCC 电荷：antechamber -c bcc（sqm），参数与 scripts/make_fixtures.py 相同（MVP-3 计划 §2.3）。"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.charge.base import ChargeError, ChargeProvider
from ffcompiler.charge.loaded import read_mol2_atoms
from ffcompiler.torsion.rotors import check_mol_topology


def antechamber_exe() -> str:
    """FFC_ANTECHAMBER 或 PATH 里的 antechamber。"""
    exe = os.environ.get("FFC_ANTECHAMBER") or shutil.which("antechamber")
    if not exe:
        raise ChargeError("antechamber not found (set FFC_ANTECHAMBER or put it on PATH)")
    return exe


def _ambertools_version(exe: str) -> str | None:
    meta = Path(exe).resolve().parent.parent / "conda-meta"
    found = sorted(meta.glob("ambertools-*.json")) if meta.is_dir() else []
    if not found:
        return None
    return json.loads(found[-1].read_text(encoding="utf-8")).get("version")


class AM1BCCCharges(ChargeProvider):
    method = "am1bcc"

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)

    def compute(self, ir, positions_nm, mol) -> np.ndarray:
        n = len(ir.atoms)
        check_mol_topology(ir, mol)  # 在调用 antechamber 之前（审阅 Important 1）
        self.workdir.mkdir(parents=True, exist_ok=True)
        mol = Chem.Mol(mol)
        conf = mol.GetConformer()
        for idx, xyz in enumerate(np.asarray(positions_nm, dtype=np.float64) * 10.0):
            conf.SetAtomPosition(idx, Point3D(*map(float, xyz)))
        Chem.MolToMolFile(mol, str(self.workdir / "input.sdf"))
        argv = [antechamber_exe(), "-i", "input.sdf", "-fi", "sdf", "-o", "am1bcc.mol2", "-fo", "mol2",
                "-c", "bcc", "-at", "gaff2", "-nc", str(ir.total_charge), "-pf", "y"]
        (self.workdir / "command.txt").write_text(" ".join(argv) + "\n", encoding="utf-8")
        out_mol2 = self.workdir / "am1bcc.mol2"
        out_mol2.unlink(missing_ok=True)
        proc = subprocess.run(argv, cwd=self.workdir, capture_output=True, text=True)
        log = self.workdir / "antechamber.log"
        log.write_text(f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}", encoding="utf-8")
        if proc.returncode != 0 or not out_mol2.is_file():
            tail = "\n".join(log.read_text(encoding="utf-8").splitlines()[-20:])
            raise ChargeError(f"antechamber failed (exit {proc.returncode}) in {self.workdir}:\n{tail}")
        rows = read_mol2_atoms(out_mol2)
        if len(rows) != n:
            raise TopologyMismatchError(f"antechamber mol2 has {len(rows)} atoms, IR has {n}")
        for idx, ((name, _), atom) in enumerate(zip(rows, ir.atoms)):
            element = re.sub(r"\d+$", "", name)
            if element.lower() != atom.element.lower():
                raise TopologyMismatchError(
                    f"antechamber atom {idx} {name!r} is not element {atom.element!r}")
        self._exe = argv[0]
        return np.array([q for _, q in rows], dtype=np.float64)

    def source_info(self) -> dict:
        exe = getattr(self, "_exe", None) or antechamber_exe()
        return {"antechamber": exe, "ambertools_version": _ambertools_version(exe),
                "workdir": str(self.workdir)}
