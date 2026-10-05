"""T4：从头重新生成固件（external：调用 antechamber/parmchk2/tleap）。"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from rdkit import Chem

from ffcompiler.formats.amber import read_prmtop

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "make_fixtures.py"


@pytest.mark.external
def test_make_fixtures_ethanol_from_scratch(tmp_path: Path) -> None:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--only", "ethanol", "--out", str(tmp_path)],
        capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    mdir = tmp_path / "ethanol"
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    entry = manifest["molecules"]["ethanol"]
    assert entry["smiles"] == "CCO" and entry["charge"] == 0

    # prmtop 原子数与 SDF（含 H）一致
    sdf = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    prm = read_prmtop(mdir / "ethanol.prmtop")
    assert sdf is not None
    assert sdf.GetNumAtoms() == prm.pointer("NATOM")

    # 电荷总和与 -nc 一致到 1e-4
    assert abs(sum(q / 18.2223 for q in prm["CHARGE"]) - 0) <= 1e-4
