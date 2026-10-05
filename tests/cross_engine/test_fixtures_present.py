"""T4：固件存在性与完整性（非 external，只读仓库内的文件）。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"

EXPECTED = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
    "acetonitrile",
)
FILES = ("{n}.sdf", "{n}.mol2", "{n}.frcmod", "{n}.prmtop", "{n}.inpcrd")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.mark.parametrize("name", EXPECTED)
def test_fixture_files_present_with_matching_sha256(name: str) -> None:
    mdir = FIXTURE_DIR / name
    assert mdir.is_dir(), f"missing fixture dir {mdir}"
    for tpl in FILES:
        path = mdir / tpl.format(n=name)
        assert path.is_file(), f"missing {path}"


def test_fixture_manifest_consistent() -> None:
    manifest = json.loads((FIXTURE_DIR / "manifest.json").read_text(encoding="utf-8"))
    assert set(manifest["molecules"]) == set(EXPECTED)
    for name, entry in manifest["molecules"].items():
        assert entry["smiles"] and isinstance(entry["charge"], int)
        for filename, digest in entry["files"].items():
            path = FIXTURE_DIR / name / filename
            assert path.is_file(), f"manifest lists missing file {path}"
            assert sha256_file(path) == digest, f"sha256 mismatch for {path}"
