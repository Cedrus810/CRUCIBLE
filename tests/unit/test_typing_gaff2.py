"""typing/gaff2.py 的可离线单元：电荷均摊、provenance 增强、执行器错误路径。"""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"

MOL2_HEAD = "@<TRIPOS>MOLECULE\nmol\n 3 2 1 0 0 0\nSMALL\nNO_CHARGES\n\n\n"


def _write_mol2(path: Path, charges: list[str]) -> None:
    rows = "\n".join(
        f"{i + 1} C{i} 0.0 0.0 0.0 c1 1 LIG {q}" for i, q in enumerate(charges)
    )
    path.write_text(
        MOL2_HEAD + "@<TRIPOS>ATOM\n" + rows + "\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n",
        encoding="utf-8")


def test_equalize_mol2_charges_sums_exact(tmp_path):
    from ffcompiler.typing.gaff2 import equalize_mol2_charges

    mol2 = tmp_path / "m.mol2"
    _write_mol2(mol2, ["0.1000", "0.2000", "0.6999"])  # 和 0.9999，目标 1
    log = equalize_mol2_charges(mol2, 1)
    in_atoms = False
    rows = []
    for line in mol2.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            continue
        if in_atoms and line.strip():
            rows.append(line.split())
    total = sum(float(r[8]) for r in rows)
    assert abs(total - 1.0) < 1e-9
    assert "adjusted" in log


def test_equalize_mol2_charges_rejects_short_row(tmp_path):
    from ffcompiler.typing.gaff2 import TypingError, equalize_mol2_charges

    mol2 = tmp_path / "m.mol2"
    mol2.write_text(MOL2_HEAD + "@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 c1 1 LIG\n",
                    encoding="utf-8")
    with pytest.raises(TypingError, match="列数不足"):
        equalize_mol2_charges(mol2, 0)


def test_equalize_mol2_charges_rejects_bad_number(tmp_path):
    from ffcompiler.typing.gaff2 import TypingError, equalize_mol2_charges

    mol2 = tmp_path / "m.mol2"
    _write_mol2(mol2, ["0.1", "0.2", "abc"])
    with pytest.raises(TypingError, match="不是数字"):
        equalize_mol2_charges(mol2, 0)


def test_run_raises_typing_error_on_missing_binary(tmp_path):
    from ffcompiler.typing.gaff2 import TypingError, _run

    with pytest.raises(TypingError, match="无法执行"):
        _run([str(tmp_path / "definitely-not-here"), "--x"], tmp_path, "x.log")


def test_ambertools_version_returns_string():
    from ffcompiler.typing.gaff2 import ambertools_version

    assert isinstance(ambertools_version(), str)


def test_file_stem_sanitized():
    from ffcompiler.typing.gaff2 import _file_stem

    assert _file_stem("lig and/x") == "lig_and_x"
    assert _file_stem("ethanol") == "ethanol"


def test_enrich_entry_provenance_keeps_sources():
    """验收要点 3 的前半：entry 块加入后，各项 source 仍是 gaff2-library。"""
    from ffcompiler.importers import from_amber
    from ffcompiler.molecule.io import MoleculeInput
    from ffcompiler.typing.gaff2 import _enrich_entry_provenance

    mdir = FIXTURES / "ethanol"
    ir, _ = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mi = MoleculeInput(path=mdir / "ethanol.sdf", name="ethanol", net_charge=0,
                       n_atoms=9, sha256="0" * 64)
    _enrich_entry_provenance(ir, mi, [["antechamber", "-i", "ethanol.sdf"]])
    assert ir.provenance["entry"]["ambertools_version"]
    assert ir.provenance["entry"]["sdf_sha256"] == "0" * 64
    assert ir.provenance["entry"]["geometry_checks"]["no_geometry_generation"] is True
    assert all(b.source == "gaff2-library" for b in ir.bonds)
    assert all(a.source == "gaff2-library" for a in ir.angles)
    assert all(t.source == "gaff2-library" for t in ir.proper_torsions)
