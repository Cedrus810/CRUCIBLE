"""原生 AMBER prmtop / rst7 读写（决策 0015）的必测用例。

ParmEd 不再是运行依赖；环境里还有 ParmEd 时，用它作为解析结果的独立对照。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.formats.amber import (
    PrmtopError,
    read_prmtop,
    read_rst7,
    write_prmtop,
    write_rst7,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "gaff2"
NAMES = sorted(p.name for p in FIXTURES.iterdir() if p.is_dir())


@pytest.mark.parametrize("name", NAMES)
def test_read_matches_parmed_raw_data(name: str) -> None:
    """每个 FLAG 的原始数组与 ParmEd 的 parm_data 完全相同（CHARGE 为 ParmEd 的 e 单位 × 18.2223）。"""
    parmed = pytest.importorskip("parmed")
    path = FIXTURES / name / f"{name}.prmtop"
    ours = read_prmtop(path)
    ref = parmed.amber.AmberParm(str(path))
    assert list(ours.flags) == list(ref.flag_list)
    for flag in ref.flag_list:
        expected = list(ref.parm_data[flag])
        got = ours.flags[flag].values
        if flag == "CHARGE":
            np.testing.assert_allclose(np.array(got) / 18.2223, expected, rtol=0, atol=1e-12)
        else:
            assert got == expected, flag


@pytest.mark.parametrize("name", NAMES)
def test_write_read_roundtrip(name: str, tmp_path: Path) -> None:
    prm = read_prmtop(FIXTURES / name / f"{name}.prmtop")
    out = tmp_path / "copy.prmtop"
    write_prmtop(prm, out)
    again = read_prmtop(out)
    assert list(again.flags) == list(prm.flags)
    for flag, section in prm.flags.items():
        assert again.flags[flag].fmt == section.fmt
        assert again.flags[flag].values == section.values, flag


def test_float_format_is_fortran_e16_8(tmp_path: Path) -> None:
    prm = read_prmtop(FIXTURES / "ethanol" / "ethanol.prmtop")
    out = tmp_path / "e.prmtop"
    write_prmtop(prm, out)
    text = out.read_text().splitlines()
    i = text.index("%FLAG CHARGE".ljust(80))
    assert text[i + 1].startswith("%FORMAT(5E16.8)")
    assert len(text[i + 2]) == 80
    assert all(field[-4] == "E" for field in (text[i + 2][k:k + 16] for k in range(0, 80, 16)))


def test_pointer_accessor() -> None:
    prm = read_prmtop(FIXTURES / "ethanol" / "ethanol.prmtop")
    assert prm.pointer("NATOM") == 9
    assert prm.pointer("NTYPES") == 5


@pytest.mark.parametrize("name", NAMES)
def test_rst7_matches_parmed(name: str) -> None:
    parmed = pytest.importorskip("parmed")
    path = FIXTURES / name / f"{name}.inpcrd"
    coords = read_rst7(path)
    ref = parmed.amber.Rst7.open(str(path)).coordinates.reshape(-1, 3)
    np.testing.assert_array_equal(coords, ref)


def test_rst7_roundtrip(tmp_path: Path) -> None:
    coords = read_rst7(FIXTURES / "ethanol" / "ethanol.inpcrd")
    out = tmp_path / "x.rst7"
    write_rst7(coords, out, title="test")
    np.testing.assert_allclose(read_rst7(out), coords, rtol=0, atol=5e-8)


def test_truncated_flag_rejected(tmp_path: Path) -> None:
    text = (FIXTURES / "ethanol" / "ethanol.prmtop").read_text().splitlines()
    i = text.index(next(line for line in text if line.startswith("%FLAG MASS")))
    del text[i + 2]  # 删掉 MASS 的第一行数据
    bad = tmp_path / "bad.prmtop"
    bad.write_text("\n".join(text) + "\n")
    with pytest.raises(PrmtopError, match="MASS"):
        read_prmtop(bad)


def test_missing_required_flag_rejected(tmp_path: Path) -> None:
    text = (FIXTURES / "ethanol" / "ethanol.prmtop").read_text()
    bad = tmp_path / "bad.prmtop"
    bad.write_text(text.replace("%FLAG POINTERS", "%FLAG POINTERZ"))
    with pytest.raises(PrmtopError, match="POINTERS"):
        read_prmtop(bad)


def test_rst7_wrong_atom_count_rejected(tmp_path: Path) -> None:
    lines = (FIXTURES / "ethanol" / "ethanol.inpcrd").read_text().splitlines()
    bad = tmp_path / "bad.rst7"
    bad.write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(PrmtopError, match="coordinates"):
        read_rst7(bad)
