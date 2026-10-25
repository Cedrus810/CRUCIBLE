"""fflib（Sobtop 数据文件解析与匹配）单元测试。

全部期望值来自 Sobtop 2026.1.16 的动态实验输出（沙箱 itp 逐字段核对），
见 docs/theory/sobtop_reverse_engineering.md 附录 A。
"""

from __future__ import annotations

import textwrap

import pytest

from ffcompiler.fflib import (
    GUESS_ANGLE_K_KJ_MOL_RAD2,
    GUESS_BOND_K_KJ_MOL_NM2,
    GUESS_IMPROPER_BARRIER_KCAL,
    GUESS_IMPROPER_PN,
    find_angle,
    find_bond,
    find_dihedrals,
    find_improper,
    find_lj,
    improper_center_types,
    lj_to_sigma_eps,
    parse_bondcrit,
    parse_bonded_param,
    parse_lj_param,
)

BONDED_PARAM = textwrap.dedent(
    """\
    ; 测试夹具：条目数值取自出厂 bonded_param.dat（GAFF 1.81 / parm19 摘录）
    $BOND
    c2-c2  569.4    1.3343
    c2-ha  343.0    1.0879
    $end
    $ANGLE
    c2-c2-ha   49.90    120.430
    ha-c2-ha   38.10    116.900
    $end
    $DIHEDRAL
    X -c3-c3-X    9    1.400       0.000           3.000
    hc-c3-c3-oh   1    0.00        0.0            -3.
    hc-c3-c3-oh   1    0.25        0.0             1.
    X -c2-c2-X    4   26.600     180.000           2.000
    $end
    $HARM_DIH
    $end
    $IMPROPER
    X -o -c -o          1.1          180.          2.
    X -X -ca-ha         1.1          180.          2.
    X -X -c -o         10.5          180.          2.
    $end
    """
)

LJ_PARAM = textwrap.dedent(
    """\
    ;This file contains LJ parameters that can be assigned by Sobtop
    Ctest    1.9255      0.105   ;My carbon
    UF_C      1.9255      0.105   ;UFF atom type of C
    UF_O      1.7500      0.060   ;UFF atom type of O
    UF_H      1.4430      0.044   ;UFF atom type of H
    """
)

BONDCRIT = textwrap.dedent(
    """\
    ;Fe Cl 2.0 2.4
    Fe Cl 2.0 2.4
    Fe C  0.0 2.1
    """
)


@pytest.fixture()
def bonded_file(tmp_path):
    p = tmp_path / "bonded_param.dat"
    p.write_text(BONDED_PARAM)
    return p


def test_parse_lj_param(tmp_path):
    p = tmp_path / "LJ_param.dat"
    p.write_text(LJ_PARAM)
    entries = parse_lj_param(p)
    assert [e.name for e in entries] == ["Ctest", "UF_C", "UF_O", "UF_H"]
    assert entries[1].r_half_a == pytest.approx(1.9255)
    assert entries[1].eps_kcal == pytest.approx(0.105)


def test_lj_to_sigma_eps_matches_sobtop_itp():
    # Sobtop itp 输出：UF_C sigma 3.430851E-01 nm、epsilon 4.393200E-01 kJ/mol
    sigma, eps = lj_to_sigma_eps(1.9255, 0.105)
    assert sigma == pytest.approx(3.430851e-1, abs=1e-7)
    assert eps == pytest.approx(4.393200e-1, abs=1e-9)


def test_parse_bonded_param_sections(bonded_file):
    lib = parse_bonded_param(bonded_file)
    assert len(lib.bonds) == 2
    assert len(lib.angles) == 2
    assert len(lib.dihedrals) == 4
    assert len(lib.impropers) == 3
    assert lib.bonds[0].types == ("c2", "c2")
    assert lib.bonds[0].k_kcal_a2 == pytest.approx(569.4)
    assert lib.bonds[0].r0_a == pytest.approx(1.3343)


def test_find_bond_exact_and_reverse(bonded_file):
    lib = parse_bonded_param(bonded_file)
    assert find_bond(lib, "c2", "ha").k_kcal_a2 == pytest.approx(343.0)
    assert find_bond(lib, "ha", "c2").k_kcal_a2 == pytest.approx(343.0)
    assert find_bond(lib, "oh", "ho") is None


def test_find_angle(bonded_file):
    lib = parse_bonded_param(bonded_file)
    entry = find_angle(lib, "ha", "c2", "ha")
    assert entry.theta0_deg == pytest.approx(116.900)
    assert find_angle(lib, "c2", "c2", "ha").theta0_deg == pytest.approx(120.430)


def test_find_dihedrals_specific_beats_wildcard(bonded_file):
    lib = parse_bonded_param(bonded_file)
    # 乙醇混合模式 itp：oh-c3-c3-hc（反序）命中 hc-c3-c3-oh 两条，而非 X-c3-c3-X
    terms = find_dihedrals(lib, "oh", "c3", "c3", "hc")
    assert [(t.kd_kj, t.phase_deg, t.pn) for t in terms] == [
        (pytest.approx(0.0), 0.0, 3.0),
        (pytest.approx(1.046, abs=1e-4), 0.0, 1.0),  # 0.25 kcal × 4.184
    ]


def test_find_dihedrals_idivf_and_wildcard(bonded_file):
    lib = parse_bonded_param(bonded_file)
    # hc-c3-c3-hc 无具体条目 → X -c3-c3-X（IDIVF=9）：1.400/9×4.184 = 0.65084 kJ
    terms = find_dihedrals(lib, "hc", "c3", "c3", "hc")
    assert len(terms) == 1
    assert terms[0].kd_kj == pytest.approx(0.650844, abs=1e-6)
    assert terms[0].pn == 3.0
    # 乙烯 X -c2-c2-X（IDIVF=4）：26.6/4×4.184 = 27.8236 kJ、phase 180、pn 2
    terms = find_dihedrals(lib, "ha", "c2", "c2", "ha")
    assert terms[0].kd_kj == pytest.approx(27.82360, abs=1e-5)
    assert terms[0].phase_deg == 180.0
    assert terms[0].pn == 2.0


def test_find_dihedrals_negative_pn(bonded_file):
    lib = parse_bonded_param(bonded_file)
    # 'hc-c3-c3-oh 1 0.00 0.0 -3.'：cos 偶函数 → pn=3、phase 0（Sobtop itp 实证）
    terms = find_dihedrals(lib, "oh", "c3", "c3", "hc")
    assert terms[0].pn == 3.0
    assert terms[0].phase_deg == 0.0


def test_improper_center_types(bonded_file):
    lib = parse_bonded_param(bonded_file)
    # 第三列（中心）类型集合：X-o-c-o → c；X-X-ca-ha → ca；X-X-c-o → c
    assert improper_center_types(lib.impropers) == {"c", "ca"}


def test_find_improper(bonded_file):
    lib = parse_bonded_param(bonded_file)
    # 苯环 improper：中心 ca、邻居 (ca, ca, ha) → X-X-ca-ha
    term = find_improper(lib, "ca", "ca", "ca", "ha")
    assert term is not None
    assert term.kd_kj == pytest.approx(1.1 * 4.184)
    assert term.phase_deg == 180.0
    assert term.pn == 2.0
    # 反序同样命中
    assert find_improper(lib, "ha", "ca", "ca", "ca") is not None
    # 中心不匹配 → None（调用方用 GUESS_IMPROPER_* 兜底）
    assert find_improper(lib, "ho", "oh", "c3", "c3") is None


def test_find_lj(tmp_path):
    p = tmp_path / "LJ_param.dat"
    p.write_text(LJ_PARAM)
    entries = parse_lj_param(p)
    assert find_lj(entries, "UF_Fe") is None
    assert find_lj(entries, " UF_C ").name == "UF_C"


def test_parse_bondcrit(tmp_path):
    p = tmp_path / "bondcrit.dat"
    p.write_text(BONDCRIT)
    crits = parse_bondcrit(p)
    assert len(crits) == 2
    assert (crits[0].element1, crits[0].element2) == ("Fe", "Cl")
    assert crits[0].lower_a == pytest.approx(2.0)
    assert crits[0].upper_a == pytest.approx(2.4)


def test_guess_constants_match_sobtop():
    # FeCl2 lazy itp：缺键 k=2.5e5 kJ/mol/nm²；缺角 k=500 kJ/mol/rad²；
    # 缺 improper 用 X-X-ca-ha 等价 1.1 kcal / 180° / pn 2
    assert GUESS_BOND_K_KJ_MOL_NM2 == 2.5e5
    assert GUESS_ANGLE_K_KJ_MOL_RAD2 == 500.0
    assert GUESS_IMPROPER_BARRIER_KCAL == pytest.approx(1.1)
    assert GUESS_IMPROPER_PN == 2.0
