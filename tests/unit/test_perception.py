"""perception（几何化学感知）单元测试。

期望值来自 Sobtop 2026.1.16 对相同几何的输出（沙箱 et/ethylene/benzene
itp 的键/角/二面/improper 列表逐条对照），见
docs/theory/sobtop_reverse_engineering.md 附录 A。
"""

from __future__ import annotations

import numpy as np
import pytest

from ffcompiler.perception import (
    COVALENT_RADII_A,
    build_angles,
    build_dihedrals,
    covalent_radius,
    find_improper_centers,
    nearest_image,
    perceive_bonds,
)

# 乙醇几何（et.mol2，Å → nm）；原子序 C,C,O,H×6
ET_POS_NM = np.array(
    [
        [0.894, -0.088, -0.105],
        [-0.578, -0.360, 0.127],
        [-1.236, 0.834, 0.522],
        [1.419, -1.003, -0.395],
        [1.358, 0.315, 0.800],
        [1.030, 0.660, -0.893],
        [-1.055, -0.740, -0.782],
        [-0.710, -1.101, 0.920],
        [-1.122, 1.483, -0.194],
    ]
) * 0.1
ET_Z = [6, 6, 8, 1, 1, 1, 1, 1, 1]
# mol2 BOND 块（0 基）
ET_BONDS = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (2, 8)]
# Sobtop itp [angles]（0 基）
ET_ANGLES = [
    (1, 0, 3), (1, 0, 4), (1, 0, 5), (3, 0, 4), (3, 0, 5), (4, 0, 5),
    (0, 1, 2), (0, 1, 6), (0, 1, 7), (2, 1, 6), (2, 1, 7), (6, 1, 7),
    (1, 2, 8),
]
# Sobtop itp [dihedrals]（0 基，规范 i<l）
ET_DIHEDRALS = sorted(
    [
        (0, 1, 2, 8),
        (2, 1, 0, 3), (2, 1, 0, 4), (2, 1, 0, 5),
        (3, 0, 1, 6), (3, 0, 1, 7),
        (4, 0, 1, 6), (4, 0, 1, 7),
        (5, 0, 1, 6), (5, 0, 1, 7),
        (6, 1, 2, 8), (7, 1, 2, 8),
    ]
)


def test_covalent_radii_table():
    # 从 Sobtop .rodata 抠出的 Cordero 2008 值
    assert COVALENT_RADII_A[1] == pytest.approx(0.31)
    assert COVALENT_RADII_A[6] == pytest.approx(0.76)
    assert COVALENT_RADII_A[8] == pytest.approx(0.66)
    assert COVALENT_RADII_A[17] == pytest.approx(1.02)
    assert COVALENT_RADII_A[26] == pytest.approx(1.32)
    assert covalent_radius(0) == 0.0


def test_perceive_bonds_ethanol():
    bonds = perceive_bonds(ET_Z, ET_POS_NM)
    assert bonds == sorted(ET_BONDS)


def test_build_angles_ethanol_matches_sobtop():
    angles = build_angles(ET_BONDS, 9)
    assert angles == ET_ANGLES
    assert len(angles) == 13


def test_build_dihedrals_ethanol_matches_sobtop():
    dih = build_dihedrals(ET_BONDS, 9, ET_POS_NM)
    assert sorted(dih) == ET_DIHEDRALS
    assert len(dih) == 12  # 甲基-甲基 9 条全出（含氢、不去重）


def test_linear_flanking_angle_cuts_dihedral():
    # CH3–C≡N：C1-C2-N 共线（180°）→ 绕 C1-C2 的二面角全部剔除；
    # 但线性角本身仍进入角列表（Sobtop 语义）
    pos = np.array(
        [
            [0.0, 0.0, 0.0],      # 0 C1
            [0.15, 0.0, 0.0],     # 1 C2
            [0.30, 0.0, 0.0],     # 2 N（共线）
            [0.05, 0.1, 0.05],    # 3 H
            [0.05, -0.05, 0.1],   # 4 H
            [0.05, -0.05, -0.1],  # 5 H
        ]
    )
    bonds = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5)]
    dih = build_dihedrals(bonds, 6, pos)
    assert dih == []
    angles = build_angles(bonds, 6)
    assert (0, 1, 2) in angles  # 线性角保留


def test_find_improper_centers_benzene():
    # 苯（1 基 itp: (2,6,1,7)(1,3,2,8)(2,4,3,9)(3,5,4,10)(4,6,5,11)(1,5,6,12)）
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
             (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)]
    types = ["ca"] * 6 + ["ha"] * 6
    out = find_improper_centers(types, bonds, {"ca", "c", "c2", "n", "na", "n2"})
    assert out == [
        (1, 5, 0, 6), (0, 2, 1, 7), (1, 3, 2, 8),
        (2, 4, 3, 9), (3, 5, 4, 10), (0, 4, 5, 11),
    ]


def test_find_improper_centers_requires_type_in_set():
    # 乙烯拓扑：UFF 类型（UF_C 不在集合）→ 无 improper（Sobtop 实证）；
    # GAFF c2 在集合 → 两条
    bonds = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
    assert find_improper_centers(["UF_C"] * 6, bonds, {"ca"}) == []
    out = find_improper_centers(["c2", "c2", "ha", "ha", "ha", "ha"], bonds, {"c2"})
    assert out == [(1, 2, 0, 3), (0, 4, 1, 5)]  # itp: (2,3,1,4)(1,5,2,6)


def test_nearest_image_orthorhombic():
    box = np.diag([2.0, 2.0, 2.0])
    d = nearest_image(np.array([1.9, 0.3, -1.9]), box)
    assert d == pytest.approx(np.array([-0.1, 0.3, 0.1]))


def test_nearest_image_triclinic():
    # FeCl2 六方晶胞：a=(1.7023,0,0) b=(-0.85115,1.47424,0) c=(0,0,2.2461)
    box = np.array(
        [[1.70230, 0.0, 0.0], [-0.85115, 1.47424, 0.0], [0.0, 0.0, 2.24610]]
    )
    d = nearest_image(box[0] + box[1], box)  # 格点向量折回原点
    assert np.allclose(d, 0.0, atol=1e-12)


def test_perceive_bonds_periodic_mirror():
    # 孤立体系下相距 1.96 nm 的两个 H 不成键；2 nm 正交盒的最近镜像下相距
    # 0.04 nm < (0.31+0.31)×1.15 Å → 成键
    pos = np.array([[0.02, 0.5, 0.5], [1.98, 0.5, 0.5]])
    assert perceive_bonds([1, 1], pos) == []
    box = np.diag([2.0, 2.0, 2.0])
    assert perceive_bonds([1, 1], pos, box_nm=box) == [(0, 1)]


def test_perceive_bonds_pair_window_override():
    # 半径判据 (1.32+1.02)×1.15 = 2.691 Å：2.75 Å 的 Fe-Cl 不成键；
    # bondcrit 窗口 [2.6, 2.8] 覆盖 → 成键（窗口优先于半径判据的并集语义）
    pos = np.array([[0.0, 0.0, 0.0], [0.275, 0.0, 0.0]])
    assert perceive_bonds([26, 17], pos) == []
    windows = {(26, 17): (2.6, 2.8)}
    assert perceive_bonds([26, 17], pos, pair_windows=windows) == [(0, 1)]
