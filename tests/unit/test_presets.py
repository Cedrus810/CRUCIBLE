"""presets.py：§17 表逐字冻结 + 编排语义的边界（2026-11-11 计划 §4）。"""

from __future__ import annotations

from ffcompiler.presets import PRESETS


def test_preset_table_matches_design_17_verbatim():
    assert PRESETS == {
        "fast": {
            "typing": "gaff2",
            "bonded": {"provider": "hessian", "hessian": "gfnff", "method": "mseminario"},
            "torsion": "mlp",
            "charge": "new_charge",
            "vdw": "gaff2",
        },
        "balanced": {
            "typing": "gaff2",
            "geometry": {"provider": "gfn2"},
            "bonded": {"provider": "hessian", "hessian": "gfn2", "method": "mseminario"},
            "torsion": "mlp",
            "charge": "new_charge",
            "vdw": "gaff2",
        },
        "high": {
            "geometry": {"provider": "orca"},
            "bonded": {"provider": "hessian", "hessian": "orca", "method": "mseminario"},
            "torsion": "mlp",
            "charge": "new_charge",
            "vdw": "gaff2",
        },
        "rigid": {
            "bonded": {"hessian": "gfn2", "method": "m2seminario"},
            "rigid_torsion": {"method": "drih"},
            "charge": "new_charge",
        },
    }


def test_rigid_preset_has_no_torsion_step():
    assert "torsion" not in PRESETS["rigid"]
    assert "rigid_torsion" in PRESETS["rigid"]  # 设计表不一致，原样保留待澄清


def test_all_default_charges_are_new_charge():
    for name, cfg in PRESETS.items():
        assert cfg["charge"] == "new_charge", name
