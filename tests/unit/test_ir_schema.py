"""T1：IR schema 与 JSON 读写的必测用例。"""

from __future__ import annotations

import json
import math
import struct
from dataclasses import asdict, FrozenInstanceError
from pathlib import Path

import pytest

from ffcompiler.ir import io as ir_io
from ffcompiler.ir import (
    SCHEMA_VERSION,
    Angle,
    Atom,
    Bond,
    ForceFieldIR,
    ImproperTorsion,
    ProperTorsion,
)


def _atoms4() -> tuple[Atom, ...]:
    return tuple(
        Atom(
            index=i,
            name=f"A{i}",
            element="C",
            atomic_number=6,
            mass=12.011 + 0.001 * i,
            atom_type="c3",
            formal_charge=0,
            partial_charge=-0.15 + 0.1 * i,
            sigma=0.3 + 0.01 * i,
            epsilon=0.4 + 0.02 * i,
        )
        for i in range(4)
    )


def make_4atom_ir() -> ForceFieldIR:
    propers = tuple(
        ProperTorsion(
            i=0, j=1, k=2, l=3, periodicity=n, phase=0.0 + 0.1 * n,
            force_constant=1.5 - 0.2 * n, source="test",
        )
        for n in (1, 2, 3)
    )
    return ForceFieldIR(
        name="t4",
        atoms=_atoms4(),
        bonds=(Bond(i=0, j=1, r0=0.152, k=250000.0, source="test"),),
        angles=(Angle(i=0, j=1, k_atom=2, theta0=1.9, force_constant=500.0, source="test"),),
        proper_torsions=propers,
        impropers=(
            ImproperTorsion(
                i=1, j=0, k=2, l=3, periodicity=2, phase=math.pi,
                force_constant=4.0, source="test",
            ),
        ),
        exclusions=((0, 1), (0, 2), (1, 2), (1, 3), (2, 3)),
        pairs14=((0, 3),),
        total_charge=0,
        forcefield_family="amber",
        combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2,
        vdw_14_scale=0.5,
        provenance={"tool": "test"},
    )


def _float_bits(value: float) -> str:
    return struct.pack(">d", value).hex()


def _assert_floats_bitwise(a: object, b: object, path: str = "") -> None:
    if isinstance(a, float) or isinstance(b, float):
        assert isinstance(a, float) and isinstance(b, float), path
        assert _float_bits(a) == _float_bits(b), f"float differs at {path}: {a!r} vs {b!r}"
    elif isinstance(a, dict):
        assert set(a) == set(b), path
        for key in a:
            _assert_floats_bitwise(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), path
        for n, (x, y) in enumerate(zip(a, b)):
            _assert_floats_bitwise(x, y, f"{path}[{n}]")


def test_json_roundtrip_exact(tmp_path: Path) -> None:
    ir = make_4atom_ir()
    assert len(ir.proper_torsions) == 3  # 同一四元组 n=1,2,3
    path = tmp_path / "ir.json"
    ir_io.to_json(ir, path)
    ir2 = ir_io.from_json(path)

    assert ir2 == ir
    _assert_floats_bitwise(asdict(ir2), asdict(ir))
    assert isinstance(ir2.atoms, tuple) and isinstance(ir2.proper_torsions, tuple)


def test_charged_ir_roundtrip_bitwise(tmp_path: Path) -> None:
    charges = [-0.8123456789012345, 0.3123456789012345, -0.5000000000000001, 1.0]
    atoms = tuple(
        Atom(
            index=i, name=f"A{i}", element="C", atomic_number=6, mass=12.011,
            atom_type="c3", formal_charge=0, partial_charge=q,
            sigma=0.3, epsilon=0.4,
        )
        for i, q in enumerate(charges)
    )
    ir = ForceFieldIR(
        name="charged", atoms=atoms, bonds=(), angles=(),
        proper_torsions=(), impropers=(),
        exclusions=((0, 1), (0, 2), (1, 2), (1, 3), (2, 3)), pairs14=((0, 3),),
        total_charge=0, forcefield_family="amber", combination_rule="lorentz-berthelot",
        coulomb_14_scale=1.0 / 1.2, vdw_14_scale=0.5,
    )
    path = tmp_path / "charged.json"
    ir_io.to_json(ir, path)
    ir2 = ir_io.from_json(path)

    for a, b in zip(ir.atoms, ir2.atoms):
        assert _float_bits(a.partial_charge) == _float_bits(b.partial_charge)
    assert ir2 == ir


def test_schema_version_mismatch_raises(tmp_path: Path) -> None:
    ir = make_4atom_ir()
    path = tmp_path / "ir.json"
    ir_io.to_json(ir, path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"] == SCHEMA_VERSION
    payload["schema"] = "ffir/0.2"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=r"ffir/0\.2.*ffir/0\.1|ffir/0\.1.*ffir/0\.2"):
        ir_io.from_json(path)


def test_replace_terms_returns_new_object() -> None:
    ir = make_4atom_ir()
    new_bond = Bond(i=0, j=1, r0=0.153, k=260000.0, source="replaced")
    new_proper = ProperTorsion(
        i=0, j=1, k=2, l=3, periodicity=1, phase=0.5, force_constant=2.0, source="replaced",
    )

    ir2 = ir.replace_terms(
        bonds=(new_bond,), proper_torsions=(new_proper,),
        provenance_update={"bonded": "mvp1"},
    )

    assert ir2.bonds == (new_bond,)
    assert ir2.proper_torsions == (new_proper,)
    assert ir2.provenance == {"tool": "test", "bonded": "mvp1"}
    assert ir2.angles == ir.angles and ir2.impropers == ir.impropers
    assert ir2.atoms == ir.atoms
    # 原对象完全不变
    assert ir.bonds == (Bond(i=0, j=1, r0=0.152, k=250000.0, source="test"),)
    assert ir.provenance == {"tool": "test"}
    assert ir2 is not ir


def test_dataclasses_are_frozen() -> None:
    ir = make_4atom_ir()
    with pytest.raises(FrozenInstanceError):
        ir.atoms[0].partial_charge = 1.0
    with pytest.raises(FrozenInstanceError):
        ir.bonds[0].k = 0.0
    with pytest.raises(FrozenInstanceError):
        ir.name = "other"
