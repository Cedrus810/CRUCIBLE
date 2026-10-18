"""坐标 sidecar（ffcompiler/ir/coords.py，C1；ffcoords/0.2 带原子名）的必测用例。

纯 Python/JSON，无外部依赖；接收方在计算节点先跑本文件再跑管线串联。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ffcompiler.ir import coords as ir_coords
from ffcompiler.ir.coords import Coordinates, CoordinatesMismatchError


def _mini_ir_like_atoms():
    # (atomic_number, element, name)
    return ((6, "C", "C1"), (1, "H", "H1"), (1, "H", "H2"))


def _bundle(positions, order=None) -> Coordinates:
    return Coordinates(
        positions_nm=np.asarray(positions, dtype=np.float64),
        source="test",
        atom_order=order or _mini_ir_like_atoms(),
    )


class _AtomLike:
    def __init__(self, z, el, name):
        self.atomic_number, self.element, self.name = z, el, name


class _IRLike:
    def __init__(self, atoms):
        self.atoms = atoms


def test_json_roundtrip_float_exact(tmp_path: Path) -> None:
    pos = np.array([[0.1234567890123456789, -1e-17, 3.141592653589793],
                    [1e-300, -2.5e10, 1.0 / 3.0], [0.0, -0.0, 7.0]])
    path = tmp_path / "c.json"
    ir_coords.to_json(_bundle(pos), path)
    back = ir_coords.from_json(path)
    assert back.source == "test"
    assert back.atom_order == _mini_ir_like_atoms()
    assert back.positions_nm.tolist() == pos.tolist()  # 逐位一致（repr 最短往返）


def test_schema_version_mismatch_rejected(tmp_path: Path) -> None:
    path = tmp_path / "c.json"
    ir_coords.to_json(_bundle(np.zeros((3, 3))), path)
    import json

    data = json.loads(path.read_text())
    data["schema"] = "ffcoords/9.9"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="schema"):
        ir_coords.from_json(path)


def test_units_must_be_nm(tmp_path: Path) -> None:
    import json

    path = tmp_path / "c.json"
    ir_coords.to_json(_bundle(np.zeros((3, 3))), path)
    data = json.loads(path.read_text())
    data["units"] = "angstrom"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="nm"):
        ir_coords.from_json(path)


def test_check_against_ir_accepts_matching_order() -> None:
    ir = _IRLike([_AtomLike(*a) for a in _mini_ir_like_atoms()])
    ir_coords.check_against_ir(_bundle(np.zeros((3, 3))), ir)


def test_check_against_ir_rejects_permutation_or_wrong_element() -> None:
    ir = _IRLike([_AtomLike(*a) for a in _mini_ir_like_atoms()])
    permuted = (_mini_ir_like_atoms()[1], _mini_ir_like_atoms()[0],
                _mini_ir_like_atoms()[2])
    with pytest.raises(CoordinatesMismatchError, match="原子顺序不一致"):
        ir_coords.check_against_ir(_bundle(np.zeros((3, 3)), order=permuted), ir)
    wrong_el = ((6, "C", "C1"), (1, "H", "H1"), (8, "O", "H2"))
    with pytest.raises(CoordinatesMismatchError):
        ir_coords.check_against_ir(_bundle(np.zeros((3, 3)), order=wrong_el), ir)
    # 同元素原子换序（H1/H2 对调）：元素序列相同，只有原子名能查出
    same_element_swap = ((6, "C", "C1"), (1, "H", "H2"), (1, "H", "H1"))
    with pytest.raises(CoordinatesMismatchError, match="原子顺序不一致"):
        ir_coords.check_against_ir(_bundle(np.zeros((3, 3)), order=same_element_swap), ir)


def test_bad_shapes_rejected() -> None:
    with pytest.raises(ValueError, match="N,3"):
        _bundle(np.zeros(3))
    with pytest.raises(ValueError, match="atom_order"):
        _bundle(np.zeros((3, 3)), order=((6, "C", "C1"), (1, "H", "H1")))


def test_from_positions_takes_order_from_ir() -> None:
    ir = _IRLike([_AtomLike(*a) for a in _mini_ir_like_atoms()])
    bundle = ir_coords.from_positions(np.zeros((3, 3)), ir, "hessian:gfn2")
    assert bundle.source == "hessian:gfn2"
    assert bundle.atom_order == _mini_ir_like_atoms()
