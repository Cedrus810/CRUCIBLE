"""ForceFieldIR 的 JSON 读写（MVP-0 计划 T1）。

浮点数必须逐位精确往返：``json`` 模块序列化 float 用 ``repr``（最短往返表示），
不加任何格式化截断。
"""

from __future__ import annotations

import json
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

from .schema import (
    SCHEMA_VERSION,
    Angle,
    Atom,
    Bond,
    ForceFieldIR,
    ImproperTorsion,
    ProperTorsion,
)


def _term_to_dict(term: Any) -> dict[str, Any]:
    d = asdict(term)
    if "meta" in d:
        d["meta"] = dict(d["meta"])
    return d


def ir_to_dict(ir: ForceFieldIR) -> dict[str, Any]:
    return {
        "name": ir.name,
        "atoms": [_term_to_dict(a) for a in ir.atoms],
        "bonds": [_term_to_dict(b) for b in ir.bonds],
        "angles": [_term_to_dict(a) for a in ir.angles],
        "proper_torsions": [_term_to_dict(t) for t in ir.proper_torsions],
        "impropers": [_term_to_dict(t) for t in ir.impropers],
        "exclusions": [list(p) for p in ir.exclusions],
        "pairs14": [list(p) for p in ir.pairs14],
        "total_charge": ir.total_charge,
        "forcefield_family": ir.forcefield_family,
        "combination_rule": ir.combination_rule,
        "coulomb_14_scale": ir.coulomb_14_scale,
        "vdw_14_scale": ir.vdw_14_scale,
        "provenance": ir.provenance,
    }


def _term_from_dict(cls: type, d: dict[str, Any]) -> Any:
    d = dict(d)
    names = {f.name for f in fields(cls)}
    if "meta" in names:
        d["meta"] = dict(d.get("meta") or {})
    else:
        d.pop("meta", None)
    return cls(**d)


def ir_from_dict(data: dict[str, Any]) -> ForceFieldIR:
    schema = data.get("schema")
    if schema != SCHEMA_VERSION:
        raise ValueError(
            f"ForceFieldIR schema version mismatch: file has {schema!r}, "
            f"this code expects {SCHEMA_VERSION!r}"
        )
    return ForceFieldIR(
        name=data["name"],
        atoms=tuple(_term_from_dict(Atom, a) for a in data["atoms"]),
        bonds=tuple(_term_from_dict(Bond, b) for b in data["bonds"]),
        angles=tuple(_term_from_dict(Angle, a) for a in data["angles"]),
        proper_torsions=tuple(_term_from_dict(ProperTorsion, t) for t in data["proper_torsions"]),
        impropers=tuple(_term_from_dict(ImproperTorsion, t) for t in data["impropers"]),
        exclusions=tuple(tuple(p) for p in data["exclusions"]),
        pairs14=tuple(tuple(p) for p in data["pairs14"]),
        total_charge=data["total_charge"],
        forcefield_family=data["forcefield_family"],
        combination_rule=data["combination_rule"],
        coulomb_14_scale=data["coulomb_14_scale"],
        vdw_14_scale=data["vdw_14_scale"],
        provenance=dict(data.get("provenance") or {}),
    )


def to_json(ir: ForceFieldIR, path: Path) -> None:
    payload = {"schema": SCHEMA_VERSION, **ir_to_dict(ir)}
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    Path(path).write_text(text + "\n", encoding="utf-8")


def from_json(path: Path) -> ForceFieldIR:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ir_from_dict(data)
