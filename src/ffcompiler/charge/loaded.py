"""载入电荷：保留 IR 现有电荷，或从 mol2 / 文本文件读（MVP-3 计划 §2.3）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ffcompiler.bonded.assign import TopologyMismatchError
from ffcompiler.charge.base import ChargeProvider


def read_mol2_atoms(path: Path) -> list[tuple[str, float]]:
    """mol2 的 @<TRIPOS>ATOM 段 → [(atom name, 电荷)]。"""
    rows: list[tuple[str, float]] = []
    in_atoms = False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            continue
        if in_atoms and line.strip():
            parts = line.split()
            rows.append((parts[1], float(parts[8])))
    return rows


class LoadedCharges(ChargeProvider):
    method = "loaded"

    def __init__(self, source: Path | None = None):
        self.source = Path(source) if source is not None else None

    def compute(self, ir, positions_nm, mol) -> np.ndarray:
        n = len(ir.atoms)
        if self.source is None:
            return np.array([a.partial_charge for a in ir.atoms], dtype=np.float64)
        if self.source.suffix.lower() == ".mol2":
            rows = read_mol2_atoms(self.source)
            if len(rows) != n:
                raise TopologyMismatchError(f"{self.source}: {len(rows)} atoms, IR has {n}")
            for idx, ((name, _), atom) in enumerate(zip(rows, ir.atoms)):
                if name != atom.name:
                    raise TopologyMismatchError(
                        f"{self.source}: atom {idx} is named {name!r}, IR has {atom.name!r}")
            return np.array([q for _, q in rows], dtype=np.float64)
        values = [float(line) for line in self.source.read_text(encoding="utf-8").splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
        if len(values) != n:
            raise TopologyMismatchError(f"{self.source}: {len(values)} charges, IR has {n} atoms")
        return np.array(values, dtype=np.float64)

    def source_info(self) -> dict:
        return {"source": "ir" if self.source is None else str(self.source)}
