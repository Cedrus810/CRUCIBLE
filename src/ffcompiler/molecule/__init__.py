"""分子层：SDF 入口的读入与检查（设计 §3.1；MVP-6）。"""

from .io import CLASH_MIN_A, MoleculeInput, MoleculeInputError, read_sdf

__all__ = ["CLASH_MIN_A", "MoleculeInput", "MoleculeInputError", "read_sdf"]
