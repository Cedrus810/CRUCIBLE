"""原子 typing 层（设计 §3.1/§21；MVP-6）：GAFF2 typing 由外部 antechamber 承担。"""

from .gaff2 import EntryResult, TypingError, ambertools_version, equalize_mol2_charges, run_entry

__all__ = ["EntryResult", "TypingError", "ambertools_version", "equalize_mol2_charges", "run_entry"]
