"""电荷写回 IR（MVP-3 计划 §2.2）。"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from ffcompiler.charge.base import ChargeError, ChargeProvider, ChargeReport, equalize_residual
from ffcompiler.ir import ForceFieldIR, validate


def assign_charges(
    ir: ForceFieldIR,
    provider: ChargeProvider,
    positions_nm: np.ndarray,
    mol,
    *,
    scan_cache: Path | None = None,
) -> tuple[ForceFieldIR, ChargeReport]:
    """provider 算电荷 → 均摊残差 → 只改 partial_charge，写 provenance["charge"]；
    IR 有 MLP 二面角时用 scan_cache 重拟合（缺缓存 → TorsionRefitError）。
    provider 返回的形状/有限性在这里统一检查；provenance["charge"]["previous"] 保留上一次的电荷来源。"""
    raw = np.asarray(provider.compute(ir, positions_nm, mol), dtype=np.float64)
    n = len(ir.atoms)
    if raw.shape != (n,):
        raise ChargeError(f"{provider.method}: provider 返回的电荷形状为 {raw.shape}，IR 有 {n} 个原子")
    if not np.all(np.isfinite(raw)):
        raise ChargeError(f"{provider.method}: provider 返回了非有限电荷（NaN/Inf）")
    charges, residual = equalize_residual(raw, ir.total_charge)
    report = ChargeReport(method=provider.method, raw_sum=float(raw.sum()), residual=residual,
                          source=provider.source_info())
    atoms = tuple(replace(a, partial_charge=float(q)) for a, q in zip(ir.atoms, charges, strict=True))
    prov = {"method": report.method, "raw_sum": report.raw_sum, "residual": report.residual,
            "source": report.source, "previous": ir.provenance.get("charge")}
    new_ir = ir.replace_terms(atoms=atoms, provenance_update={"charge": prov})
    validate(new_ir)
    if any(t.source == "mlp" for t in new_ir.proper_torsions):
        # 决策 0018：先电荷后二面角；已有 MLP 二面角就用缓存的扫描重拟合
        from ffcompiler.torsion.mlp import refit_torsions

        new_ir, refit = refit_torsions(new_ir, positions_nm, mol, scan_cache=scan_cache)
        report.torsion_refit = refit.to_dict()
    return new_ir, report
