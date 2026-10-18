"""MLPTorsionProvider（MVP-2 计划 §2.3–§2.7，T28）。

每根可转键：约束松弛扫描 → y = E_ref − E_MM[IR \\ b] → Fourier 拟合 → 替换该键上的 proper torsion。
失败退回 GAFF2（保留原项），逐键记录原因（§2.6）。
"""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Mapping

import numpy as np
from rdkit import Chem

from ..ir.schema import ForceFieldIR, ProperTorsion
from ..ir.validation import validate
from .base import TorsionProvider
from .fourier import fit_fourier, to_ir_terms
from .mace import MODEL_NAME, molecule_info
from .mm import mm_energy, mm_energy_without_bond
from .rotors import find_rotors
from .scan import ScanError, ScanResult, relaxed_scan

SCAN_DEFAULTS = {"step_deg": 10.0, "fmax": 0.005, "max_steps": 2000}


class TorsionRefitError(RuntimeError):
    """换电荷后无法重拟合 MLP 二面角（缺缓存、缺 replaced_terms 记录）；不静默退回 GAFF2。"""


@dataclass
class TorsionReport:
    rotors: list[dict]
    molecule_fallback: str | None
    scans: dict[tuple[int, int], ScanResult] = field(default_factory=dict, repr=False)
    refinement: dict | None = None

    def to_dict(self) -> dict:
        out = {"molecule_fallback": self.molecule_fallback, "rotors": self.rotors}
        if self.refinement is not None:
            out["refinement"] = self.refinement
        return out


def _validate_refinement_rounds(rounds: int) -> None:
    if type(rounds) is not int or not 0 <= rounds <= 3:
        raise ValueError("refinement_rounds must be an integer from 0 to 3 (extra sweeps after single-pass)")


def _round_metrics(ir: ForceFieldIR, rows: list[dict], scans: Mapping) -> dict:
    """同一轮末 IR 的样本内剖面；每条扫描独立平移，逐帧等权。"""
    bonds, residuals = [], []
    for row in rows:
        bond = tuple(row["bond"])
        scan = scans[bond]
        energy = np.asarray(mm_energy(ir, scan.positions), dtype=np.float64)
        delta = energy - scan.energy
        if not np.all(np.isfinite(delta)):
            raise TorsionRefitError(f"non-finite refinement profile for bond {bond}")
        centered = delta - delta.mean()
        residuals.append(centered)
        min_shift = (energy - energy.min()) - (scan.energy - scan.energy.min())
        bonds.append({"bond": list(bond), "n_points": int(len(delta)),
                      "fit_rmse_at_update": row["fit_rmse"],
                      "profile_rmse_best_shift": float(np.sqrt(np.mean(centered ** 2))),
                      "profile_rmse_min_shift": float(np.sqrt(np.mean(min_shift ** 2))),
                      "max_hysteresis": row["max_hysteresis"],
                      "max_bond_change_pm": row["max_bond_change_pm"]})
    pooled = np.concatenate(residuals) if residuals else np.empty(0)
    return {"bonds": bonds, "n_frames": int(pooled.size),
            "pooled_profile_rmse": float(np.sqrt(np.mean(pooled ** 2))) if pooled.size else None}


def _default_factory(ir: ForceFieldIR):
    from .mace import polar_calculator

    return polar_calculator()


def _default_supported() -> frozenset[int]:
    from .mace import supported_atomic_numbers

    return supported_atomic_numbers()


PACKAGES = ("mace-torch", "graph_longrange", "ase", "torch")


def package_versions() -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    out = {}
    for name in PACKAGES:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = "not installed"
    return out


def model_file_identity(path: Path | None) -> dict | None:
    """模型文件身份：路径、大小、修改时间（不用 hash，全局规则 §8.4）。"""
    if path is None:
        return None
    st = Path(path).stat()
    return {"path": str(path), "size": st.st_size, "mtime": st.st_mtime}


def scan_cache_key(model: str, quartet, info: Mapping, scan_kwargs: Mapping, *,
                   model_file: Mapping | None, numbers, bonds, reference: Mapping | None = None) -> str:
    """扫描缓存键：模型名与文件身份、四元组、电荷/spin、扫描参数、原子序数与键（计划 §2.7）；
    reference：参考势依赖 IR 时（mm-self）的身份，逐项字面值。为 None 时不写入，已有 MACE 缓存键不变。"""
    payload = {"model": model, "model_file": dict(model_file) if model_file else None,
               "quartet": list(quartet), "info": dict(info),
               "scan": {k: scan_kwargs[k] for k in sorted(scan_kwargs)},
               "numbers": [int(z) for z in numbers],
               "bonds": sorted([min(int(a), int(b)), max(int(a), int(b))] for a, b in bonds)}
    if reference is not None:
        payload["reference"] = dict(reference)
    return json.dumps(payload, sort_keys=True)


def save_scan(path: Path, key: str, scan: ScanResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, key=np.array(key), quartet=np.array(scan.quartet), phi_target_rad=scan.phi_target_rad,
             phi_rad=scan.phi_rad, energy=scan.energy, positions=scan.positions,
             hysteresis=scan.hysteresis, converged=scan.converged,
             max_bond_change_pm=np.array(scan.max_bond_change_pm), params=np.array(json.dumps(scan.params)))
    tmp.replace(path)


def load_scan(path: Path, key: str) -> ScanResult | None:
    """文件存在且记录的扫描参数逐项相等才复用（不用 hash，计划 §2.7）。"""
    if not path.is_file():
        return None
    with np.load(path) as d:
        if str(d["key"]) != key:
            return None
        return ScanResult(
            quartet=tuple(int(a) for a in d["quartet"]), phi_target_rad=d["phi_target_rad"],
            phi_rad=d["phi_rad"], energy=d["energy"], positions=d["positions"],
            hysteresis=d["hysteresis"], converged=d["converged"].astype(bool),
            max_bond_change_pm=float(d["max_bond_change_pm"]), params=json.loads(str(d["params"])),
        )


def molecule_fallback_reason(ir: ForceFieldIR, supported: frozenset[int]) -> str | None:
    """整分子兜底原因（计划 §2.6）；None 表示可以扫描。"""
    if any(a.atomic_number not in supported for a in ir.atoms):
        return "unsupported_element"
    if (sum(a.atomic_number for a in ir.atoms) - ir.total_charge) % 2:
        return "odd_electrons"
    return None


class MLPTorsionProvider(TorsionProvider):
    def __init__(
        self,
        calculator_factory: Callable[[ForceFieldIR], object] | None = None,
        *,
        supported: frozenset[int] | None = None,
        scan_kwargs: Mapping | None = None,
        scan_cache: Path | None = None,
        model_name: str = MODEL_NAME,
        model_path: Path | None = None,
        cache_only: bool = False,
        reference: Mapping | None = None,
        fit_bonds: Iterable[Iterable[int]] | None = None,
        refinement_rounds: int = 0,
    ) -> None:
        """model_path：参考势的模型文件（写进 provenance 与缓存键）；默认工厂时取 default_model_path()。
        cache_only：只从缓存读扫描，缺失记为兜底 `no_cached_scan`，不在本进程重扫（Gate 拟合阶段用）。
        reference：参考势依赖 IR 时的身份（mm-self 记扫描时的电荷），进缓存键与 provenance；参考势在扫描时冻结。
        fit_bonds：只处理这些键（重拟合用）；其余键保留原项，记为兜底 `not_selected`。
        refinement_rounds：单遍之后额外做 0–3 轮顺序更新（默认 0），不重扫；
        每次仍仅使用本键扫描，不保证全部扫描的联合目标单调下降。"""
        _validate_refinement_rounds(refinement_rounds)
        self.refinement_rounds = refinement_rounds
        if model_path is None and calculator_factory is None:
            from .mace import default_model_path

            model_path = default_model_path()
        self.model_path = Path(model_path) if model_path is not None else None
        self.cache_only = cache_only
        self.calculator_factory = calculator_factory or _default_factory
        self._supported = supported
        self.scan_kwargs = {**SCAN_DEFAULTS, **dict(scan_kwargs or {})}
        self.scan_cache = Path(scan_cache) if scan_cache is not None else None
        self.model_name = model_name
        self.reference = dict(reference) if reference is not None else None
        self.fit_bonds = None if fit_bonds is None else {frozenset(b) for b in fit_bonds}

    def parameterize(self, ir: ForceFieldIR, positions_nm: np.ndarray, mol: Chem.Mol):
        pos = np.asarray(positions_nm, dtype=np.float64)
        rotors = sorted(find_rotors(ir, mol, pos), key=lambda rotor: rotor.bond)
        supported = self._supported if self._supported is not None else _default_supported()
        info = molecule_info(ir)
        rows = [{"bond": list(r.bond), "rotor_class": r.rotor_class,
                 "scan_quartet": list(r.scan_quartet), "status": None, "reason": None}
                for r in rotors]

        molecule_fallback = molecule_fallback_reason(ir, supported)
        if molecule_fallback:
            for row in rows:
                row.update(status="fallback", reason=molecule_fallback)
            report = TorsionReport(rows, molecule_fallback)
            torsions = list(ir.proper_torsions)
            if self.refinement_rounds:
                torsions = self._refine(ir, torsions, report)
            return self._finish(ir, torsions, report, info), report

        numbers = [a.atomic_number for a in ir.atoms]
        bonds = [(b.i, b.j) for b in ir.bonds]
        torsions = list(ir.proper_torsions)
        report = TorsionReport(rows, None)
        model_file = model_file_identity(self.model_path)
        for rotor, row in zip(rotors, rows):
            j, k = rotor.bond
            if self.fit_bonds is not None and frozenset(rotor.bond) not in self.fit_bonds:
                row.update(status="fallback", reason="not_selected")
                continue
            key = scan_cache_key(self.model_name, rotor.scan_quartet, info, self.scan_kwargs,
                                 model_file=model_file, numbers=numbers, bonds=bonds,
                                 reference=self.reference)
            cache = self.scan_cache / f"{j}-{k}.npz" if self.scan_cache else None
            scan = load_scan(cache, key) if cache else None
            if scan is None and self.cache_only:
                row.update(status="fallback", reason="no_cached_scan")
                continue
            if scan is None:
                try:
                    calc = self.calculator_factory(ir)
                    scan = relaxed_scan(calc, numbers, pos, rotor.scan_quartet, bonds,
                                        info=info, **self.scan_kwargs)
                except ScanError as exc:
                    row.update(status="fallback", reason="mlp_error", error=str(exc))
                    continue
                if cache:
                    save_scan(cache, key, scan)
            report.scans[rotor.bond] = scan
            row.update(max_hysteresis=float(np.max(scan.hysteresis)),
                       max_bond_change_pm=scan.max_bond_change_pm)
            if not scan.converged.all():
                row.update(status="fallback", reason="scan_not_converged")
                continue
            y = scan.energy - mm_energy_without_bond(ir, rotor.bond, scan.positions)
            fit = fit_fourier(scan.phi_rad, y)
            meta = {"model": self.model_name, "scan_quartet": "-".join(map(str, rotor.scan_quartet)),
                    "fit_rmse": fit.rmse, "n_points": int(len(y)), "rotor_class": rotor.rotor_class}
            torsions = [t for t in torsions if {t.j, t.k} != {j, k}]
            torsions += to_ir_terms(fit, rotor.scan_quartet, source="mlp", meta=meta)
            replaced = [{"index": idx, **asdict(t), "meta": dict(t.meta)}
                        for idx, t in enumerate(ir.proper_torsions) if {t.j, t.k} == {j, k}]
            row.update(status="mlp", fit_rmse=fit.rmse, fit_max_abs_residual=fit.max_abs_residual,
                       replaced_terms=replaced)
        if self.refinement_rounds:
            torsions = self._refine(ir, torsions, report)
        return self._finish(ir, torsions, report, info), report

    def _refine(self, ir: ForceFieldIR, torsions: list, report: TorsionReport) -> list:
        """只更新本次单遍成功的键，固定键号顺序；原 replaced_terms 不改写。"""
        rows = sorted((r for r in report.rotors if r["status"] == "mlp"),
                      key=lambda r: tuple(r["bond"]))
        current = ir.replace_terms(proper_torsions=tuple(torsions))
        history = [{"round": 0, **_round_metrics(current, rows, report.scans)}]
        for rnd in range(1, (self.refinement_rounds if rows else 0) + 1):
            for row in rows:
                bond = tuple(row["bond"])
                scan = report.scans[bond]
                y = scan.energy - mm_energy_without_bond(current, bond, scan.positions)
                fit = fit_fourier(scan.phi_rad, y)
                meta = {"model": self.model_name,
                        "scan_quartet": "-".join(map(str, scan.quartet)),
                        "fit_rmse": fit.rmse, "n_points": int(len(y)),
                        "rotor_class": row["rotor_class"], "refinement_round": rnd}
                kept = [t for t in current.proper_torsions if {t.j, t.k} != set(bond)]
                kept += to_ir_terms(fit, scan.quartet, source="mlp", meta=meta)
                current = current.replace_terms(proper_torsions=tuple(kept))
                row.update(fit_rmse=fit.rmse, fit_max_abs_residual=fit.max_abs_residual)
            history.append({"round": rnd, **_round_metrics(current, rows, report.scans)})
        report.refinement = {
            "extra_rounds": self.refinement_rounds,
            "completed_rounds": self.refinement_rounds if rows else 0,
            "objective": "per_bond_scan", "evaluation": "in_sample",
            "order": [list(row["bond"]) for row in rows], "history": history,
        }
        return list(current.proper_torsions)

    def _finish(self, ir: ForceFieldIR, torsions: list, report: TorsionReport, info: Mapping) -> ForceFieldIR:
        from ffcompiler.bonded.staleness import fixed_proper_snapshot, mm_subtraction_snapshot

        order = sorted([list(r["bond"]) for r in report.rotors if r["status"] == "mlp"])
        prov = {"provider": "mlp", "model": self.model_name, "info": dict(info),
                "scan": dict(self.scan_kwargs), "reference": deepcopy(self.reference), "model_file": model_file_identity(self.model_path),
                "fit_protocol": {"name": "sequential" if self.refinement_rounds else "single_pass",
                                 "refinement_rounds": self.refinement_rounds, "order": order},
                "charge_method": (ir.provenance.get("charge") or {}).get("method", "prmtop"),
                # C4（2026-10-15 计划 §5）：MM 扣除用的 bond/angle/电荷/LJ 快照，供 staleness 检查
                "mm_dependency": mm_subtraction_snapshot(ir),
                "package_versions": package_versions(), **deepcopy(report.to_dict())}
        if self.refinement_rounds:
            prov["fixed_proper_dependency"] = fixed_proper_snapshot(ir, order)
        new_ir = ir.replace_terms(proper_torsions=tuple(torsions), provenance_update={"torsion": prov})
        validate(new_ir)
        return new_ir


def _no_reference(ir: ForceFieldIR):
    raise TorsionRefitError("refit reads cached scans only; the reference potential must not be called")


def refit_torsions(ir: ForceFieldIR, positions_nm: np.ndarray, mol: Chem.Mol, *, scan_cache: Path | None,
                   model_path: Path | None = None, calculator_factory=None):
    """换电荷后重拟合 MLP 二面角（MVP-3 计划 §2.4）：还原被替换的原项（按原下标），
    再用缓存按记录的单遍/顺序更新协议重做；旧 IR 缺协议时沿用单遍。
    没有 mlp 项时原样返回。"""
    if not any(t.source == "mlp" for t in ir.proper_torsions):
        return ir, TorsionReport([], None)
    prov = ir.provenance.get("torsion") or {}
    protocol = prov.get("fit_protocol", {})
    if not isinstance(protocol, dict):
        raise TorsionRefitError("invalid torsion fit_protocol; re-run derive-torsions")
    rounds = protocol.get("refinement_rounds", 0)
    _validate_refinement_rounds(rounds)
    if protocol.get("name", "single_pass") != ("sequential" if rounds else "single_pass"):
        raise TorsionRefitError("unsupported torsion fit_protocol; re-run derive-torsions")
    rows = [r for r in prov.get("rotors", []) if r.get("status") == "mlp"]
    if not rows or any("replaced_terms" not in r for r in rows):
        raise TorsionRefitError("IR has MLP torsions but no replaced_terms record; "
                                "re-run derive-torsions from the GAFF2 IR (with --scan-cache)")
    if protocol and protocol.get("order") != sorted([list(r["bond"]) for r in rows]):
        raise TorsionRefitError("torsion fit_protocol order does not match the original mlp bonds")
    if scan_cache is None:
        raise TorsionRefitError("IR has MLP torsions: changing charges needs scan_cache to refit them")
    mlp_bonds = [frozenset(r["bond"]) for r in rows]
    terms = [t for t in ir.proper_torsions if frozenset((t.j, t.k)) not in mlp_bonds]
    for rec in sorted((rec for r in rows for rec in r["replaced_terms"]), key=lambda rec: rec["index"]):
        fields_ = {k: v for k, v in rec.items() if k != "index"}
        terms.insert(rec["index"], ProperTorsion(**{**fields_, "meta": dict(fields_.get("meta") or {})}))
    base = ir.replace_terms(proper_torsions=tuple(terms))
    if model_path is None and prov.get("model_file"):
        model_path = Path(prov["model_file"]["path"])
    provider = MLPTorsionProvider(calculator_factory or _no_reference, supported=frozenset(range(1, 119)),
                                  scan_kwargs=prov.get("scan"), scan_cache=scan_cache,
                                  model_name=prov.get("model", MODEL_NAME), model_path=model_path,
                                  cache_only=True, reference=prov.get("reference"), fit_bonds=mlp_bonds,
                                  refinement_rounds=rounds)
    new_ir, report = provider.parameterize(base, positions_nm, mol)
    lost = [r["bond"] for r in report.rotors if r["status"] != "mlp" and frozenset(r["bond"]) in mlp_bonds]
    if lost:
        raise TorsionRefitError("no cached scan for bond(s) " + ", ".join(f"{a}-{b}" for a, b in lost)
                                + f" in {scan_cache}")
    # 原来就不是 mlp 的键（mlp_error、scan_not_converged 等）只重拟合原 mlp 键，其余原样沿用原行与原项；
    # 缓存里后来补齐了扫描也不升级为 mlp（2026-10-13 用户定）
    original = {frozenset(r["bond"]): r for r in prov.get("rotors", [])}
    kept = [r if frozenset(r["bond"]) in mlp_bonds else deepcopy(original.get(frozenset(r["bond"]), r))
            for r in report.rotors]
    report.rotors = kept
    torsion_prov = {**new_ir.provenance["torsion"], "rotors": deepcopy(kept)}
    new_ir = new_ir.replace_terms(provenance_update={"torsion": torsion_prov})
    return new_ir, report
