# Changelog

All notable changes to ffcompiler. Dates are Asia/Tokyo. Decision records live in
[`docs/decisions/`](docs/decisions/README.md).

## [0.2.0] — 2026-10-05 — MVP-1: Hessian-derived bond/angle parameters

Gate G0 and Gate G1 both pass (`jobs/gate_mvp0_20261005_3/`, `jobs/gate_mvp1_20261005_4/`).

### Added
- **Hessian subsystem** (`ffcompiler.hessian`): `CartesianHessian` container (kJ mol⁻¹ nm⁻²),
  quality checks (asymmetry, translation residual, stationarity, imaginary modes), harmonic
  frequencies, Gaussian fchk reader, xtb 6.7.1 runner/reader (GFN-FF, GFN2).
- **Bond/angle derivation** (`ffcompiler.bonded`): mSeminario (Allen 2018 eigen-projection with
  adjacent-angle scaling), m2Seminario (candidate C2), DRIH (redundant internal coordinates in
  atomic units), global Hessian fit with H_fixed subtraction (MVP-1b), graph-automorphism
  symmetrization, `derive_bonded` driver with atom-order/topology checks.
- **CLI**: `ffcompiler derive-bonded` and `ffcompiler --version`.
- **Gate G1** (`scripts/gate_mvp1.py`): 11 molecules × {GFN-FF, GFN2} × {gaff2, mseminario,
  m2seminario, drih, global_fit}; acceptance derived criterion by criterion (G1.1–G1.7).
- **Native AMBER I/O** (`ffcompiler.formats.amber`): prmtop / rst7 reader and writer.
- Sobtop oracle data (12 cases) and `scripts/sobtop_oracle.py`.

### Changed
- AMBER importer and exporter no longer depend on ParmEd (decision 0015). The importer now
  rejects non-Lorentz–Berthelot LJ tables and takes elements only from `ATOMIC_NUMBER`.
- Mixed-precision GROMACS is judged with an fp32-aware tolerance `1e-4 + 1e-5·|E|` kJ/mol;
  double builds keep the strict file tolerance (decision 0013, supersedes 0006).
- Acetonitrile (near-linear angle) is a diagnostic case, not a Gate fixture (decision 0013).
- Linear angles (> 175°) raise `LinearAngleError` by default; `library` policy keeps GAFF2 values.

### Fixed
- `derive_bonded` reports non-positive force constants with `NonPositiveForceConstantError`;
  Sobtop-confirmed cases are recorded as known method limitations (decision 0014).
- `RunContext` probes git with `--no-optional-locks` so a killed `git status` cannot leave
  `index.lock` behind.
- Theory document: AMBER conversion now states the factor 2 (K = k/2); SF6 check value
  corrected to 350.02 kcal mol⁻¹ Å⁻².

## [0.1.0] — 2026-10-04 — MVP-0: engine-independent IR and three exporters

### Added
- `ForceFieldIR` schema, validation, NumPy reference energies, AMBER prmtop importer (GAFF2).
- Exporters for OpenMM (System/XML), AMBER (prmtop/inpcrd) and GROMACS (top/itp/gro).
- Cross-engine comparison (`numpy-ir`, OpenMM in-memory/XML/prmtop/top, sander, GROMACS) and
  Gate G0 (`scripts/gate_mvp0.py`) with fail-closed pytest/JUnit accounting.
