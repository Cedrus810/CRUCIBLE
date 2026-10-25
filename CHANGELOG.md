# Changelog

All notable changes to ffcompiler. Dates are Asia/Tokyo. Decision records live in
[`docs/decisions/`](docs/decisions/README.md).

## [Unreleased] — README overhaul (2026-10-25; docs only)

- All three READMEs (EN / zh-CN / ja) fully rewritten from the MVP-progress style into a
  "what the tool is" structure: pipeline diagram up front, capabilities by section, real
  quickstart, verification and known-limitations sections kept, roadmap rewritten with
  explicit status tiers (done & gated / framework only / designed not implemented /
  pending validation / open decisions). Implemented features only; the MVP-6–9 design
  amendment appears solely in the roadmap as "designed, not implemented".
- The authoritative fine-grained project state remains HANDOVER.md.

## [Unreleased] — Next-phase design amendment (2026-10-25; decision 0024; design only, no code)

- Design doc amended with the four next-phase directions, all as **new sections — existing
  numbering untouched** (acceptance reports cite § numbers): §3.1 molecule entry,
  §6.4.4–§6.4.6 joint fitting and alternation, §10.4–§10.5 scan protocol, §23 MVP-6–MVP-9,
  §26 boundaries, §27 references; header carries a revision note.
- §3.1 SDF → skeleton IR: antechamber as the external GAFF2 typer (typing-only, no sqm;
  charges stay with the §11 provider layer), GAFF2 library values as prior / fallback /
  `policy: library` answer, geometry contract with no silent conformer generation
  (SDF-only first version; consistent with decision 0023's mol2 rejection).
- §6.4.4 equilibrium-geometry joint fit (the interface reserved by §6.4.3-1): q⁰ becomes a
  fit variable plus a gradient-matching term; linear in k given q⁰ (curvature + gradient +
  prior rows stacked), opt-in with bit-exact degeneracy to §6.4.1. §6.4.5 multi-conformer
  joint global fit shares the same objective (C = 1 degenerates; shared q⁰; ETKDG →
  same-provider optimization → dedup → energy-window selection; hold-out conformer
  validation institutionalized). §6.4.6 closes the §6.4.2 alternation gap:
  charge → derive-torsions → global-fit rounds (≤ N = 3) with ε-tolerant completion
  criteria and an ε-aware check-staleness.
- §10.4 native multi-start relaxed scans with per-grid-point branch bookkeeping: branches
  are diagnosis-only and never fit targets; Fourier fitting consumes the per-point minimum
  envelope; default scan protocol unchanged and caches backward compatible; no automatic
  hysteresis thresholds (user ruling of 2026-10-17 stands). §10.5 reserves TorsionDrive
  for opt-in 2-D coupled scans — installing it into the production env is an explicit user
  decision, not a construction detail.
- No code, tests, or Gates touched this round.

## [Unreleased] — Sobtop absorption round 2 (2026-10-24; decisions 0022/0023; unit suite green)

- New `perception` and `fflib` packages (decision 0022): geometry bonding/angle/dihedral/improper
  perception with Cordero 2008 radii ×1.15, bondcrit.dat distance windows, triclinic PBC nearest
  image; readers/matchers for Sobtop's `LJ_param.dat` / `bonded_param.dat` / `bondcrit.dat`
  (X-wildcard precedence specific-first/forward-first, IDIVF, negative periodicity, LJ conversion,
  guess constants). 24 tests against sandbox ground truth.
- `--method localhess` (experimental, decision 0023): local-compliance bond/angle derivation
  reverse-engineered from Sobtop's hidden k_method=5 (`k = 1/(b^T pinv(H_loc) b)`); no reference
  oracle exists (Sobtop has no menu entry), so validation is internal-consistency only.
- `--hessian-provider cp2k-out` (decision 0023): CP2K `.out` VIB Hessian reader; bit-checked
  against the official FeCl2 example. The printed matrix is a forward-difference raw Jacobian —
  symmetrized on read with the raw asymmetry recorded in `meta`; Hessian quality checks still
  gate it.
- IR periodic box `box_nm` (decision 0023): JSON round-trip, validation (shape/finite/det>0),
  NumPy reference energies use bond-graph unwrap + minimum image (wrap-invariant), GROMACS gro
  box line (orthorhombic 3 / triclinic 9 values), OpenMM periodic box vectors (lower-triangular),
  AMBER IFBOX=1 + BOX_DIMENSIONS + rst7 box line (orthorhombic only; triclinic raises).
- Full reverse-engineering report in-tree: `docs/theory/sobtop_reverse_engineering.md`.
- mol2 as a structure input is explicitly rejected (user decision): the `charge/loaded` provider
  already consumes mol2 with a different semantics; structure input stays prmtop/inpcrd/sdf.

## [Unreleased] — Optional sequential refinement (2026-10-18; verified, Gates G0–G3 re-run and passed: `jobs/gate_mvp*_20261018_2/`)

- `derive-torsions --refinement-rounds {0,1,2,3}` adds bounded sweeps after the original single pass;
  default 0 keeps single-pass behavior. Charge changes replay the same recorded protocol from the
  original replaced terms and cached scans; original fallback bonds remain unchanged.
- Reports include per-round in-sample profile metrics and the canonical update order. Fixed proper
  dependencies are checked for staleness; reports cannot mutate recorded IR provenance.
- Analytic coupled-profile tests and zero/three-sweep pipeline coverage. Written without running;
  verified by the receiver on node-A: related files 66 passed, `pytest -m "not external"` 524 passed;
  `refinement_rounds=3` reproduces the B3 experiment bit for bit and the default 0 reproduces the
  previous single pass bit for bit (`experiments/20261018_refinement_verify/`). See
  [`2026-10-18-sequential-refinement-handover.md`](docs/superpowers/plans/2026-10-18-sequential-refinement-handover.md).

## [Unreleased] — 2026-10-17 build: ORCA .hess reader, C1–C5 pipeline contract, MVP-5 fallback framework, MVP-4 skeleton. Written on the login node untested; **accepted and fixed the same evening on node-A: `pytest -m "not external"` 502 passed** (`docs/superpowers/plans/2026-10-17-unbuilt-modules.md` §0)

### Fixed (2026-10-17; decision 0021)
- `RunContext` and `scripts/run_logged.sh` no longer probe Git in the NFS workspace.
  Manifest Git fields remain null; source/input identity records and command exit-code logging remain.

### Added (2026-10-18)
- `scripts/b3_order_dependence.py`: all bond orders per molecule for the 3-round sequential refit
  (`jobs/torsion_b3_order_20261018/`): every order tends to the same fixed point; the order-induced
  difference in fitted torsion profiles is ≤ 2.5e-2 kJ/mol after round 1, ≤ 1.2e-3 after round 2,
  ≤ 1e-4 after round 3 (in-sample; not checked against the joint least-squares solution).
- `tests/unit/test_pipeline_chain.py::test_local_route_chain_with_torsions`: C2 chain in the design order
  (Hessian bond/angle → charges → torsions → three backends) over `--ir-in`/`--coords-in`; all parameters
  and the Hessian geometry survive, `check-staleness` is clean, and a later charge change refits the
  torsions and refreshes the dependency record.

### Verified (2026-10-18)
- Gates G0–G3 re-run on the current source (`jobs/gate_mvp{0,1,2,3}_20261018/`), all passed with the same
  numbers as 2026-10-13; pytest 502 + 59 passed in each; `git_commit` null (decision 0021).

### Added (2026-10-17, B1/B2 diagnosis)
- `scripts/diagnose_torsion_profiles.py`: read-only diagnosis of a G2 run (no MACE/xtb). Reproduces the
  G2 metrics and decomposes each bond's profile on its own scan frames for GAFF2 / only-this-bond /
  final IRs, with three energy identities checked. Run `jobs/torsion_diagnosis_20261017/`: the methyl
  benzoate 1-2 defect comes entirely from replacing the other bonds' torsions (mainly the 0-1 methyl
  rotor) after single-pass fitting; nma 1-3 scan shows 57 kJ/mol hysteresis.
- `scripts/refit_torsions_sequential.py`: B3 mechanism experiment (read-only, provider defaults unchanged).
  3-round sequential refit, each bond on its own scan against the others' current fits (not block
  coordinate descent of the joint least squares). Run `jobs/torsion_b3_sequential_20261017/`: methyl
  benzoate 1-2 3.36 → 0.44 kJ/mol, pooled residual 0.83 → 0.22, no bond notably worse (in-sample).

### Fixed (2026-10-17, acceptance of the build below)
- CLI `check-staleness` crashed on every call (`AttributeError: coords_in`).
- `staleness` compared the whole `h_fixed_dependency` record instead of its `inputs` snapshot, so every
  global_fit IR was reported stale right after fitting.
- ORCA `.hess` reader rewritten against a real ORCA 6.1.1 file (`$orca_hessian_file` header, `$hessian`
  column-blocked full matrix with dimension 3N, `$atoms` lines `El mass x y z`); the original guessed
  layout could not read real files. Frequencies from the reader match ORCA's own to 5.7e-8 relative
  (`experiments/20261017_orca_hess_format/`, fixture `tests/fixtures/hessian/orca/water.hess`).
- `check-equivalence --ir-in` no longer attaches the sander-tleap leg (the original prmtop is not the
  Hamiltonian of a modified IR; the leg would always fail).
- Coordinates sidecar `ffcoords/0.2`: `atom_order` carries atom names so same-element permutations are
  rejected; `--ir-in` with `--prmtop` but without the full prmtop/inpcrd/sdf set now says on stderr that
  the topology cross-check was not done.
- MLP MM-subtraction snapshot includes impropers (other bonds' propers intentionally excluded — B2).
- `LearnedBCCCharges` requires an explicit `base_provider` (defaulting to the IR's current charges
  stacked corrections on repeated assignment).
- `global-fit` lower-bound failures write `ir_failed_lower_bound.json` / `coords_failed_lower_bound.json`
  instead of the standard names; the fallback DFT tier checks the file's atomic numbers.

### Added (2026-10-17)
- **ORCA .hess reader** (`ffcompiler.hessian.orca`, roadmap D3): see Fixed above (rewritten against a
  real ORCA 6.1.1 file).
- **Coordinates sidecar** (`ffcompiler.ir.coords`, 2026-10-15 plan C1): `ffcoords/0.2` JSON travelling
  with the IR — explicit atom-order pairing (rejects permutations, never pairs by filename), source
  provenance, nm units. CLI gains `--ir-in`/`--coords-in` everywhere; `derive-bonded` emits
  `ir.json`/`coords.json`/`bonded_report.json`; `convert`/`check-equivalence` accept the final IR.
- **Dependency provenance & staleness** (C3/C4): `global_fit` records `h_fixed_dependency`
  (charges/torsions/LJ/scales snapshot + Hessian identity); MLP torsion records `mm_dependency`.
  CLI `check-staleness` exits 1 on stale fits; CLI `global-fit` exposes the QMDFF-style fit with an
  explicit `failed_lower_bound` artifact state (no engine export, exit 1) per decision 0014.
- **Hessian fallback framework** (`ffcompiler.hessian.fallback`, MVP-5 design §8.1/§16/§23):
  injectable tiers, `check_hessian` quality gate, `gfnff → gfn2 → dft(file)` ladder; warnings do not
  escalate; no CLI auto mode.
- **learned BCC skeleton** (`ffcompiler.charge.learned_bcc`, MVP-4 design §13–15): charge-conserving
  bond-correction provider plumbing with locked sign convention; missing artifacts raise; research
  proposal `2026-10-17-mvp4-learned-bcc-proposal.md` leaves S1–S6 open, training forbidden until
  approval.
- 49 new/rewritten unit/combination tests (coords, staleness, fallback, learned BCC, ORCA, pipeline chain
  driven by a synthetic harmonic fchk — no xtb needed).

## [Unreleased] — MVP-3: charge providers; MVP-2: MLP torsions (sander legs report-only, decision 0017). Gates G0–G3 re-run on one source tree 2026-10-13, all passed (`jobs/gate_mvp{0,1,2,3}_20261013/`, `docs/history/ACCEPTANCE_20261013_A3.md`)

### Fixed (2026-10-13; decision 0020)
- `refit_torsions`: only bonds that were `mlp` are refitted; original fallback bonds keep their row and
  GAFF2 terms even if a matching scan later appears in the cache (`MLPTorsionProvider(fit_bonds=...)`).
- Gate G2 `--resume DIR`: exits with an argument error if `DIR/scans/` does not exist (no silent full
  rescan); warns if it holds no scans for the selected molecules.
- Gate G3 `--scan-cache` is now required and must be a directory (no default pointing at an old G2 run).

### Fixed (2026-10-11, A1/A2 close-out; decision 0019)
- Gates G1/G2/G3: runs over a subset (`--providers` / `--molecules`) record `scope: "subset"` and
  never `passed: true`; full runs record `scope: "full"`.
- Gate G2 `--resume`: copies the old run's `scans/` into the new run directory and never writes to the
  old one (`scan_cache_seeded_from` in the manifest).
- `assign_charges`: rejects provider output with the wrong shape or non-finite values (`ChargeError`)
  instead of truncating; `provenance["charge"]["previous"]` keeps the earlier charge source.
- `refit_torsions`: bonds that were fallbacks keep their original reason (e.g. `mlp_error`) instead of
  being relabelled `no_cached_scan`; the scan cache is only read.
- `derive-torsions --reference mm`: the scan-cache key and torsion provenance record the reference
  potential identity (charges at scan time, literal values); MACE cache keys are unchanged.

### Added
- **Torsion subsystem** (`ffcompiler.torsion`): rotatable-bond perception, two-pass ASE constrained
  relaxed scans (10°, FixInternals + LBFGS), MM subtraction, Fourier fit to AMBER periodic terms
  (n = 1..6, arbitrary phase), `MLPTorsionProvider` with GAFF2 fallback reasons and a resumable scan
  cache; reference potential MACE-POLAR-1-M via mace's ASE calculator (total charge passed; not via
  openmm-ml, decision 0016).
- **CLI**: `ffcompiler derive-torsions` (`--reference mace|mm`).
- **Gate G2** (`scripts/gate_mvp2.py`).
- **Charge subsystem** (`ffcompiler.charge`, decision 0018): `ChargeProvider` with three providers —
  `LoadedCharges` (keep IR charges or read mol2/text), `AM1BCCCharges` (antechamber/sqm, same
  settings as the fixtures), `NAGLCharges` (openff-gnn-am1bcc-1.0.0); residual equalized so the total
  charge is exact; topology checked by atomic numbers and bond set before computing. RESP deferred.
- **Torsion refit on charge change**: MLP torsion provenance records the replaced GAFF2 terms (with
  their original index); `assign_charges` refits MLP torsions from the scan cache (`refit_torsions`),
  never falling back silently (`TorsionRefitError`).
- **CLI**: `ffcompiler assign-charges --method loaded|am1bcc|nagl`; `derive-torsions --scan-cache`
  (default `<outdir>/scans`) and `--ir-in` (checked against `--prmtop`/`--sdf`).
- **Gate G3** (`scripts/gate_mvp3.py`).

### Changed
- `bonded.symmetry.term_classes` accepts dihedrals.
- `ForceFieldIR.replace_terms` accepts `atoms=`.
- Gate G2: criterion G2.3 is judged with OpenMM as the reference; sander legs still run but are
  report-only (fitted phases near 0/π may be snapped by sander; known limitation, decision 0017).
  Each molecule's cross-engine report is saved as `compare/<molecule>/report.json`.

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
