# Changelog

All notable changes to ffcompiler. Dates are Asia/Tokyo. Decision records live in
[`docs/decisions/`](docs/decisions/README.md).

## [Unreleased] — Batch 1 code layer, continued: xtb_sp stage and multistart report (2026-11-24)

### Added
- `hessian/xtb.py`: `run_xtb_singlepoint` (GFN2/GFN-FF single point + gradient at the given
  geometry, `--grad --chrg Q --uhf 0`, stale gradient files removed before the run, same
  closed-shell pre-check as the Hessian runner) and `read_xtb_gradient` (Turbomole gradient
  file, last cycle, Eh and Eh/Bohr). Format checked against real xtb 6.7.1 output.
  `run_xtb_hessian` is unchanged.
- `scripts/coverage_diagnose.py`: the `xtb_sp` stage is now implemented (opt-in,
  `--xtb-threads`, default 4); failures are recorded as gaps under the existing verdict rules.
- `scripts/torsion_multistart_report.py`: read-only multistart diagnostic report (plan §3.1)
  from a multistart `torsion_report.json`, optionally the scan cache and a single-start
  baseline report; writes report.md, per_point.csv, summary.json. Lists numbers only; no
  thresholds, no verdicts.

- `scripts/coverage_diagnose.py`: opt-in stages `fallback` (GFN-FF --ohess → mSeminario,
  counts linear angles kept at GAFF2 library values and notes that the fast preset would reject
  them; non-positive k recorded as a failure) and `scan_probe` (first rotor, current
  MLPTorsionProvider with MACE at 10°, timed, cached per molecule; skip without rotors);
  new `--mace-threads` and `--device`. `scan_probe` joins the verdict stage list.

### Changed
- Ions and zwitterions use xtb GFN2 as the torsion-scan reference potential (decision 0028):
  `torsion/xtb_ref.py` (`is_ion_or_zwitterion`, ASE `XtbGFN2Calculator` calling `xtb --gfn 2 --grad`
  per step); `derive-torsions --reference` gains `auto` (new default: GFN2 for ions/zwitterions,
  MACE otherwise) and `xtb-gfn2`; `run` uses `auto`; coverage refcheck/scan_probe follow the same
  rule. Neutral molecules and the G2/G3 gates are unchanged.

### Fixed (2026-11-25)
- `run` passed derive-torsions a Namespace without `multistart`/`multistart_threshold` (since
  2026-11-22), so `run` failed at the torsion step; a new test checks every `args.X` the handler
  reads.

### Changed (continued)
- `run_xtb_hessian` restarts `--ohess` from xtb's `xtbhess.xyz` when the optimisation stops at a
  saddle point (lowest frequency below the same −20 cm⁻¹ line that `check_hessian` rejects), up to
  `max_imaginary_restarts=2`; restarts are recorded in `meta["imaginary_restarts"]` and in
  `provenance.bonded.hessian_imaginary_restarts` only when they happen (decision 0027). Without a
  restart the result is bit-identical (ethanol hessian file byte-equal to the earlier run).
  Acetone, methyl acetate, styrene and glycine zwitterion now pass the fast-preset bonded step.
- Coverage verdict (user decision 2026-11-24, same reasoning as decision 0013): when every
  failing export row comes from the mixed-precision GROMACS leg, the export stage records a
  `gmx_mixed_precision` gap instead of blocking. Failures on any other leg, leg errors and force
  failures still block. Acetonitrile now reads ok_with_gaps.

- `scripts/make_coverage_set.py` and `benchmarks/coverage_set/` (molecules.json + 22 new SDFs):
  the 34-molecule coverage table approved by the user; new SDFs use the fixture recipe
  (`make_fixtures.build_sdf`, ETKDGv3 seed 20261004 + MMFF94); the 12 fixtures are referenced in
  place. `tests/unit/test_coverage_set.py` checks count, plan categories, 3D SDFs and charges.

### Fixed
- `scripts/coverage_diagnose.py` (2026-11-22 code, never run on real molecules before): the entry
  stage reported ir.json/coords.json paths without writing them; `_write_record` was called but
  not defined; downstream stages used `read_rst7` output (Å) as nm, scaling geometries 10×.
  Downstream positions now come from the entry stage's nm coordinate sidecar. Covered by the new
  external test `test_driver_light_stages_on_real_ethanol`.

### Verified
- Related tests on this machine: test_hessian_xtb + test_coverage 30 passed (non-external);
  new external xtb single-point test 1 passed; test_torsion_multistart_report 4 passed.
- Smoke `experiments/20261124_multistart_report_smoke/` (ethanol, mm reference, 20° grid,
  `--multistart 2`): report produced; the single-start baseline sharing the cache finished in
  1.9 s with all cache files untouched, confirming the default {x0} cache key matches the
  multistart x0 leg. Full regression not run.
- Coverage-driver smoke `experiments/20261124_coverage_stages_smoke/` (ethanol, acetonitrile,
  benzene; all nine stages): ethanol ok (scan_probe 250 s), benzene ok, acetonitrile blocked by
  the gmx leg on its near-linear angle (the mixed-precision effect of decision 0013).
  test_coverage 25 passed including 2 external tests.
- Full coverage diagnostic over 34 molecules (`jobs/coverage_20261124/`): 26 ok / 8 gaps / 0 blocked
  before decision 0027, 30 ok / 4 gaps / 0 blocked after it (remaining gaps: two ions pending the
  reference-potential decision, acetonitrile, oxirane). Pre-fix results kept in `prefix_0027/`.
  After decision 0028 the three ion/zwitterion entries were rerun with GFN2: 32 ok / 2 gaps /
  0 blocked (acetonitrile, oxirane); previous records and MACE scans kept in `prefix_0028/`.
- Light-stage smoke over all 34 molecules (`experiments/20261124_coverage_light_smoke/`): 0 blocked,
  31 ok, 3 ok_with_gaps (acetonitrile gmx mixed precision; acetate and methylammonium pending the
  ion reference-potential decision); no frcmod ATTN lines.

## [Unreleased] — Batch 1 code layer: coverage diagnostics, MVP-7 multistart scans, independent-validation assembly (2026-11-22)

Plan: `docs/superpowers/plans/2026-11-20-coverage-validation-batch1-plan.md`. Code and
unit tests only, per the user's instruction; no scans, no MD, no external tools, no
environment changes, no git operations.

### Added
- Coverage diagnostics (§A): `validation/coverage.py` (record schema, parmchk2 frcmod
  ATTN parsing, locked verdict rules — blocked on entry/typing failure or non-equivalent
  export, ok_with_gaps on frcmod gaps / reference-applicability / fit fallbacks — table
  loading, summary rendering) and the driver `scripts/coverage_diagnose.py` (light stages
  entry/rotors/frcmod/export/refcheck by default; mace_smoke opt-in; xtb_sp deliberately
  not implemented this round and errors when named; per-molecule record.json enables
  resume-by-existence). No molecule table built yet (34-molecule draft pending decision 1).
- MVP-7 multistart scans (§B, design §10.4; default {x0} path bit-identical — unit-checked:
  no `start` block in the default cache key, collect_candidates leaves the main arrays
  untouched, no multistart provenance block on the default path):
  `clustering.py` (reflection-free Kabsch, aligned heavy-atom RMSD, energy-ordered greedy
  clustering shared by §B branches and §C1 basins); `torsion/scan.py` grows optional
  `collect_candidates` (per-direction candidate arrays) and `grid_origin_rad` (aligns a
  non-x0 start's grid to x0's φ0); `torsion/multistart.py` (`etkdg_starts`: embed →
  reference-potential optimization → dedup, identities carry seed/params/RDKit version into
  the cache key; `merge_start_scans`: per-grid-point envelope + branch diagnostics —
  branches are diagnostics only, never enter the fit objective); `MLPTorsionProvider`
  accepts `multistart=` (fresh tuning spec, or refit form replaying archived identities
  from the cache without calling RDKit/the reference); candidate blocks persist as `cand_*`
  (old caches stay readable); `refit_torsions` restores multistart identities; CLI gains
  `--multistart N --multistart-threshold nm` on derive-torsions (both mace and mm-self).
  Branch threshold defaults to 0.05 nm (build-time ruling, overridable).
- Independent validation (§C): `validation/holdout.py` (basin-wise whole-basin-in/out
  holdout split with a holdout-fraction cap; single-basin scans yield an empty holdout
  with an explanatory note); `validation/variants.py` (`substitute_terms` /
  `comparison_matrix`: bonded / torsion / charge / vdw single swaps plus full, each
  variant passes `validate`); `validation/liquid.py` (`pack_box`: deterministic cubic-grid
  packing at a target density, errors instead of silently lowering the density;
  `to_liquid_system`: PME + LJ cutoff + dispersion correction with exporter-consistent
  exception semantics — the NoCutoff exporter path is untouched, decision 0023 semantics;
  `liquid_dynamics`: NPT Simulation on the CPU platform).

### Verified
- Per-file unit runs on this machine: test_coverage.py 16 passed; test_torsion_multistart.py
  16 passed (142 s, GAFF2 MM as the reference potential — no MACE); test_holdout_variants.py
  + test_liquid.py 17 passed.
- The full `pytest -m "not external"` regression was NOT run (stopped by user instruction,
  2026-11-22); external tests and the G0–G3 gates were not touched. The latest full-gate
  pass remains the 2026-10-28 round 3.

## [Unreleased] — MVP-6 molecule entry & single-command CLI (2026-11-11)

### Added
- `ffcompiler/molecule/io.py`: SDF input contract (design §3.1) — parse, require 3D
  coordinates (no silent conformer generation), formal-charge vs net-charge check,
  0.5 Å atom-overlap gate, sha256 for provenance.
- `ffcompiler/typing/gaff2.py`: antechamber `-c bcc -at gaff2` (the §3.1-sanctioned fallback;
  this AmberTools has no `-c zero`) → charge equalization → parmchk2 → tleap → the existing
  prmtop importer. Command lines, AmberTools version and input hashes recorded under
  `provenance["entry"]`; term sources stay `gaff2-library`.
  `equalize_mol2_charges`/`ambertools_version` moved here from scripts/make_fixtures.py
  (math unchanged line-for-line; fixtures not regenerated).
- `ffcompiler/presets.py`: design §17 presets verbatim (fast/balanced/high/rigid).
- CLI: `enter` (SDF → skeleton IR + coords sidecar) and `run` (design §22 single command:
  entry → derive-bonded → assign-charges → derive-torsions → one-shot export, orchestrating
  the existing subcommand handlers). The literal form `ffcompiler ligand.sdf --preset fast ...`
  works via an argv shim. Unavailable preset providers error explicitly: `new_charge` is
  MVP-4 (not chartered); `--charge entry` keeps the skeleton AM1-BCC charges (decision 0024
  ruling 1); design name `orca` maps to the `orca-hess` provider. The rigid preset's
  `rigid_torsion: drih` is a design-table inconsistency (DRIH is a bonded method) — not
  executed, recorded in run_manifest as pending clarification.
- A4 cleanups (2026-10-11 plan §3): malformed-mol2 errors carry file:line (charge/loaded.py);
  early step_deg point-count check on both CLI scan entries (scan.py itself untouched);
  TorsionRefitError message names `--scan-cache`; scan failure messages carry the exception
  class (torsion/scan.py); duplicate find_rotors call removed (scripts/gate_mvp3.py);
  G1 `--providers` help text (scripts/gate_mvp1.py).

### Verified
- `pytest -m "not external"`: 629 passed (was 603; 26 new unit tests).
- Design §3.1 acceptance point 2: all 12 committed fixtures rebuilt from their original SDFs
  match the prmtop-imported IRs on names/elements/GAFF2 types, LJ, exclusions/1-4 and every
  library bond/angle/torsion value (charges within 1e-6); external test 12/12 (5.3 s).
- End-to-end smoke of the single command on ethanol (`--preset fast --charge entry`): entry
  (real antechamber/sqm) → bonded (xtb) → torsion (MACE) all execute on this machine; the
  formal G4 gate run needs the gate machine (resource report pending).

## [Unreleased] — Constrained relaxation check, phase 1 (2026-11-11)

### Added
- `scripts/torsion_relaxed_check.py`: pre-fit, reference-free test of whether the re-relaxed
  torsion profile differs from the fixed-geometry one, and whether the torsion-independent
  background relaxation predicts it before fitting (definitions and criteria locked in
  `docs/superpowers/plans/2026-11-01-relaxed-check-plan.md`). Read-only on the 10-13 G2 cache,
  MM only (no MACE, no xtb, no rescan); resume per bond via existing frames CSVs.
- Custom `MultiDihedral` ASE constraint in the script: simultaneous Gauss–Newton position
  projection plus Euclidean force projection. Needed because ASE FixInternals (sequential
  mass-weighted corrections, epsilon 1e-7 deg) does not converge for dmso's four co-axial
  dihedral constraints at any step size or tolerance that meets the 1e-6 rad drift gate.
  Feasible set, fmax criterion and optimizer unchanged; applied uniformly to all 12 bonds.
- `jobs/torsion_relaxed_check_20261111/report.md` and `manifest.json` (only these job files
  are tracked): acceptance, per-bond X/Y/P0/P1, confusion tables, per-frame mechanism ratios,
  stop-condition verdict.

### Verified
- Acceptance vs prior artifacts: P0 matches the B2 best-shift metrics to 2.7e-11 kJ/mol on all
  12 bonds; methyl benzoate background relaxations match the 2026-11-01 report (FixInternals,
  different machine) to 4.6e-6 kJ/mol; all 24 rows pass the self-checks (worst constraint drift
  9.7e-13 rad).
- Locked stop conditions (2026-11-01 plan §3): none triggered. Y > sigma = 0.41 kJ/mol on 5/11
  bonds (only basis) — including ethanol 0-1/1-2 and methyl benzoate 2-4, whose fixed-geometry
  metrics are near-perfect (0.02–0.04). X classifies with 0 false positives and more true
  positives than P0/P1 in both bases; per frame it explains ≤5% of the re-relaxation change.
  nma 1-3 reported separately (bad scan). Phase 2 (new molecules, new MACE scans) is a pending
  user decision; nothing enters the pipeline or the Gates.

## [Unreleased] — External tool configuration (2026-10-28; decision 0026)

### Added
- `ffcompiler/toolconfig.py`: one lookup for every external tool — `FFC_<TOOL>` → per-user
  `~/.config/crucible/tools.toml` (or `FFC_CONFIG`; never committed) → `PATH` → legacy `~` default.
  `xtb` is pinned and never taken from `PATH`. A written-but-unusable path is an error; a missing
  tool only disables the features that need it.
- `ffcompiler init-config` (fills in what it can detect, placeholders for the rest; no overwrite
  without `--force`) and `ffcompiler check-env` (tool table with source/path/version plus a
  feature-availability table; exit 1 only on configuration errors).
- `crucible.example.toml` template with placeholder paths; `scripts/env.sh` layout overridable via
  `FFC_MAMBA_ROOT` / `FFC_ENV_NAME` / `FFC_MODULEFILES` / `FFC_GMX_MODULE` (defaults unchanged).

### Changed
- xtb, gmx (incl. extra library dirs), antechamber, sander, MACE model, `make_fixtures.py`,
  `make_sobtop_refs.py` resolve through toolconfig; `XTB_BINARY` and `mace.DEFAULT_MODEL` removed.
- CLI `check-equivalence` skips the sander legs when sander is not found (`report.skipped`);
  `compare()` keeps failing by default, so Gates still require it.
- `make_sobtop_refs.prepare_sobtop` needs the Sobtop source directory only when `tools/sobtop/`
  has not been copied yet.
- Tests: `pytest -m "not external"` 603 passed; unit tests are isolated from the developer's own
  config via an autouse fixture in `tests/conftest.py`.

### Verified
- Config-only run (second account, no HOME change, no `FFC_<TOOL>`): full `pytest -m external` 60 passed;
  Gates G0–G3 run 3 (`jobs/gate_mvp{0,1,2,3}_20261028_3/`) all passed, every report.md byte-identical to
  run 2; pytest 603 + 60, 0 failed.

## [Unreleased] — Review fixes for the 2026-10-24 code (2026-10-28; decision 0025)

### Fixed
- Periodic box is simulation-cell metadata only: `ir/energy.evaluate` no longer unwraps or
  applies minimum image (drops the `box_nm` argument). The previous NumPy semantics matched no
  engine — OpenMM export uses NoCutoff (no PBC), sander checks use `ntb=0`, GROMACS checks use an
  oversized box. Measured on ethanol in a 3 nm box with wrapped coordinates: NumPy −8.45 vs
  OpenMM 8.8e6 kJ/mol before; identical after.
- AMBER periodic prmtop (IFBOX=1) now writes `SOLVENT_POINTERS` and `ATOMS_PER_MOLECULE`;
  sander previously aborted with `Flag "SOLVENT_POINTERS" not found`. rst7 box line now carries
  6 values (a b c 90 90 90).
- GROMACS export rejects non-lower-triangular boxes (same check as OpenMM).
- CLI `cp2k-out` keeps the file's QUICKSTEP geometry so the topology check compares it against
  the input geometry (it was overwritten by the input geometry, making the check vacuous); only a
  file without a coordinate table falls back to the input geometry, with a stderr note.
- `perceive_bonds`: bondcrit windows replace the radius criterion for listed element pairs
  (was a union), per `setbondcrit` in the decompilation.
- `fflib`: negative periodicity keeps its phase (AMBER continuation flag; was phase → −phase);
  `find_improper` follows `findffimp` (fixed center, order-free outer atoms, 4th column exact).
- `test_default_model_path_is_polar_m` no longer hard-codes `/home/user/...`.
- Tests: `pytest -m "not external"` 579 passed; the touched external tests (sander periodic prmtop,
  MACE, CLI; 34) pass; see `experiments/20261028_review_pytest/`.

### Verified (2026-10-28)
- Gates G0–G3 run 2 (`jobs/gate_mvp{0,1,2,3}_20261028_2/`), after the gfnff CLI test was changed to pass
  HOME and FFC_XTB to its subprocess: **all four passed**; every report.md is byte-identical to 2026-10-18
  run 2; pytest 579 + 60 passed, 0 failed in each.
- Run 1 (`_20261028`, superseded) as a second (non-author) account: every gate-specific criterion
  matches 2026-10-18 run 2 bit for bit (G0 report.md identical, G2 metrics.csv identical). All four
  gates are marked not passed solely because `test_cli_derive_bonded_ethanol_gfnff` failed
  (59 passed, 1 failed in each external pytest): the test strips HOME from the CLI subprocess, so
  `~` resolved to the second account's home and xtb was not found. The same command exits 0 once
  `HOME` points at the author's home is set.

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
