# Crucible

**English** | [中文](README.zh-CN.md) | [日本語](README.ja.md)

Crucible is an engine-independent force-field compiler for small molecules. One
parameterization lives in a single intermediate representation, `ForceFieldIR`, and is
written natively to **GROMACS**, **AMBER** and **OpenMM** — with an automated energy/force
equivalence check across all three engines on every export.

```text
AMBER prmtop/inpcrd + SDF
        │  import (native reader, no ParmEd)
        ▼
   ForceFieldIR ◄──── coordinates sidecar (every command: --ir-in / --coords-in)
        │
        ├─ bond/angle ← mSeminario · m2Seminario · DRIH · global Hessian fit
        │                Hessians: Gaussian fchk · ORCA .hess · CP2K .out · xtb (GFN-FF / GFN2)
        ├─ torsions   ← MACE-POLAR-1-M relaxed scans → Fourier terms (GAFF2 fallback, recorded)
        └─ charges    ← loaded · AM1-BCC (antechamber) · NAGL
        │
        ▼
   ForceFieldIR (complete — every parameter carries provenance)
        │
   ┌────┼─────────┐
   ▼    ▼         ▼
GROMACS  AMBER   OpenMM   → cross-engine energy / force comparison
```

**Version 0.2.0.** MVP-0 through MVP-3 are complete and gated (G0–G3 re-run on the current
source on 2026-10-18, all pass). Since 2026-10-13 the code base has additionally gained an
ORCA `.hess` reader, a CP2K `.out` reader, the IR + coordinates pipeline shared by all
commands, `global-fit` and `check-staleness` commands, periodic-box support, a Hessian
fallback framework (MVP-5, Python API only), a learned-BCC skeleton (MVP-4, no trained
model), and an optional sequential torsion refinement (`--refinement-rounds {0,1,2,3}`,
default 0 = single pass). The 34-point acceptance against the design document is
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md); the
authoritative fine-grained project state is [HANDOVER.md](HANDOVER.md); changes are listed
in [CHANGELOG](CHANGELOG.md).

## Capabilities

### One IR, three engines

`ForceFieldIR` holds atoms, bonds, angles, proper/improper torsions, LJ + charges,
exclusions and 1-4 pairs, and an optional periodic box; it validates itself and ships with
a NumPy reference energy. The three exporters write GROMACS top/itp/gro, AMBER
prmtop/inpcrd and OpenMM System/XML natively — no ParmEd, no intermediate formats. The
same IR, written to all three engines, must reproduce the same per-term energies and
atomic forces; `check-equivalence` runs that comparison (8 legs: NumPy, OpenMM ×4, sander
×2, GROMACS) and the G0 gate runs it over 11 molecules × 4 geometries.

### Bonded parameters from a Hessian

Four methods derive bond/angle parameters from a Cartesian Hessian, sharing one provider
abstraction and one unit system (atomic units internally):

- **mSeminario** — local projection, the general-purpose default candidate;
- **m2Seminario** — symmetry- and direction-consistent variant for rigid/high-symmetry systems;
- **DRIH** — redundant internal-coordinate diagonal, for rigid and special-chemistry fragments;
- **global Hessian fit** (QMDFF-inspired) — all force constants fitted against the full
  Hessian with the fixed part (LJ + Coulomb + 1-4 + torsions) subtracted, a prior
  regularization, and reported condition number, residuals and residual force at x₀;
- **localhess** (experimental) — a local-compliance method, validated for internal
  consistency only (decision 0023; not a Gate criterion).

Graph-symmetry averaging over equivalent terms is on by default (`--no-symmetrize` to
disable). Non-positive force constants raise instead of silently falling back.

### Torsions from an MLP

Rotatable bonds are perceived from the molecular graph, scanned with constrained relaxed
scans against MACE-POLAR-1-M (via ASE), the MM contribution of the current field is
subtracted, and the profile is fitted to AMBER periodic terms. Fallbacks to GAFF2 are
recorded per bond with reasons and are never silently upgraded on refits (decision 0020).
An optional bounded sequential refinement (`--refinement-rounds {0,1,2,3}`) refits each
bond against the current fits of the others; charge changes replay the recorded protocol
from cached scans bit-for-bit.

### Charges as providers

Loaded charges, AM1-BCC (antechamber/sqm) and NAGL (`openff-gnn-am1bcc-1.0.0`) share one
provider interface with exact total charge. Changing providers refits MLP torsions from
cached scans without re-scanning. A learned-BCC research skeleton (bond-antisymmetric,
charge-conserving corrections) is in place; no trained model exists yet. RESP is deferred
(decision 0018).

### Perception and force-field libraries

A geometry-based perception layer (Cordero 2008 radii ×1.15 with data-file distance
windows, angle/dihedral generation with a linear cutoff, triclinic PBC nearest image,
data-driven improper centers) and readers/matchers for Sobtop-format `LJ_param.dat` /
`bonded_param.dat` / `bondcrit.dat` parameter files (X-wildcard precedence, IDIVF,
negative periodicity, LJ conversion) support the reference-data workflows; see decision
0022 and the reverse-engineering report in `docs/theory/`.

### Provenance and staleness

Every parameter carries `source` + metadata; fits record their inputs (which charges,
which torsion versions, which Hessian identity) as per-term snapshots.
`check-staleness` fails with exit 1 if anything a fit recorded has changed since. The
IR + coordinates sidecar passes geometry faithfully between commands — atom order,
elements and names are never re-ordered.

## Getting started

Use the existing scientific environment; nothing is installed by this project.

- Python ≥ 3.12 with `numpy`, `scipy`, `openmm`, `rdkit`; for the MLP torsion and NAGL
  charge features also ASE + MACE and OpenFF NAGL
- External programs (only for the corresponding features/tests):
  - `xtb` 6.7.1 — default `~/software/xtb-dist/bin/xtb`, override with `FFC_XTB`
  - GROMACS — found via `FFC_GMX`, `PATH`, or `~/gmx26.3/bin/gmx`; mixed precision is fine
  - AmberTools `sander` / `tleap` — cross-engine checks and fixture generation; AM1-BCC
    charges additionally use `antechamber` (`FFC_ANTECHAMBER` or `PATH`)
  - MACE-POLAR-1-M model for torsion scans — default
    `~/MLP/mace/MACE-POLAR-1-M.model`, override with `FFC_MACE_MODEL`
  - Sobtop 2026.1.16 — only to regenerate reference data (`scripts/make_sobtop_refs.py`)

```bash
source scripts/env.sh
```

Convert an AMBER-parameterized molecule to all three engines:

```bash
PYTHONPATH=src python -m ffcompiler.cli convert \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output gromacs amber openmm --outdir jobs/ethanol_conversion
```

Re-derive bonds and angles from a GFN2 Hessian (xtb runs in `<outdir>/xtb/`):

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --hessian-provider gfn2 --method mseminario \
  --output gromacs amber openmm --outdir jobs/ethanol_mseminario
```

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess,cp2k-out}` (`fchk`, `orca-hess` and
  `cp2k-out` need `--hessian-file`)
- `--method {mseminario,m2seminario,drih}`; `--method localhess` selects the experimental
  method above
- `--no-symmetrize` disables averaging over symmetry-equivalent terms
- `--linear-angles {error,library}`: default `error`; `library` keeps the original
  parameters of angles > 175°

Fit rotatable-bond torsions from MACE-POLAR relaxed scans (cached under `<outdir>/scans/`):

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg` (default 10°), `--device` (default cpu), `--scan-cache` (default
  `<outdir>/scans`)
- `--refinement-rounds {0,1,2,3}`: extra sequential sweeps after the single pass;
  default 0 keeps the single pass

Check that the exported files reproduce the IR energy in every engine:

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

Chain steps on one IR — each step writes `ir.json` + `coords.json`, the next step reads
both:

```bash
# 1) charges
PYTHONPATH=src python -m ffcompiler.cli assign-charges --prmtop P --inpcrd C --sdf S \
  --method am1bcc --output openmm --outdir jobs/mol/q
# 2) bonds/angles from a Hessian, continuing from step 1
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --ir-in jobs/mol/q/ir.json --coords-in jobs/mol/q/coords.json \
  --hessian-provider gfn2 --method mseminario --output openmm --outdir jobs/mol/b
# 3) export the final IR + coordinates
PYTHONPATH=src python -m ffcompiler.cli convert \
  --ir-in jobs/mol/b/ir.json --coords-in jobs/mol/b/coords.json \
  --output gromacs amber openmm --outdir jobs/mol/final
# exit 1 if a fit's recorded inputs changed afterwards
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` exposes the QMDFF-style global Hessian fit (same Hessian options; `--prior`,
`--alpha`); if any class hits the k = 0 bound it writes `ir_failed_lower_bound.json`,
exports nothing and exits 1.

From Python:

```python
from ffcompiler.importers import from_amber
from ffcompiler.hessian.xtb import run_xtb_hessian
from ffcompiler.bonded import derive_bonded

ir, pos = from_amber("mol.prmtop", "mol.inpcrd", "mol.sdf")
h = run_xtb_hessian([a.atomic_number for a in ir.atoms], pos,
                    method="gfn2", charge=ir.total_charge, workdir="jobs/xtb")
new_ir, geom, report = derive_bonded(ir, h, method="drih", input_positions_nm=pos)
```

## Conventions

- **Units in the IR**: nm, rad, kJ mol⁻¹, elementary charge.
- **Harmonic terms** are stored as E = ½ k Δq², i.e. k is the Hessian diagonal. GROMACS
  and OpenMM use k directly; AMBER uses K = k/2 (E = K Δq²), plus kcal and Å conversion.
- **DRIH** works internally in atomic units (Bohr, E_h/a₀²); stretch diagonals convert
  with E_h/a₀² ≈ 9.37583 × 10⁵ kJ mol⁻¹ nm⁻², angle diagonals with
  E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻² (no extra a₀² factor).
- Supported Hamiltonian: gas phase, no cutoff, Lorentz–Berthelot mixing, uniform AMBER 1-4
  scaling. Anything else is rejected at import.

## Verification status

All four Gates were re-run on the current source on 2026-10-18 (after the optional
sequential-refinement integration) and pass, with the same numbers as the 2026-10-13
runs; pytest 524 + 59 passed in each. The Gates judge the default single pass; the
sequential mode is covered by unit tests and bit-for-bit replay checks
(`experiments/20261018_refinement_verify/`).

| Gate | Run | Result |
|---|---|---|
| G0 — export equivalence | `jobs/gate_mvp0_20261018_2/` | passed: 11 molecules × 4 geometries, 8 comparison legs, 2596 rows, 0 failures |
| G1 — Hessian methods | `jobs/gate_mvp1_20261018_2/` | passed: all three methods match Sobtop on 12/12 reference cases; global-fit self-consistency 6.0 × 10⁻⁹ |
| G2 — MLP torsions | `jobs/gate_mvp2_20261018_2/` | passed: 12/12 rotatable bonds fitted from MLP, 0 fallbacks; all non-sander cross-engine legs pass. Known limitation (decision 0017): sander snaps a fitted phase near π (nma dihedral 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol), reported only |
| G3 — charge providers | `jobs/gate_mvp3_20261018_2/` | passed: AM1-BCC reproduces all 12 fixture charge sets (max 4.3 × 10⁻¹¹ e); total charge conserved for 3 providers × 11 molecules (≤ 8.3 × 10⁻¹⁷ e), all cross-engine comparisons pass; MLP torsions refit after a charge change: bit-identical when charges are unchanged, 12/12 bonds stay MLP with NAGL charges |

Frequency RMSE of the derived force fields vs. the reference Hessian (median over 11
molecules, cm⁻¹; reported, not gated):

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | Global fit |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Run the Gates yourself (each writes a new `jobs/gate_*_<date>[_N]/` directory; keep old
ones):

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 parallel xtb runs × 4 threads
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume OLD_G2_RUN_DIR]   # MACE scans
python scripts/gate_mvp3.py --scan-cache jobs/<this G2 run>/scans          # required
PYTHONPATH=src python -m pytest -m "not external" -q   # unit tests (no external programs)
```

## Known limitations

- **Known defect (MVP-2)**: on methyl benzoate bond 1-2 the MLP-fitted torsion reproduces
  the MACE profile worse than GAFF2 (profile RMSE 3.36 vs. 2.70 kJ/mol; fit RMSE 0.60).
  It is the only one of the 12 rotatable bonds where this happens. Diagnosed 2026-10-17
  (`jobs/torsion_diagnosis_20261017/report.md`): bonds are fitted one at a time against
  GAFF2 for the others, so the 1-2 term absorbed the stiff GAFF2 methyl (0-1) torsion;
  once 0-1 is replaced the compensation is wrong. `--refinement-rounds 3` brings it to
  0.44 kJ/mol and the pooled residual from 0.83 to 0.22 (in-sample;
  `jobs/torsion_b3_sequential_20261017/report.md`); every update order converges to the
  same fixed point, with order-induced differences ≤ 1e-4 kJ/mol after round 3
  (`jobs/torsion_b3_order_20261018/report.md`). It stays opt-in: the default is still the
  single pass, hold-out validation (B4) is pending, and the sequential mode is not a Gate
  criterion.
- **Scan quality (MVP-2)**: the nma 1-3 (amide C–N) forward/backward scan hysteresis
  reaches 57 kJ/mol (14/36 frames above 1 kJ/mol, concentrated in the ±60–120° barrier
  region). This is a scan defect rather than a fit defect; whether to rescan is undecided
  and no automatic fallback threshold is set.
- A bond's torsion energy is placed on a single dihedral quartet (MVP-2 plan §7 item 3);
  its range of validity is still open.
- The ORCA `.hess` reader has been verified against one closed-shell ORCA 6.1.1 sample
  (frequencies match to 5.7 × 10⁻⁸ relative).
- Linear angles (> 175°) are not derived; acetonitrile is a diagnostic case only.
- m2Seminario can give negative ring-angle force constants on GFN-FF Hessians
  (cyclopropane); Sobtop gives the same values. `derive_bonded` raises
  `NonPositiveForceConstantError` instead of silently falling back — use mSeminario or
  DRIH there.
- mSeminario bond constants follow the written bond direction (as Sobtop does);
  symmetrization (on by default) removes the resulting asymmetry.
- Export writes a single residue. The IR carries an optional periodic box (`box_nm`,
  decision 0023): NumPy energies use minimum image + bond-graph unwrap, GROMACS/OpenMM
  export full triclinic boxes, AMBER prmtop only orthorhombic; periodic cross-engine legs
  are not yet part of check-equivalence.

## Repository layout

```text
src/ffcompiler/
  ir/          ForceFieldIR schema, validation, NumPy energies, JSON I/O, coordinates sidecar
  perception/  geometry bonding/angle/dihedral/improper perception, Cordero radii, PBC mirror
  fflib/       Sobtop-format LJ/bonded/bondcrit library readers and wildcard matching
  importers/   AMBER prmtop importer
  exporters/   GROMACS, AMBER, OpenMM writers
  formats/     native prmtop / rst7 reader and writer
  hessian/     CartesianHessian, fchk / ORCA .hess / CP2K .out / xtb providers, fallback ladder
  bonded/      internal coordinates, mSeminario, m2Seminario, DRIH, localhess, symmetry, global fit, staleness
  torsion/     rotatable bonds, MACE-POLAR relaxed scans, MM subtraction, Fourier fit, MLP provider
  charge/      charge providers (loaded, AM1-BCC, NAGL), learned-BCC skeleton, write-back and torsion refit
  validation/  cross-engine comparison, GROMACS adapter, vibrations, pytest accounting
scripts/       Gates, fixture generation, Sobtop reference data
tests/         unit and cross-engine tests, fixtures
docs/          plans, decision records, theory
jobs/          formal runs (only manifest/acceptance/report are tracked)
```

## Documentation

- [forcefield_compiler_design.md](forcefield_compiler_design.md) — the design document
  (acceptance baseline; next-phase sections §3.1, §6.4.4–6.4.6, §10.4–10.5 marked
  "designed, not implemented")
- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md) —
  closure acceptance of the current version (34 design points; known defects, pending
  decisions and missing validation in §2)
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md) — derivations, units,
  conventions
- [docs/decisions/](docs/decisions/README.md) — decision records 0001–0024
- [docs/superpowers/plans/](docs/superpowers/plans/) — roadmap and MVP plans
- [HANDOVER.md](HANDOVER.md) — current verified state and next steps
- [CHANGELOG.md](CHANGELOG.md)

## Roadmap

**Done and gated**: MVP-0 (IR + exporters), MVP-1/1b (Hessian methods + global fit),
MVP-2 (MLP torsions), MVP-3 (charge providers) — G0–G3 pass on the current source.

**Framework only, not researched/completed**: MVP-4 learned BCC (skeleton + research
proposal, training not approved), MVP-5 fallback scheduling (ladder in place, no
uncertainty-driven dispatch).

**Designed, not implemented** (design amendment 2026-10-25, decision 0024 — sections
§3.1, §6.4.4–6.4.6, §10.4–10.5 of the design document):

- **MVP-6 molecule entry** — SDF in, skeleton IR out: antechamber as the external GAFF2
  typer (typing-only), GAFF2 library values as prior/fallback, no silent conformer
  generation; unlocks the single-command CLI and presets (§17/§22).
- **MVP-7 scan quality** — multi-start relaxed scans with per-grid-point branch
  diagnosis; branches are diagnosis-only and never fitted; the nma 1-3 rescan is the
  first application.
- **MVP-8 equilibrium-geometry joint fit** — q⁰ as a fit variable plus gradient matching
  (the interface reserved by §6.4.3-1), opt-in, bit-exact degeneracy to the current
  global fit.
- **MVP-9 multi-conformer joint fit + alternation** — same objective with C ≥ 1 and
  shared q⁰, hold-out conformer validation, and the charge → torsions → global-fit
  alternation with ε-tolerant completion criteria.
- TorsionDrive is reserved for opt-in 2-D coupled scans; installing it into the
  production environment is an explicit pending decision.

**Pending validation**: B4 hold-out validation of the torsion refinement, the formal
external pipeline acceptance (C5) with the real MACE reference, and the Level 3–4
benchmarks (design §20).

**Open decisions**: default bonded method, default refinement rounds, version number,
TorsionDrive installation — tracked in HANDOVER.md and the closure report §2.
