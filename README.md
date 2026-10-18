# Crucible

**English** | [中文](README.zh-CN.md) | [日本語](README.ja.md)

An engine-independent force-field compiler for small molecules. A parameterized molecule is
read into a single intermediate representation (`ForceFieldIR`), optionally re-parameterized
from a quantum-chemical Hessian, and written natively to **GROMACS**, **AMBER** and **OpenMM**.

**Version 0.2.0** — MVP-0 (IR + three exporters), MVP-1 (Hessian-derived bond/angle parameters),
MVP-2 (MLP torsions, with OpenMM as the reference engine per decision 0017) and MVP-3 (charge
providers: loaded, AM1-BCC, NAGL) are complete; Gates G0–G3 were re-run on the current source on
2026-10-18 and all pass. Since 2026-10-13 the code base also has an ORCA `.hess` reader, an
IR + coordinates pipeline shared by all commands (`--ir-in` / `--coords-in`), `global-fit` and
`check-staleness` commands, a Hessian fallback framework (MVP-5, Python API only), a learned-BCC
skeleton (MVP-4, no trained model), and an optional sequential torsion refinement
(`derive-torsions --refinement-rounds {0,1,2,3}`, default 0 = single pass). The 2026-10-18
closure acceptance against the design document (34 points; the current deviation list) is in
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md); see
[CHANGELOG](CHANGELOG.md).

## Features

| Area | What it does |
|---|---|
| IR | `ForceFieldIR`: atoms, bonds, angles, proper/improper torsions, LJ + charges, exclusions and 1-4 pairs; validation and a NumPy reference energy |
| Import | AMBER prmtop/inpcrd + SDF (formal charges), native reader, no ParmEd |
| Export | GROMACS top/itp/gro, AMBER prmtop/inpcrd, OpenMM System/XML — all written natively |
| Hessian | Gaussian fchk, ORCA `.hess`, xtb 6.7.1 (GFN-FF, GFN2) with convergence, charge and quality checks; GFN-FF → GFN2 → DFT-file fallback framework (Python API only) |
| Bonded derivation | mSeminario, m2Seminario, DRIH, global Hessian fit; graph-symmetry averaging |
| Torsions (MVP-2) | Rotatable-bond perception, relaxed scans with MACE-POLAR-1-M (via ASE), MM subtraction, fit to AMBER periodic terms; GAFF2 fallback with recorded reasons; optional sequential refinement rounds (`ffcompiler derive-torsions`) |
| Charges (MVP-3) | Loaded charges, AM1-BCC (antechamber/sqm), NAGL (openff-gnn-am1bcc-1.0.0); exact total charge; MLP torsions refit from cached scans when charges change (`ffcompiler assign-charges`). RESP not yet |
| Pipeline | IR + coordinates sidecar passed between commands (`--ir-in` / `--coords-in`); recorded fit dependencies with a staleness check (`ffcompiler check-staleness`) |
| Validation | Cross-engine energy/force comparison (NumPy, OpenMM ×4, sander ×2, GROMACS) and four acceptance Gates (G0–G3) |

## Requirements

Use the existing scientific environment; nothing is installed by this project.

- Python ≥ 3.12 with `numpy`, `scipy`, `openmm`, `rdkit`; for the MLP torsion and NAGL charge
  features also ASE + MACE and OpenFF NAGL (model `openff-gnn-am1bcc-1.0.0`)
- External programs (only for the corresponding features/tests):
  - `xtb` 6.7.1 — default `/home/user/software/xtb-dist/bin/xtb`, override with `FFC_XTB`
  - GROMACS — found via `FFC_GMX`, `PATH`, or `/home/user/gmx26.3/bin/gmx`; mixed precision is fine
  - AmberTools `sander` / `tleap` — cross-engine checks and fixture generation; AM1-BCC charges
    additionally use `antechamber` (`FFC_ANTECHAMBER` or `PATH`)
  - MACE-POLAR-1-M model for torsion scans — default `/home/user/MLP/mace/MACE-POLAR-1-M.model`, override with `FFC_MACE_MODEL`
  - Sobtop 2026.1.16 — only to regenerate reference data (`scripts/make_sobtop_refs.py`)

```bash
source scripts/env.sh
```

## Quick start

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

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess}` (`fchk` and `orca-hess` need `--hessian-file`)
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize` disables averaging over symmetry-equivalent terms
- `--linear-angles {error,library}`: default `error`; `library` keeps the original parameters of angles > 175°

Fit rotatable-bond torsions from MACE-POLAR relaxed scans (scans are cached under `<outdir>/scans/`):

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg` (default 10°), `--device` (default cpu), `--scan-cache` (default `<outdir>/scans`)
- `--refinement-rounds {0,1,2,3}`: extra sequential sweeps after the single pass, each bond refit
  against the current fits of the others; default 0 keeps the single pass. Charge changes
  (`assign-charges`) replay the recorded protocol from the cached scans.

Check that the exported files reproduce the IR energy in every engine:

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

Chain steps on one IR: each step writes `ir.json` and `coords.json` (the geometry it used or
returned, e.g. the Hessian geometry); the next step reads both with `--ir-in` / `--coords-in`
(atom numbers, elements and names must match; nothing is re-ordered):

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
# exit 1 if a fit's recorded inputs (charges, torsions, ...) changed afterwards
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` exposes the QMDFF-style global Hessian fit (same Hessian options; `--prior`,
`--alpha`); if any class hits the k = 0 bound it writes `ir_failed_lower_bound.json`, exports
nothing and exits 1.


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
- **Harmonic terms** are stored as E = ½ k Δq², i.e. k is the Hessian diagonal. GROMACS and
  OpenMM use k directly; AMBER uses K = k/2 (E = K Δq²), plus kcal and Å conversion.
- **DRIH** works internally in atomic units (Bohr, E_h/a₀²); stretch diagonals convert with
  E_h/a₀² ≈ 9.37583 × 10⁵ kJ mol⁻¹ nm⁻², angle diagonals with E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻²
  (no extra a₀² factor).
- Supported Hamiltonian: gas phase, no cutoff, Lorentz–Berthelot mixing, uniform AMBER 1-4
  scaling. Anything else is rejected at import.

## Validation status

All four Gates were re-run on the current source on 2026-10-18 (after the optional
sequential-refinement integration) and pass, with the same numbers as the 2026-10-13 runs;
pytest 524 + 59 passed in each. The Gates judge the default single pass; the sequential mode is
covered by unit tests and bit-for-bit replay checks (`experiments/20261018_refinement_verify/`).

| Gate | Run | Result |
|---|---|---|
| G0 — export equivalence | `jobs/gate_mvp0_20261018_2/` | passed: 11 molecules × 4 geometries, 8 comparison legs, 2596 rows, 0 failures |
| G1 — Hessian methods | `jobs/gate_mvp1_20261018_2/` | passed: all three methods match Sobtop on 12/12 reference cases; global-fit self-consistency 6.0 × 10⁻⁹ |
| G2 — MLP torsions | `jobs/gate_mvp2_20261018_2/` | passed: 12/12 rotatable bonds fitted from MLP, 0 fallbacks; all non-sander cross-engine legs pass. Known limitation (decision 0017): sander snaps a fitted phase near π (nma dihedral 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol), reported only |
| G3 — charge providers | `jobs/gate_mvp3_20261018_2/` | passed: AM1-BCC reproduces all 12 fixture charge sets (max 4.3 × 10⁻¹¹ e); total charge conserved for 3 providers × 11 molecules (≤ 8.3 × 10⁻¹⁷ e), all cross-engine comparisons pass; MLP torsions refit after a charge change: bit-identical when charges are unchanged, 12/12 bonds stay MLP with NAGL charges |

Frequency RMSE of the derived force fields vs. the reference Hessian (median over 11 molecules,
cm⁻¹; reported, not gated):

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | Global fit |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Run the Gates yourself (each writes a new `jobs/gate_*_<date>[_N]/` directory and keep old ones):

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 parallel xtb runs × 4 threads
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume OLD_G2_RUN_DIR]   # MACE scans (reused from OLD_G2_RUN_DIR/scans if given)
python scripts/gate_mvp3.py --scan-cache jobs/<this G2 run>/scans          # required: the scans of the G2 run being accepted
PYTHONPATH=src python -m pytest -m "not external" -q   # unit tests (no external programs)
```

## Known limitations

- **Known defect (MVP-2)**: on methyl benzoate bond 1-2 the MLP-fitted torsion reproduces the MACE
  profile worse than GAFF2 (profile RMSE 3.36 vs. 2.70 kJ/mol; fit RMSE 0.60). It is the only one of
  the 12 rotatable bonds where this happens. Diagnosed 2026-10-17 (`jobs/torsion_diagnosis_20261017/report.md`):
  bonds are fitted one at a time against GAFF2 for the others, so the 1-2 term absorbed the stiff GAFF2
  methyl (0-1) torsion; once 0-1 is replaced the compensation is wrong. `--refinement-rounds 3`
  brings it to 0.44 kJ/mol and the pooled residual from 0.83 to 0.22 (in-sample;
  `jobs/torsion_b3_sequential_20261017/report.md`); every update order converges to the same fixed
  point, with order-induced differences ≤ 1e-4 kJ/mol after round 3
  (`jobs/torsion_b3_order_20261018/report.md`). It stays opt-in: the default is still the single
  pass, hold-out validation (B4) is pending, and the sequential mode is not a Gate criterion.
- **Scan quality (MVP-2)**: the nma 1-3 (amide C–N) forward/backward scan hysteresis reaches
  57 kJ/mol (14/36 frames above 1 kJ/mol, concentrated in the ±60–120° barrier region). This is a
  scan defect rather than a fit defect; whether to rescan is undecided and no automatic fallback
  threshold is set.
- A bond's torsion energy is placed on a single dihedral quartet (MVP-2 plan §7 item 3); its range
  of validity is still open.
- The ORCA `.hess` reader has been verified against one closed-shell ORCA 6.1.1 sample
  (frequencies match to 5.7 × 10⁻⁸ relative).
- Linear angles (> 175°) are not derived; acetonitrile is a diagnostic case only.
- m2Seminario can give negative ring-angle force constants on GFN-FF Hessians (cyclopropane);
  Sobtop gives the same values. `derive_bonded` raises `NonPositiveForceConstantError` instead
  of silently falling back — use mSeminario or DRIH there.
- mSeminario bond constants follow the written bond direction (as Sobtop does); symmetrization
  (on by default) removes the resulting asymmetry.
- Export writes a single residue; periodic systems and other combination rules are out of scope.

## Repository layout

```text
src/ffcompiler/
  ir/          ForceFieldIR schema, validation, NumPy energies, JSON I/O, coordinates sidecar
  importers/   AMBER prmtop importer
  exporters/   GROMACS, AMBER, OpenMM writers
  formats/     native prmtop / rst7 reader and writer
  hessian/     CartesianHessian, fchk / ORCA .hess / xtb providers, fallback ladder
  bonded/      internal coordinates, mSeminario, m2Seminario, DRIH, symmetry, global fit, staleness
  torsion/     rotatable bonds, MACE-POLAR relaxed scans, MM subtraction, Fourier fit, MLP provider
  charge/      charge providers (loaded, AM1-BCC, NAGL), learned-BCC skeleton, write-back and torsion refit
  validation/  cross-engine comparison, GROMACS adapter, vibrations, pytest accounting
scripts/       Gates, fixture generation, Sobtop reference data
tests/         unit and cross-engine tests, fixtures
docs/          plans, decision records, theory
jobs/          formal runs (only manifest/acceptance/report are tracked)
```

## Documentation

- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md) —
  closure acceptance of the current version (34 design points; known defects, pending decisions
  and missing validation in §2)
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md) — derivations, units, conventions
- [docs/decisions/](docs/decisions/README.md) — decision records 0001–0021
- [docs/superpowers/plans/](docs/superpowers/plans/) — roadmap and MVP plans
- [HANDOVER.md](HANDOVER.md) — current verified state and next steps
- [CHANGELOG.md](CHANGELOG.md)

## Roadmap

MVP-2 MLP torsions — done, Gate G2 passed. MVP-3 charge providers — done, Gate G3 passed
(RESP deferred, decision 0018). MVP-4 learned BCC — provider skeleton and research proposal
only, training not approved. MVP-5 — fallback framework in place, scheduling research not
started. Next, per the closure report §2: B4 hold-out validation of the torsion refinement, the
formal external pipeline acceptance (C5) with the real MACE reference, and the MVP-4 proposal;
later, the molecule/typing entry point and presets/single-command CLI (design §1/§17/§22).
Open decisions — default bonded method, default refinement rounds, version number — are listed
there. The repository is hosted privately as CRUCIBLE.
