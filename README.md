# Crucible

**English** | [中文](README.zh-CN.md) | [日本語](README.ja.md)

An engine-independent force-field compiler for small molecules. A parameterized molecule is
read into a single intermediate representation (`ForceFieldIR`), optionally re-parameterized
from a quantum-chemical Hessian, and written natively to **GROMACS**, **AMBER** and **OpenMM**.

**Version 0.2.0** — MVP-0 (IR + three exporters) and MVP-1 (Hessian-derived bond/angle
parameters) are complete; Gate G0 and Gate G1 both pass. See [CHANGELOG](CHANGELOG.md).

## Features

| Area | What it does |
|---|---|
| IR | `ForceFieldIR`: atoms, bonds, angles, proper/improper torsions, LJ + charges, exclusions and 1-4 pairs; validation and a NumPy reference energy |
| Import | AMBER prmtop/inpcrd + SDF (formal charges), native reader, no ParmEd |
| Export | GROMACS top/itp/gro, AMBER prmtop/inpcrd, OpenMM System/XML — all written natively |
| Hessian | Gaussian fchk, xtb 6.7.1 (GFN-FF, GFN2) with convergence, charge and quality checks |
| Bonded derivation | mSeminario, m2Seminario, DRIH, global Hessian fit; graph-symmetry averaging |
| Validation | Cross-engine energy/force comparison (NumPy, OpenMM ×4, sander ×2, GROMACS) and two acceptance Gates |

## Requirements

Use the existing scientific environment; nothing is installed by this project.

- Python ≥ 3.12 with `numpy`, `scipy`, `openmm`, `rdkit`
- External programs (only for the corresponding features/tests):
  - `xtb` 6.7.1 — default `/home/user/software/xtb-dist/bin/xtb`, override with `FFC_XTB`
  - GROMACS — found via `FFC_GMX`, `PATH`, or `/home/user/gmx26.3/bin/gmx`; mixed precision is fine
  - AmberTools `sander` / `tleap` — cross-engine checks and fixture generation
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

- `--hessian-provider {gfnff,gfn2,fchk}` (`fchk` needs `--hessian-file`)
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize` disables averaging over symmetry-equivalent terms
- `--linear-angles {error,library}`: default `error`; `library` keeps the original parameters of angles > 175°

Check that the exported files reproduce the IR energy in every engine:

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

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

| Gate | Run | Result |
|---|---|---|
| G0 — export equivalence | `jobs/gate_mvp0_20261005_3/` | passed: 11 molecules × 4 geometries, 8 comparison legs, 2596 rows, 0 failures |
| G1 — Hessian methods | `jobs/gate_mvp1_20261005_4/` | passed: all three methods match Sobtop on 12/12 reference cases; global-fit self-consistency 6.0 × 10⁻⁹ |

Frequency RMSE of the derived force fields vs. the reference Hessian (median over 11 molecules,
cm⁻¹; reported, not gated):

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | Global fit |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Run the Gates yourself (both write a new `jobs/gate_*_<date>[_N]/` directory and keep old ones):

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 parallel xtb runs × 4 threads
PYTHONPATH=src python -m pytest -m "not external" -q   # unit tests (no external programs)
```

## Known limitations

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
  ir/          ForceFieldIR schema, validation, NumPy energies, JSON I/O
  importers/   AMBER prmtop importer
  exporters/   GROMACS, AMBER, OpenMM writers
  formats/     native prmtop / rst7 reader and writer
  hessian/     CartesianHessian, fchk and xtb providers
  bonded/      internal coordinates, mSeminario, m2Seminario, DRIH, symmetry, global fit
  validation/  cross-engine comparison, GROMACS adapter, vibrations, pytest accounting
scripts/       Gates, fixture generation, Sobtop reference data
tests/         unit and cross-engine tests, fixtures
docs/          plans, decision records, theory
jobs/          formal runs (only manifest/acceptance/report are tracked)
```

## Documentation

- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md) — derivations, units, conventions
- [docs/decisions/](docs/decisions/README.md) — decision records 0001–0015
- [docs/superpowers/plans/](docs/superpowers/plans/) — roadmap and MVP plans
- [HANDOVER.md](HANDOVER.md) — current verified state and next steps
- [CHANGELOG.md](CHANGELOG.md)

## Roadmap

MVP-2 MLP torsions · MVP-3 charge providers · MVP-4 learned BCC · MVP-5 fallback scheduling —
not started. The repository is hosted privately as CRUCIBLE.
