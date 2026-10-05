"""Real GROMACS gas-phase single points, without changing shell configuration.

compare() applies the file tolerance; mixed builds use the fp32-aware policy of decision 0013.
Missing executables may be skipped; execution or parsing failures never may.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from ffcompiler.ir.energy import evaluate

TERMS = {
    "Bond": "bond", "Angle": "angle", "Proper-Dih.": "proper",
    "Per.-Imp.-Dih.": "improper", "LJ-14": "lj14", "Coulomb-14": "coulomb14",
    "LJ-(SR)": "lj", "Coulomb-(SR)": "coulomb", "Potential": "total",
}


def find_gromacs() -> str | None:
    explicit = os.environ.get("FFC_GMX")
    if explicit:
        if not Path(explicit).is_file() or not os.access(explicit, os.X_OK):
            raise ValueError(f"FFC_GMX is not executable: {explicit}")
        return str(Path(explicit).resolve())
    if found := shutil.which("gmx"):
        return found
    # User-provided installation; the ordinary AVX2 variant is the safe fallback.
    known = Path("/home/user/gmx26.3/bin/gmx")
    return str(known) if known.is_file() and os.access(known, os.X_OK) else None


def _environment(executable: str) -> dict[str, str]:
    env = dict(os.environ)
    base = Path(executable).resolve().parents[1]
    roots = [base / "lib64", base / "lib", Path("/home/user/library/HDF5/lib")]
    libraries = [str(p) for p in roots if p.is_dir()]
    if env.get("LD_LIBRARY_PATH"):
        libraries.append(env["LD_LIBRARY_PATH"])
    env["LD_LIBRARY_PATH"] = ":".join(libraries)
    env["OMP_NUM_THREADS"] = "1"
    return env


def _run(executable, args, workdir, stem, *, stdin=None, check=True):
    argv = [executable, *args]
    proc = subprocess.run(argv, cwd=workdir, env=_environment(executable),
                          input=stdin, text=True, capture_output=True, timeout=180)
    (workdir / f"{stem}.stdout").write_text(proc.stdout)
    (workdir / f"{stem}.stderr").write_text(proc.stderr)
    (workdir / f"{stem}.command.json").write_text(json.dumps(
        {"argv": argv, "stdin": stdin, "returncode": proc.returncode}, indent=2))
    if check and proc.returncode:
        raise RuntimeError(f"GROMACS {stem} exit={proc.returncode}: {(proc.stderr + proc.stdout)[-1800:]}")
    return proc


def parse_xvg(text: str, categories: list[str]) -> dict[str, float]:
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith(("#", "@"))]
    if len(lines) != 1:
        raise ValueError(f"GROMACS energy output must have one frame, got {len(lines)}")
    values = [float(x) for x in lines[0].split()]
    if len(values) != len(categories) + 1 or not np.all(np.isfinite(values)):
        raise ValueError("GROMACS energy columns missing or non-finite")
    if abs(values[0]) > 1e-9 or len(set(categories)) != len(categories):
        raise ValueError("GROMACS energy time or category mapping is invalid")
    return dict(zip(categories, values[1:]))


def energies(ir, positions_nm, topology_files, workdir, *, executable=None):
    executable = executable or find_gromacs()
    if executable is None:
        raise FileNotFoundError("gmx not installed")
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    version = _run(executable, ["--version"], workdir, "version")
    version_text = version.stdout + version.stderr
    precision = re.search(r"^Precision:\s*(\S+)", version_text, re.M)
    if precision is None:
        raise ValueError("GROMACS build precision is missing in version output")

    pos = np.asarray(positions_nm, dtype=np.float64)
    reference = evaluate(ir, pos)
    shift = 0.05 - pos.min(axis=0)
    pos = pos + shift
    translated = evaluate(ir, pos)
    if any(abs(reference[k] - translated[k]) > 1e-9 for k in reference):
        raise ValueError("rigid translation changed reference energy")
    span = float(np.max(np.linalg.norm(pos[:, None] - pos[None, :], axis=2)))
    cutoff = max(2.0, span + 0.5)
    box = 2.0 * cutoff + 2.0
    # All intramolecular pairs fit inside cutoff, all images lie outside it.
    lines = [f"ffcompiler {ir.name} t= 0", str(len(ir.atoms))]
    for n, (atom, xyz) in enumerate(zip(ir.atoms, pos), 1):
        lines.append(f"{1:5d}{atom.residue_name[:5]:<5}{atom.name[:5]:>5}{n:5d}"
                     + "".join(f"{v:15.10f}" for v in xyz))
    lines.append(" ".join(f"{box:.10f}" for _ in range(3)))
    (workdir / "frame.gro").write_text("\n".join(lines) + "\n")
    mdp = f"""integrator = md
nsteps = 0
dt = 0.001
continuation = yes
cutoff-scheme = Verlet
nstlist = 1
verlet-buffer-tolerance = -1
rlist = {cutoff:.10f}
rcoulomb = {cutoff:.10f}
rvdw = {cutoff:.10f}
coulombtype = Cut-off
coulomb-modifier = None
vdwtype = Cut-off
vdw-modifier = None
DispCorr = no
pbc = xyz
constraints = none
tcoupl = no
pcoupl = no
comm-mode = none
gen-vel = no
nstenergy = 1
nstcalcenergy = 1
nstlog = 1
"""
    (workdir / "singlepoint.mdp").write_text(mdp)
    _run(executable, ["grompp", "-f", "singlepoint.mdp", "-c", "frame.gro", "-p",
                      str(Path(topology_files["top"]).resolve()), "-o", "run.tpr",
                      "-po", "processed.mdp"], workdir, "grompp")
    _run(executable, ["mdrun", "-s", "run.tpr", "-rerun", "frame.gro", "-deffnm", "rerun",
                      "-ntmpi", "1", "-ntomp", "1", "-nb", "cpu", "-bonded", "cpu"], workdir, "mdrun")
    menu = _run(executable, ["energy", "-f", "rerun.edr", "-o", "menu.xvg", "-dp"],
                workdir, "energy_menu", stdin="0\n", check=False)
    selected = sorted(set((int(n), label) for n, label in re.findall(
        r"(?<!\S)(\d+)\s+([A-Za-z][\w()./\-]+)", menu.stdout + menu.stderr) if label in TERMS))
    if not any(label == "Potential" for _, label in selected):
        raise ValueError("GROMACS energy menu lacks Potential")
    _run(executable, ["energy", "-f", "rerun.edr", "-o", "energy.xvg", "-dp"],
         workdir, "energy", stdin=" ".join(str(n) for n, _ in selected) + "\n0\n")
    observed = parse_xvg((workdir / "energy.xvg").read_text(), [TERMS[label] for _, label in selected])
    for category in reference:
        if category not in observed:
            if abs(reference[category]) > 1e-12:
                raise ValueError(f"GROMACS lacks nonzero energy category {category}")
            observed[category] = 0.0
    metadata = {"executable": executable, "precision": precision.group(1),
                "version": version_text, "rigid_translation_nm": shift.tolist(),
                "max_pair_distance_nm": span, "cutoff_nm": cutoff, "box_nm": box,
                "tolerance_policy": (
                    "mixed_fp32 (decision 0013): 1e-4 + 1e-5*abs(E_ref) kJ/mol"
                    if precision.group(1) == "mixed"
                    else "strict_file: 1e-4 + 1e-7*abs(E_ref) kJ/mol")}
    (workdir / "evaluation.json").write_text(json.dumps({"energies": observed, **metadata}, indent=2))
    return observed, metadata
