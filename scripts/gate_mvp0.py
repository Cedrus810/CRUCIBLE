#!/usr/bin/env python
"""Gate G0：全部固件 × 4 构型的跨引擎等价验收（MVP-0 计划 T10）。

在 RunContext("gate_mvp0") 下运行，写 report.md（每个固件一节）和
acceptance.json。pytest 摘要经子进程收集（非 external 与 external 各一轮）。

    python scripts/gate_mvp0.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import numpy as np  # noqa: E402

from ffcompiler.importers import from_amber  # noqa: E402
from ffcompiler.runtime import RunContext  # noqa: E402
from ffcompiler.validation.engines import compare  # noqa: E402
from ffcompiler.validation.pytest_results import summarize_pytest  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "gaff2"
NAMES = (
    "ethanol", "benzene", "nma", "cyclopropane", "cyclobutane", "norbornane",
    "acetate", "methylammonium", "chlorobenzene", "methyl_benzoate", "dmso",
)  # 乙腈（近线性角）不进 Gate，只作诊断（决策 0013）
N_GEOMETRIES = 4
SEED = 20261004


def _pytest_summary(run_dir: Path) -> dict:
    summary = {}
    for marker in ("not external", "external"):
        stem = "pytest_unit" if marker == "not external" else "pytest_external"
        xml = (run_dir / f"{stem}.xml").resolve()
        proc = subprocess.run(
            ["bash", str(REPO_ROOT / "scripts/run_logged.sh"), sys.executable,
             "-m", "pytest", "-m", marker, "-q", "--durations=0", f"--junitxml={xml}"],
            capture_output=True, text=True, cwd=REPO_ROOT, timeout=3600,
        )
        (run_dir / f"{stem}.log").write_text(proc.stdout + proc.stderr)
        summary[marker] = summarize_pytest(xml, proc.returncode)
    return summary


def main() -> int:
    with RunContext("gate_mvp0", root=REPO_ROOT / "jobs",
                    seed=SEED, inputs=[
                        FIXTURE_DIR / n / f"{n}.{ext}" for n in NAMES
                        for ext in ("prmtop", "inpcrd", "sdf")
                    ]) as ctx:
        run_dir = ctx.dir
        reports: dict[str, list] = defaultdict(list)
        rng = np.random.default_rng(SEED)

        for name in NAMES:
            mdir = FIXTURE_DIR / name
            ir, positions = from_amber(
                mdir / f"{name}.prmtop", mdir / f"{name}.inpcrd", mdir / f"{name}.sdf",
            )
            configs = [positions] + [
                positions + rng.normal(0.0, 0.005, size=positions.shape)
                for _ in range(N_GEOMETRIES - 1)
            ]
            for n, pos in enumerate(configs):
                ctx.log.info("compare %s config %d", name, n)
                report = compare(
                    ir, pos, run_dir / "work" / f"{name}_cfg{n}",
                    tleap_prmtop=mdir / f"{name}.prmtop",
                    tleap_inpcrd=mdir / f"{name}.inpcrd",
                )
                reports[name].append(report)
            ctx.log.info("%s done: passed=%s", name, all(r.passed for r in reports[name]))

        # ---- report.md ---------------------------------------------------
        lines = ["# Gate G0 报告", ""]
        for name in NAMES:
            lines.append(f"## {name}")
            for n, report in enumerate(reports[name]):
                lines.append(f"### config {n}  passed={report.passed}")
                lines.append("")
                lines.append(report.to_markdown())
                lines.append("")
        (run_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")

        # ---- acceptance.json ---------------------------------------------
        legs: dict[str, dict[str, int]] = defaultdict(
            lambda: {"passed": 0, "failed": 0}
        )
        skipped_legs: dict[str, str] = {}
        max_rel_by_category: dict[str, float] = defaultdict(float)
        for name in NAMES:
            for report in reports[name]:
                for row in report.rows:
                    legs[row["leg"]]["passed" if row["passed"] else "failed"] += 1
                    if row["leg"] == "numpy-ir":
                        continue
                    rel = row["delta"] / max(abs(row["E_ref"]), 1.0)
                    max_rel_by_category[row["category"]] = max(
                        max_rel_by_category[row["category"]], rel,
                    )
                for row in report.force_rows:
                    legs[row["leg"]]["passed" if row["passed"] else "failed"] += 1
                skipped_legs.update(report.skipped)

        all_reports = [r for name in NAMES for r in reports[name]]
        acceptance = {
            "gate": "G0",
            "passed": all(r.passed for r in all_reports),
            "n_fixtures": len(NAMES),
            "n_geometries": N_GEOMETRIES,
            "legs": dict(legs),
            "max_rel_error_by_category": dict(max_rel_by_category),
            "skipped_legs": skipped_legs,
            "pytest_summary": {},
            "engine_errors": [r.errors for r in all_reports if r.errors],
            "gromacs_builds": [r.metadata["gmx"] for r in all_reports if "gmx" in r.metadata],
        }
        acceptance["pytest_summary"] = _pytest_summary(run_dir)
        acceptance["passed"] = acceptance["passed"] and all(
            s["accepted"] for s in acceptance["pytest_summary"].values())
        (run_dir / "acceptance.json").write_text(
            json.dumps(acceptance, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        ctx.log.info("G0 acceptance: passed=%s", acceptance["passed"])
        print(f"G0 passed={acceptance['passed']}  dir={run_dir}")
        rc = 0 if acceptance["passed"] else 1
        ctx.finish(rc)
        return rc


if __name__ == "__main__":
    raise SystemExit(main())
