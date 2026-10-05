"""Fail-closed pytest JUnit summaries for formal acceptance runs."""
from pathlib import Path
import xml.etree.ElementTree as ET


def summarize_pytest(path: Path, returncode: int) -> dict:
    summary = dict(passed=0, failed=0, errors=0, skipped=0,
                   unexpected_skips=0, returncode=returncode, accepted=False)
    try:
        cases = list(ET.parse(path).getroot().iter("testcase"))
        if not cases:
            raise ValueError("JUnit contains no testcases")
        for case in cases:
            if case.find("error") is not None:
                summary["errors"] += 1
            elif case.find("failure") is not None:
                summary["failed"] += 1
            elif (skip := case.find("skipped")) is not None:
                summary["skipped"] += 1
                allowed = (case.get("classname") == "tests.unit.test_export_gromacs"
                           and case.get("name") == "test_real_gromacs_rerun"
                           and skip.get("message") == "gmx not installed")
                summary["unexpected_skips"] += not allowed
            else:
                summary["passed"] += 1
        summary["accepted"] = (returncode == 0 and summary["passed"] > 0
                               and summary["failed"] == summary["errors"] == 0
                               and summary["unexpected_skips"] == 0)
    except (OSError, ET.ParseError, ValueError) as exc:
        summary["parse_error"] = str(exc)
    return summary
