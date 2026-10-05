from pathlib import Path
import pytest
from ffcompiler.validation.pytest_results import summarize_pytest


@pytest.mark.parametrize("xml,rc,accepted", [
    ('<testsuite><testcase name="ok"/></testsuite>', 0, True),
    ('<testsuite><testcase name="ok"/></testsuite>', 1, False),
    ('<testsuite><testcase><failure/></testcase><testcase/></testsuite>', 1, False),
    ('<testsuite><testcase><error/></testcase></testsuite>', 2, False),
    ('<testsuite/>', 5, False),
    ('not xml', 0, False),
    ('<testsuite><testcase><skipped message="anything"/></testcase><testcase/></testsuite>', 0, False),
    ('<testsuite><testcase classname="tests.unit.test_export_gromacs" name="test_real_gromacs_rerun"><skipped message="gmx not installed"/></testcase><testcase/></testsuite>', 0, True),
])
def test_gate_fail_closed(tmp_path: Path, xml, rc, accepted):
    path = tmp_path / "pytest.xml"
    path.write_text(xml)
    assert summarize_pytest(path, rc)["accepted"] is accepted


def test_missing_results_fail(tmp_path):
    assert not summarize_pytest(tmp_path / "missing.xml", 0)["accepted"]
