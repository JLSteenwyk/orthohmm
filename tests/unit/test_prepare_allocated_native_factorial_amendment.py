import xml.etree.ElementTree as ET

import pytest

from benchmark_tools import prepare_allocated_native_factorial_amendment as preparation
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def report(tmp_path):
    root = ET.Element("testsuites")
    suite = ET.SubElement(root, "testsuite", tests=str(len(preparation.TEST_NAMES)), failures="0", errors="0", skipped="0")
    for name in preparation.TEST_NAMES:
        ET.SubElement(suite, "testcase", classname="tests.unit." + name, name="synthetic_test")
    path = tmp_path / "report.xml"
    ET.ElementTree(root).write(path)
    return path


def test_requires_every_route_and_scientific_handoff_test_module(report):
    verdict = preparation.test_report(record(report))
    assert verdict["tests"] == len(preparation.TEST_NAMES)
    assert "not a live production" in verdict["scope"]


@pytest.mark.parametrize("fault", ["root", "empty", "count", "failure", "error", "skip", "case_failure", "missing_module"])
def test_incomplete_or_failed_integration_evidence_not_frozen(report, fault):
    root = ET.parse(report).getroot()
    suite = root[0]
    if fault == "root": root.tag = "testsuite"
    elif fault == "empty": suite.clear()
    elif fault == "count": suite.set("tests", "9999")
    elif fault in {"failure", "error", "skip"}:
        suite.set({"failure":"failures", "error":"errors", "skip":"skipped"}[fault], "1")
    elif fault == "case_failure": ET.SubElement(suite[0], "failure")
    else: suite[0].set("classname", "tests.unit.unrelated_test")
    ET.ElementTree(root).write(report)
    with pytest.raises((ValueError, KeyError)):
        preparation.test_report(record(report))


def test_changed_junit_report_refuses_original_binding(report):
    ref = record(report)
    report.write_text(report.read_text() + "\n")
    with pytest.raises(ValueError):
        preparation.test_report(ref)


def test_freeze_does_not_overwrite_existing_path(report):
    with pytest.raises(ValueError, match="fresh direct"):
        preparation.prepare(report, record(report))


def test_freeze_never_accepts_relative_path(report):
    with pytest.raises(ValueError, match="fresh direct"):
        preparation.prepare("amendment.json", record(report))
