"""Freeze the tested prospective execution route; never submit or release jobs."""

import argparse
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import native_factorial_allocated_execution as contract
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save


TEST_NAMES = (
    "test_allocated_native_factorial_qfo", "test_allocated_native_factorial_contract",
    "test_run_allocated_native_factorial_cost", "test_review_allocated_native_factorial_attempt",
    "test_validate_allocated_native_factorial_outputs", "test_native_factorial_cost",
    "test_review_native_factorial_attempt", "test_prepare_native_factorial_request",
    "test_native_factorial_outputs", "test_allocated_threadripper_scaling",
    "test_native_factorial_adapter", "test_prepare_native_factorial_qfo_pairs",
    "test_run_native_factorial_qfo_assessment", "test_prepare_allocated_native_factorial_amendment")


def test_report(ref):
    check(ref)
    document = ET.parse(ref["path"]).getroot()
    if document.tag != "testsuites":
        raise ValueError("Require pytest JUnit testsuites report")
    suites = list(document)
    if not suites or any(s.tag != "testsuite" for s in suites):
        raise ValueError("Unexpected JUnit suite structure")
    cases = [case for suite in suites for case in suite.findall("testcase")]
    if (not cases or sum(int(s.attrib["tests"]) for s in suites) != len(cases)
            or any(int(s.attrib[k]) != 0 for s in suites for k in ("failures", "errors", "skipped"))
            or any(list(case) for case in cases)):
        raise ValueError("Integration report is failed, skipped, incomplete or malformed")
    modules = {case.attrib["classname"].removeprefix("tests.unit.").split(".")[0] for case in cases}
    if not set(TEST_NAMES).issubset(modules):
        raise ValueError("Integration report omits required execution or scientific handoff tests")
    check(ref)
    return dict(tests=len(cases), failures=0, errors=0, skipped=0, required_test_modules=list(TEST_NAMES),
        scope="synthetic integration and unchanged scientific kernels; not a live production inference pass")


def committed(ref, commit):
    path = Path(ref["path"])
    relative = path.relative_to(contract.ROOT)
    if path.read_bytes() != subprocess.check_output(
            ["git", "-C", str(contract.ROOT), "show", commit + ":" + str(relative)]):
        raise ValueError("Freeze requires committed source/test evidence: " + str(relative))
    check(ref)


def prepare(output, tests_ref):
    output = Path(output)
    if (not output.is_absolute() or output.resolve() != output or not output.is_relative_to(contract.ROOT)
            or output.exists() or output.is_symlink() or not output.parent.is_dir()):
        raise ValueError("Require a fresh direct amendment path in the repository")
    validation = test_report(tests_ref)
    new_sources = contract.sources()
    tests = [record(contract.ROOT / "tests/unit" / (name + ".py")) for name in TEST_NAMES]
    source_ref = record(__file__)
    commit = subprocess.check_output(["git", "-C", str(contract.ROOT), "rev-parse", "HEAD"], text=True).strip()
    for ref in [source_ref, *new_sources, *tests, tests_ref]:
        committed(ref, commit)
    for name in contract.SOURCE_NAMES:
        if name.endswith(".py"):
            importlib.import_module("benchmark_tools." + name[:-3])
    value = dict(schema="allocated_native_factorial_execution_amendment_v1", root=str(contract.ROOT),
        execution_scope=contract.SCOPE, allowed_indices=[10, 11, 12],
        placement_policy="actual_allocated_physical_cores_v1", resources=dict(contract.RESOURCES),
        historical_plan=record(contract.PLAN), placement_fixture=record(contract.FIXTURE),
        historical_prefix=contract.historical_prefix(), new_sources=new_sources,
        automatic_retry=False, scientific_method_modified=False, historical_evidence_modified=False,
        production_execution_authorized=True, qfo_conversion_scoring_ready=True,
        scientific_timings_admitted=False, publication_ready=False,
        source=source_ref, source_commit=commit, prepared_unix_ns=time.time_ns(),
        tests_report=tests_ref, test_sources=tests, validation=validation,
        authorization_basis="Active goal; prospectively bound allocation route for different unrun identities",
        timing_disclosure="Shared Threadripper observations with unknown, potentially tool-dependent contention; not isolated performance.",
        limitations=["Validated engineering fixture and synthetic integration are not a live production success.",
            "Held-job ownership, prefix outcomes, safe capacity and source/runtime accounting gates remain required.",
            "No inference rerun, scheduler change, scientific-result admission or publication readiness is authorized here."])
    contract.validate_amendment(value)
    for ref in [source_ref, *new_sources, *tests, tests_ref]:
        check(ref)
    save(output, value)
    return record(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tests", type=Path, required=True)
    parser.add_argument("--tests-sha256", required=True)
    args = parser.parse_args()
    ref = record(args.tests)
    if ref["sha256"] != args.tests_sha256:
        raise ValueError("Integration report checksum differs")
    print(json.dumps(prepare(args.output.absolute(), ref), sort_keys=True))
