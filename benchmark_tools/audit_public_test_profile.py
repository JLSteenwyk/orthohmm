"""Record and check the explicit public/private test boundary without running tests."""

import argparse
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import platform
import sys

import pytest


RAW_NODES = frozenset({
    "tests/test_export_swiss_duplication_strata.py::test_export_and_no_overwrite",
    "tests/test_export_swiss_duplication_strata.py::test_independent_rational_reproduction_of_export",
    "tests/unit/test_swiss_strata_count_selection.py::test_explicit_counts_preserve_previous_rows_and_reject_wrong_hash[fragment]",
    "tests/unit/test_swiss_strata_count_selection.py::test_explicit_counts_preserve_previous_rows_and_reject_wrong_hash[duplication]",
})


class CollectionAudit:
    def __init__(self):
        self.collected = []
        self.marked = []
        self.selected = []
        self.deselected = []
        self.errors = []
        self.skipped = []

    @pytest.hookimpl(tryfirst=True)
    def pytest_collection_modifyitems(self, items):
        self.collected = [item.nodeid for item in items]
        self.marked = [item.nodeid for item in items if item.get_closest_marker("raw_benchmark")]

    def pytest_deselected(self, items):
        self.deselected.extend(item.nodeid for item in items)

    def pytest_collection_finish(self, session):
        self.selected = [item.nodeid for item in session.items]

    def pytest_collectreport(self, report):
        if report.failed:
            self.errors.append(str(report.longrepr))
        elif report.skipped:
            self.skipped.append(dict(nodeid=report.nodeid, reason=str(report.longrepr)))


def evaluate(collection, exit_code):
    problems = []
    lists = (collection.collected, collection.marked, collection.selected, collection.deselected)
    if exit_code != 0 or collection.errors:
        problems.append("collection did not complete successfully")
    if collection.skipped:
        problems.append("collection skipped modules; public inventory is incomplete")
    if any(len(values) != len(set(values)) for values in lists):
        problems.append("duplicate node identifiers")
    if set(collection.marked) != RAW_NODES:
        problems.append("raw marker inventory differs from the four declared private cases")
    if set(collection.deselected) != RAW_NODES:
        problems.append("deselection differs from the declared private cases")
    if not collection.selected or set(collection.selected) != set(collection.collected) - RAW_NODES:
        problems.append("public selection is not the complete collected inventory minus private cases")
    if set(collection.selected) & RAW_NODES:
        problems.append("private cases were selected")
    return dict(
        schema="orthohmm-public-test-profile-v1",
        profile="public (not raw_benchmark); collection only, no slow filtering",
        complete_raw_benchmark_gate=False,
        tests_executed=False,
        pytest_exit_code=int(exit_code),
        collected_count=len(collection.collected),
        selected_count=len(collection.selected),
        deselected_count=len(collection.deselected),
        declared_private_nodes=sorted(RAW_NODES),
        raw_marked_nodes=sorted(collection.marked),
        selected_nodes=sorted(collection.selected),
        deselected_nodes=sorted(collection.deselected),
        collection_errors=collection.errors,
        collection_skips=collection.skipped,
        problems=problems,
        passed=not problems,
    )


def audit(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    if Path.cwd().resolve() != root:
        raise ValueError("Run from the checkout root so collected node IDs remain relative")
    collection = CollectionAudit()
    log = io.StringIO()
    arguments = ["tests", "--collect-only", "-q", "-m", "not raw_benchmark"]
    with redirect_stdout(log), redirect_stderr(log):
        status = pytest.main(arguments, plugins=[collection])
    report = evaluate(collection, status)
    report.update(checkout_root=str(root), python=platform.python_version(),
                  python_executable=sys.executable, pytest_version=pytest.__version__,
                  arguments=arguments, collection_log=log.getvalue())
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.root, args.output)
    print(json.dumps({key: report[key] for key in
                      ("profile", "collected_count", "selected_count", "deselected_count", "passed", "problems")}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
