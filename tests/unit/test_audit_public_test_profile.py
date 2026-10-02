import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.audit_public_test_profile import CollectionAudit, RAW_NODES, audit, evaluate


ROOT = Path(__file__).resolve().parents[2]
SOURCE_TESTS = ["tests/test_export_swiss_duplication_strata.py",
                "tests/unit/test_swiss_strata_count_selection.py"]


def complete_collection():
    result = CollectionAudit()
    result.marked = result.deselected = sorted(RAW_NODES)
    result.selected = ["tests/public.py::test_public"]
    result.collected = result.selected + result.marked
    return result


def test_explicit_public_boundary_is_not_execution_or_full_regression():
    result = evaluate(complete_collection(), 0)
    assert result["passed"] is True and result["deselected_count"] == 4
    assert result["tests_executed"] is result["complete_raw_benchmark_gate"] is False
    assert result["declared_private_nodes"] == result["deselected_nodes"] == sorted(RAW_NODES)


@pytest.mark.parametrize("problem", ["unmarked", "extra_mark", "wrong_exclusion", "empty", "missing_public",
                                   "private_selected", "duplicate", "error", "skip", "exit"])
def test_changed_inventory_fails_closed(problem):
    collection = complete_collection()
    # Keep lists independent when testing mutations.
    collection.marked = collection.marked.copy()
    if problem == "unmarked":
        collection.marked.pop()
    elif problem == "extra_mark":
        collection.marked.append(collection.selected[0])
    elif problem == "wrong_exclusion":
        collection.deselected = collection.deselected[:-1]
    elif problem == "empty":
        collection.selected = []
    elif problem == "missing_public":
        collection.collected.append("tests/public.py::test_other")
    elif problem == "private_selected":
        collection.selected.append(collection.marked[0])
    elif problem == "duplicate":
        collection.collected.append(collection.collected[0])
    elif problem == "error":
        collection.errors.append("import failed")
    elif problem == "skip":
        collection.skipped.append(dict(nodeid="tests/missing.py", reason="dependency missing"))
    report = evaluate(collection, 1 if problem == "exit" else 0)
    assert report["passed"] is False and report["problems"]


@pytest.mark.parametrize("expression", [None, "not raw_benchmark", "raw_benchmark"])
def test_real_collection_keeps_default_private_cases_and_never_executes(expression, tmp_path):
    output = tmp_path / "collection.json"
    program = """
import json, pathlib, sys
import pytest
from benchmark_tools.audit_public_test_profile import CollectionAudit
class NoExecution(CollectionAudit):
    def pytest_runtest_protocol(self, item, nextitem):
        raise AssertionError('Collection must not execute tests')
plugin = NoExecution()
args = json.loads(sys.argv[1])
code = pytest.main(args, plugins=[plugin])
pathlib.Path(sys.argv[2]).write_text(json.dumps(dict(code=int(code), **vars(plugin))))
"""
    args = SOURCE_TESTS + ["--collect-only", "-q"]
    if expression:
        args += ["-m", expression]
    result = subprocess.run([sys.executable, "-c", program, json.dumps(args), str(output)],
                            cwd=ROOT, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads(output.read_text())
    assert data["code"] == 0 and not data["errors"] and not data["skipped"]
    assert len(data["collected"]) == 14 and set(data["marked"]) == RAW_NODES
    if expression is None:
        assert len(data["selected"]) == 14 and not data["deselected"]
        assert RAW_NODES <= set(data["selected"])
    elif expression == "not raw_benchmark":
        assert len(data["selected"]) == 10 and set(data["deselected"]) == RAW_NODES
    else:
        assert set(data["selected"]) == RAW_NODES and len(data["deselected"]) == 10


def test_no_overwrite_before_collection(tmp_path, monkeypatch):
    output = tmp_path / "existing.json"
    output.write_text("retained")
    monkeypatch.setattr(pytest, "main", lambda *a, **k: pytest.fail("must not collect"))
    with pytest.raises(FileExistsError):
        audit(ROOT, output)
    assert output.read_text() == "retained"


def test_wrong_working_directory_before_collection(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(pytest, "main", lambda *a, **k: pytest.fail("must not collect"))
    with pytest.raises(ValueError, match="checkout root"):
        audit(ROOT, tmp_path / "receipt.json")
    assert not (tmp_path / "receipt.json").exists()


@pytest.mark.parametrize("target,expression", [
    ("test.public.fast", 'not slow and not raw_benchmark'),
    ("coverage.public.unit", 'not raw_benchmark'),
    ("coverage.public.integration", 'not raw_benchmark'),
])
def test_public_make_targets_are_named_and_record_execution(target, expression):
    text = (ROOT / "Makefile").read_text()
    commands = text.split(target + ":\n", 1)[1].split("\n\n", 1)[0]
    assert f'-m "{expression}"' in commands and "--junitxml=" in commands
    assert "--ignore=tests/integration" in commands or "pytest tests/integration" in commands


@pytest.mark.parametrize("target", ["test.unit", "test.integration", "test.fast", "coverage.unit", "coverage.integration"])
def test_existing_full_targets_do_not_exclude_private_cases(target):
    commands = (ROOT / "Makefile").read_text().split(target + ":\n", 1)[1].split("\n\n", 1)[0]
    assert "raw_benchmark" not in commands
    assert "$(PYTEST_ARGS)" in commands


def test_workflow_declares_profile_and_retains_selection_and_execution():
    text = (ROOT / ".github/workflows/ci.yml").read_text()
    assert text.count("python -m benchmark_tools.audit_public_test_profile") == 2
    assert "make test.public.fast" in text and "make test.public.coverage" in text
    assert "make test.fast\n" not in text and "make test.coverage\n" not in text
    assert text.count("${{ runner.temp }}/public-test-selection.json") == 2
    assert text.count("${{ runner.temp }}/public-unit.xml") == 2
    assert text.count("${{ runner.temp }}/public-integration.xml") == 2
    assert "public-coverage (includes slow; excludes four raw-source cases)" in text
