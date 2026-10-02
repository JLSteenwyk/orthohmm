import importlib
from importlib import metadata
import json
from pathlib import Path

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("distribution,module", [
    ("biopython", "Bio"), ("psutil", "psutil"), ("dendropy", "dendropy"),
    ("ijson", "ijson"), ("openpyxl", "openpyxl"), ("matplotlib", "matplotlib"),
    ("pillow", "PIL"), ("PyMuPDF", "fitz"), ("scipy", "scipy"),
    ("setuptools", "setuptools"), ("sqlglot", "sqlglot"), ("numpy", "numpy"),
])
def test_workflow_test_dependencies_are_declared_and_importable(distribution, module):
    rows = [Requirement(line) for line in (ROOT / "tests/requirements.txt").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    selected = [r for r in rows if canonicalize_name(r.name) == canonicalize_name(distribution)]
    assert len(selected) == 1
    assert selected[0].specifier.contains(metadata.version(distribution))
    assert importlib.import_module(module) is not None


@pytest.mark.parametrize("relative,nested", [
    ("qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json", False),
    ("qfo_sequence_swiss_bootstrap_20260918.json", False),
    ("matched_graph_scores_20260926/results.json", True),
])
def test_numpy_pin_matches_retained_numerical_reports(relative, nested):
    rows = [Requirement(line) for line in (ROOT / "tests/requirements.txt").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    pin = next(r for r in rows if canonicalize_name(r.name) == "numpy")
    report = json.loads((ROOT / "benchmark_tools/results" / relative).read_text())
    version = (report["bootstrap"] if nested else report)["numpy_version"]
    assert str(pin.specifier) == "==" + version


@pytest.mark.parametrize("target", ["test.unit", "test.fast", "coverage.unit"])
def test_unit_make_targets_do_not_omit_top_level_cases(target):
    text = (ROOT / "Makefile").read_text()
    commands = text.split(target + ":\n", 1)[1].split("\n\n", 1)[0]
    assert "-m pytest tests --ignore=tests/integration" in commands
    assert "-m pytest tests/unit " not in commands
