import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import export_three_kingdoms_comparison as module

BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def documents():
    return [json.loads((BASE / name).read_text()) for name in module.SOURCES]


def test_explicit_matched_and_historical_sonic_rows():
    rows = module.assemble(*documents())
    assert len(rows) == 9
    primary = [r for r in rows if r["use"] == "comparison"]
    assert len(primary) == 8
    sonic = next(r for r in primary if r["key"] == module.SONIC)
    assert sonic["run"] == "contemporary matched input"
    assert sonic["counts"]["true_positive_gene_pairs"] == 7272
    assert rows[-1]["counts"]["true_positive_gene_pairs"] == 7265
    assert rows[-1]["use"] == "historical diagnostic only"
    assert all("historical consumption not proven" in r["input_status"]
               for r in primary if r["key"] != module.SONIC)


@pytest.mark.parametrize("problem", ["missing", "duplicate", "unadmitted", "reference", "inputs"])
def test_invalid_panel_rejected(problem):
    old, matched, inputs = copy.deepcopy(documents())
    if problem == "missing":
        old["methods"].pop()
    elif problem == "duplicate":
        old["methods"].append(old["methods"][0])
    elif problem == "unadmitted":
        matched["accuracy_admitted"] = False
    elif problem == "reference":
        matched["counts"]["reference_genes"] = 1
    else:
        inputs["methods"].pop()
    with pytest.raises(ValueError):
        module.assemble(old, matched, inputs)


def test_preserves_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        module.export(BASE.parents[1], tmp_path)


def test_retained_export_exact_rows_and_limits():
    result = json.loads((BASE / "three_kingdoms_comparison_matched_20260918/comparison.json").read_text())
    assert result["rows"] == module.assemble(*documents())
    assert result["publication_ready"] is False
    assert result["uniform_historical_input_consumption_proven"] is False
    assert any("not penalized" in value for value in result["limitations"])
