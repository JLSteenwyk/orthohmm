import copy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools.export_current_benchmark_scores import SOURCES, METRICS, assemble, export

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def reports():
    return [json.loads((RESULTS / name).read_text()) for name in SOURCES]


def test_current_sources_and_exact_units(tmp_path):
    ob, qfo, kingdoms = reports()
    result = export(RESULTS, tmp_path / "scores")
    assert len(result["rows"]) == 8
    for row, old, corrected in zip(result["rows"], ob["methods"], qfo["methods"]):
        assert row["scores"]["OrthoBench"] == old["orthobench"]["f_score_percent"] / 100
        assert all(row["scores"][metric] == corrected["scores"][metric] for metric in METRICS)
    sonic = next(r for r in result["rows"] if r["key"] == "sonicparanoid_2_0_9")
    assert sonic["three_kingdoms"]["run"] == "contemporary matched input"
    selected = next(r for r in kingdoms["rows"] if r["key"] == sonic["key"])
    assert sonic["scores"]["ThreeKingdoms"] == selected["counts"]["f_score"]
    with (tmp_path / "scores/scores.tsv").open() as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 8
    with pytest.raises(FileExistsError):
        export(RESULTS, tmp_path / "scores")


@pytest.mark.parametrize("change", ["original", "duplicate", "missing", "admission", "mean", "nan", "diagnostic"])
def test_rejects_mixed_or_invalid_inputs(change):
    ob, qfo, kingdoms = reports()
    if change == "original":
        qfo["status"] = "historical"
    elif change == "duplicate":
        qfo["methods"].append(copy.deepcopy(qfo["methods"][0]))
    elif change == "missing":
        ob["methods"].pop()
    elif change == "admission":
        qfo["methods"][0]["status"] = "missing"
    elif change == "mean":
        qfo["methods"][0]["secondary_mean"] = 0
    elif change == "nan":
        qfo["methods"][0]["scores"]["GO"] = float("nan")
    else:
        kingdoms["rows"][-1]["use"] = "comparison"
    with pytest.raises(ValueError):
        assemble(ob, qfo, kingdoms)


def test_wrong_source_hash_fails_before_output(tmp_path):
    source = tmp_path / "inputs"
    source.mkdir()
    (source / next(iter(SOURCES))).write_text("{}")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Changed retained"):
        export(source, output)
    assert not output.exists()
