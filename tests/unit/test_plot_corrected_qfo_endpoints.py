import json
from pathlib import Path

import pytest

from benchmark_tools.plot_corrected_qfo_endpoints import plotting_data, render

SOURCE = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json"


def test_complete_points_match_retained_details():
    report = json.loads(SOURCE.read_text())
    data = plotting_data(report)
    assert sum(len(row["qfo"]) for row in data) == 48
    for row, original in zip(data, report["methods"]):
        assert row["prediction_semantics"] == original["prediction_semantics"]
        assert row["qfo"]["SwissTrees F"]["x"] == original["details"]["SwissTrees"]["recall"]
        assert row["qfo"]["FAS"]["y"] == original["scores"]["FAS"]


@pytest.mark.parametrize("change", ["status", "inventory", "admission", "f1", "axis", "nonfinite"])
def test_invalid_endpoint_evidence_rejected(change):
    report = json.loads(SOURCE.read_text())
    row = report["methods"][0]
    if change == "status":
        report["status"] = "historical"
    elif change == "inventory":
        report["methods"].pop()
    elif change == "admission":
        row["status"] = "missing"
    elif change == "f1":
        row["scores"]["SwissTrees"] = 0
    elif change == "axis":
        row["details"]["GO"]["statistic"] = "F1"
    else:
        row["scores"]["FAS"] = float("nan")
    with pytest.raises(ValueError):
        plotting_data(report)


def test_render_and_refuse_overwrite(tmp_path):
    output = tmp_path / "figure"
    manifest = render(SOURCE, output)
    assert len(manifest["outputs"]) == 3
    assert all(item["bytes"] > 1000 for item in manifest["outputs"])
    assert not manifest["publication_ready"]
    with pytest.raises(FileExistsError):
        render(SOURCE, output)
