import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.report_qfo_recovered_stages import CONTRASTS, STAGES, markdown, summarize


@pytest.fixture
def admission():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_recovered_assessment_admitted_20260917.json"
    return json.loads(path.read_text())


def test_complete_fixed_contrasts(admission):
    report = summarize(admission)
    assert len(report["contrasts"]) == 4
    for row, (a, b) in zip(report["contrasts"], CONTRASTS):
        assert (row["candidate"], row["reference"]) == (STAGES[a], STAGES[b])
        assert len(row["differences"]) == 6
        assert row["paired_uncertainty"] == "not_yet_established"
    text = markdown(report)
    assert "Native Coordinates" in text
    assert all(stage in text for stage in STAGES)
    assert not report["publication_ready"]


@pytest.mark.parametrize("change", ["order", "failed", "nan", "missing", "mean", "counts"])
def test_reject_invalid_report(admission, change):
    data = copy.deepcopy(admission)
    row = data["records"][0]
    if change == "order":
        data["records"].reverse()
    elif change == "failed":
        row["status"] = "failed"
    elif change == "nan":
        row["assessment"]["endpoints"]["GO"]["score"] = float("nan")
    elif change == "missing":
        del row["assessment"]["endpoints"]["GO"]
    elif change == "mean":
        row["assessment"]["secondary_six_metric_mean"] += 0.1
    else:
        row["conversion"]["retained_pairs"] += 1
    with pytest.raises(ValueError):
        summarize(data)
