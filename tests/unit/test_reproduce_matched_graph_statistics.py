import json
from pathlib import Path

import pytest

from benchmark_tools.reproduce_matched_graph_statistics import reproduce


def result():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/matched_graph_scores_20260926/results.json"
    return json.loads(path.read_text())


def test_reproduction():
    reproduced = reproduce(result())
    assert reproduced["metric_contrasts"] == 24 and reproduced["f1_adjusted_intervals"] == 8


@pytest.mark.parametrize("mutation", ["count", "effect", "interval", "missing", "duplicate", "wins"])
def test_corruption_rejected(mutation):
    data = result()
    if mutation == "count":
        data["records"][0]["score"]["tp"] += 1
    elif mutation == "effect":
        data["contrasts"]["overall"]["f1"]["difference_percentage_points"] += .1
    elif mutation == "interval":
        data["contrasts"]["divergent"]["f1"]["bonferroni_8_ci"][0] += .1
    elif mutation == "missing":
        data["records"].pop()
    elif mutation == "duplicate":
        data["records"].append(data["records"][0])
    else:
        data["contrasts"]["overall"]["f1"]["wins"] = 0
    with pytest.raises(ValueError):
        reproduce(data)
