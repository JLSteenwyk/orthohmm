from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.reproduce_swiss_comparators import compare_results


def fixture():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_swiss_comparator_intervals_20260917.json"
    return json.loads(path.read_text())


def test_relocation_only_is_allowed():
    expected = fixture()
    observed = deepcopy(expected)
    for item in [observed[k] for k in ("counts", "protocol", "source")] + observed["helpers"]:
        item["path"] = "/new/location/" + Path(item["path"]).name
    compare_results(expected, observed)


@pytest.mark.parametrize("mutation", ["interval", "family", "version", "input", "helper", "extra"])
def test_reproduction_cannot_hide_changed_evidence(mutation):
    expected = fixture()
    observed = deepcopy(expected)
    if mutation == "interval":
        observed["comparisons"][0]["metrics"]["F1"]["bonferroni_percentile_ci"][0] += .001
    elif mutation == "family":
        observed["comparisons"][0]["family_differences"].pop()
    elif mutation == "version":
        observed["numpy_version"] = "different"
    elif mutation == "input":
        observed["counts"]["sha256"] = "0"*64
    elif mutation == "helper":
        observed["helpers"][0]["bytes"] += 1
    else:
        observed["unexpected_field"] = True
    with pytest.raises(ValueError):
        compare_results(expected, observed)
