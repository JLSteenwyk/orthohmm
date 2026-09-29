from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.audit_publication_figures import validate_ob_complete, validate_ob_intervals


def manifest():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/ob_complete_strata_figure_20260928/manifest.json").read_text())


def test_complete_panel():
    validate_ob_complete(manifest())


@pytest.mark.parametrize("key,value", [("displayed_finite_values", 287), ("missing_values", 0), ("new_intervals", True)])
def test_reject_partial_or_inferential_panel(key, value):
    data = deepcopy(manifest())
    data[key] = value
    with pytest.raises(ValueError):
        validate_ob_complete(data)


def test_reject_unpinned_source():
    data = manifest()
    data["source"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_ob_complete(data)


def interval_manifest():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/figures_ob_complete_uncertainty_20260928/manifest.json").read_text())


def test_complete_intervals():
    validate_ob_intervals(interval_manifest())


@pytest.mark.parametrize("key,value", [("endpoints", 6), ("units", "proportion"),
    ("new_statistics", True), ("publication_ready", True), ("status", "incomplete")])
def test_reject_wrong_interval_scope(key, value):
    data = interval_manifest()
    data[key] = value
    with pytest.raises(ValueError):
        validate_ob_intervals(data)


def test_reject_wrong_interval_source():
    data = interval_manifest()
    data["inputs"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_ob_intervals(data)
