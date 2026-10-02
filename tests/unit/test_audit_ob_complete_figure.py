from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.audit_publication_figures import validate_ob_complete, validate_ob_intervals


def manifest(retained_record_at_path):
    root = Path(__file__).resolve().parents[2]
    data = json.loads((root / "benchmark_tools/results/ob_complete_strata_figure_20260928/manifest.json").read_text())
    data["source"] = retained_record_at_path(
        data["source"], root / "benchmark_tools/results/ob_complete_strata_20260928/report.json")
    return data


def test_complete_panel(retained_record_at_path):
    validate_ob_complete(manifest(retained_record_at_path))


@pytest.mark.parametrize("key,value", [("displayed_finite_values", 287), ("missing_values", 0), ("new_intervals", True)])
def test_reject_partial_or_inferential_panel(key, value, retained_record_at_path):
    data = deepcopy(manifest(retained_record_at_path))
    data[key] = value
    with pytest.raises(ValueError):
        validate_ob_complete(data)


def test_reject_unpinned_source(retained_record_at_path):
    data = manifest(retained_record_at_path)
    data["source"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_ob_complete(data)


@pytest.mark.parametrize("mutation", ["same_size", "extra_bytes"])
def test_complete_panel_rejects_changed_source_bytes(tmp_path, mutation, retained_record_at_path):
    data = manifest(retained_record_at_path)
    root = Path(__file__).resolve().parents[2]
    content = (root / "benchmark_tools/results/ob_complete_strata_20260928/report.json").read_bytes()
    changed = content.replace(b"  ", b"\t ", 1) if mutation == "same_size" else content + b"\n"
    assert changed != content
    assert json.loads(changed) == json.loads(content)
    source = tmp_path / "report.json"
    source.write_bytes(changed)
    data["source"]["path"] = str(source)
    with pytest.raises(ValueError, match="Changed complete descriptive result"):
        validate_ob_complete(data)


def interval_manifest(retained_record_at_path):
    root = Path(__file__).resolve().parents[2]
    data = json.loads((root / "benchmark_tools/results/figures_ob_complete_uncertainty_20260928/manifest.json").read_text())
    data["inputs"][0] = retained_record_at_path(
        data["inputs"][0], root / "benchmark_tools/results/ob_complete_uncertainty_20260928.json")
    return data


def test_complete_intervals(retained_record_at_path):
    validate_ob_intervals(interval_manifest(retained_record_at_path))


@pytest.mark.parametrize("key,value", [("endpoints", 6), ("units", "proportion"),
    ("new_statistics", True), ("publication_ready", True), ("status", "incomplete")])
def test_reject_wrong_interval_scope(key, value, retained_record_at_path):
    data = interval_manifest(retained_record_at_path)
    data[key] = value
    with pytest.raises(ValueError):
        validate_ob_intervals(data)


def test_reject_wrong_interval_source(retained_record_at_path):
    data = interval_manifest(retained_record_at_path)
    data["inputs"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_ob_intervals(data)


def test_complete_panel_rejects_wrong_source_size(retained_record_at_path):
    data = manifest(retained_record_at_path)
    data["source"]["bytes"] += 1
    with pytest.raises(ValueError, match="Changed complete descriptive result"):
        validate_ob_complete(data)
