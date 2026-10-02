import json
from copy import deepcopy
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools import plot_matched_graph as module
from benchmark_tools.plot_matched_graph import plot


def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/matched_graph_scores_20260926/results.json"
    return json.loads(path.read_text())


def test_plot_preserves_all_conditions_and_aligned_axes():
    fig = plot(report(), replay_policy="count-level")
    assert len(fig.axes) == 2
    assert len(fig.axes[0].get_yticklabels()) == 8
    assert fig.axes[0].get_ylim() == fig.axes[1].get_ylim()
    fig.canvas.draw()
    plt.close(fig)


def test_exact_mode_accepts_exact_current_runtime_report():
    result = report()
    computed = module.summarize(result["records"])
    result.update({key: computed[key] for key in ("contrasts", "bootstrap")})
    validation = module.validate_replay(result)
    assert validation["policy"] == "exact"
    assert validation["absolute_tolerance"] == 0
    assert validation["exact_differences"] == []
    fig = plot(result)
    plt.close(fig)


def test_reject_altered_point_estimate():
    result = report()
    result["contrasts"]["overall"]["f1"]["hmm_mean"] = .99
    with pytest.raises(ValueError):
        plot(result)


@pytest.mark.parametrize("direction", [-np.inf, np.inf])
def test_even_one_ulp_change_is_rejected_and_reported(monkeypatch, direction):
    result = report()
    expected = deepcopy(result)
    old = result["contrasts"]["overall"]["f1"]["hmm_mean"]
    changed = float(np.nextafter(old, direction))
    result["contrasts"]["overall"]["f1"]["hmm_mean"] = changed
    monkeypatch.setattr(module, "summarize", lambda _: expected)
    with pytest.raises(ValueError) as error:
        plot(result)
    differences = json.loads(str(error.value).split(": ", 1)[1])
    assert differences == [dict(path=["contrasts", "overall", "f1", "hmm_mean"],
                                reported=changed, recomputed=old)]


@pytest.mark.parametrize("reported,recomputed,expected", [
    ({"a": 1}, {"a": 1}, []),
    ({"a": [1, 2]}, {"a": [1, 3]}, [dict(path=["a", 1], reported=2, recomputed=3)]),
    ({"a": [1]}, {"a": [1, 2]}, [dict(path=["a"], reported=[1], recomputed=[1, 2])]),
    ({"a": None}, {}, [dict(path=["a"], reported_present=True, recomputed_present=False,
                           reported=None, recomputed=None)]),
    ({}, {"a": None}, [dict(path=["a"], reported_present=False, recomputed_present=True,
                           reported=None, recomputed=None)]),
    ({"a": "old"}, {"a": "new"}, [dict(path=["a"], reported="old", recomputed="new")]),
    ({"a": []}, {"a": {}}, [dict(path=["a"], reported=[], recomputed={})]),
])
def test_exact_diagnostics_cover_nested_values_and_structure(reported, recomputed, expected):
    assert module._replay_differences(reported, recomputed) == expected


def test_metadata_mismatch_is_not_suppressed(monkeypatch):
    result = report()
    expected = deepcopy(result)
    result["bootstrap"]["numpy_version"] = "different"
    monkeypatch.setattr(module, "summarize", lambda _: expected)
    with pytest.raises(ValueError) as error:
        plot(result)
    differences = json.loads(str(error.value).split(": ", 1)[1])
    assert differences == [dict(path=["bootstrap", "numpy_version"], reported="different",
                                recomputed=expected["bootstrap"]["numpy_version"])]


def test_count_level_uses_independent_counts_and_records_observed_remote_drift(monkeypatch):
    result = report()
    expected = deepcopy(result)
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/ci_exact_replay_difference_20261001.json"
    difference = json.loads(path.read_text())["guard_differences"][0]
    expected["contrasts"]["uneven_taxa"]["recall"]["marginal_95_percent_ci"][0] = difference["recomputed"]
    monkeypatch.setattr(module, "summarize", lambda _: expected)
    validation = module.validate_replay(result, "count-level")
    assert validation["exact_differences"] == [difference]
    assert validation["absolute_tolerance"] == 1e-12
    assert validation["relative_tolerance"] == 0
    assert validation["count_validation"]["metric_contrasts"] == 24
    assert validation["count_validation"]["f1_adjusted_intervals"] == 8
    assert result == report()
    with pytest.raises(ValueError):
        module.validate_replay(result)


@pytest.mark.parametrize("condition", [*module.CONDITIONS, "overall"])
@pytest.mark.parametrize("metric", ["f1", "precision", "recall"])
def test_count_level_rejects_changed_numeric_effects(condition, metric):
    result = report()
    result["contrasts"][condition][metric]["difference_percentage_points"] += 1e-9
    with pytest.raises(ValueError):
        module.validate_replay(result, "count-level")


@pytest.mark.parametrize("condition", [*module.CONDITIONS, "overall"])
def test_count_level_rejects_changed_adjusted_f1_interval(condition):
    result = report()
    result["contrasts"][condition]["f1"]["bonferroni_8_ci"][0] += 1e-9
    with pytest.raises(ValueError):
        module.validate_replay(result, "count-level")


@pytest.mark.parametrize("mutation", ["metadata", "replicates", "role", "wins", "extra", "shape",
                                      "nan", "inf", "count", "missing", "duplicate"])
def test_count_level_still_rejects_metadata_structure_and_count_corruption(mutation):
    result = report()
    row = result["contrasts"]["overall"]["f1"]
    if mutation == "metadata":
        result["bootstrap"]["numpy_version"] = "different"
    elif mutation == "replicates":
        result["bootstrap"]["replicates"] += 1
    elif mutation == "role":
        row["inference_role"] = "exploratory"
    elif mutation == "wins":
        row["wins"] = float(np.nextafter(float(row["wins"]), np.inf))
    elif mutation == "extra":
        row["unexpected"] = 0.0
    elif mutation == "shape":
        row["bonferroni_8_ci"].append(row["bonferroni_8_ci"][-1])
    elif mutation in ("nan", "inf"):
        row["hmm_mean"] = float(mutation)
    elif mutation == "count":
        result["records"][0]["score"]["tp"] += 1
    elif mutation == "missing":
        result["records"].pop()
    else:
        result["records"].append(deepcopy(result["records"][0]))
    with pytest.raises(ValueError):
        module.validate_replay(result, "count-level")


def test_unknown_policy_fails_before_drawing(monkeypatch):
    monkeypatch.setattr(module, "_draw", lambda _: pytest.fail("Invalid policy reached drawing"))
    with pytest.raises(ValueError, match="Unknown replay policy"):
        plot(report(), replay_policy="relaxed")


def test_count_failure_cannot_be_bypassed_by_equal_summary(monkeypatch):
    result = report()
    monkeypatch.setattr(module, "summarize", lambda _: result)
    result["records"][0]["score"]["tp"] += 1
    with pytest.raises(ValueError):
        module.validate_replay(result, "count-level")


def test_count_level_cli_export_records_policy_and_keeps_input_bytes(tmp_path):
    source = tmp_path / "results.json"
    source.write_text(json.dumps(report()))
    before = source.read_bytes()
    output = tmp_path / "export"
    subprocess.run([sys.executable, "-m", "benchmark_tools.plot_matched_graph", "--results", str(source),
                    "--output", str(output), "--replay-policy", "count-level"],
                   cwd=Path(__file__).resolve().parents[2], check=True, timeout=120, capture_output=True, text=True)
    manifest = json.loads((output / "manifest.json").read_text())
    validation = manifest["replay_validation"]
    assert validation["policy"] == "count-level"
    assert validation["absolute_tolerance"] == 1e-12
    assert validation["count_validation"]["metric_contrasts"] == 24
    assert manifest["input"] == module.record(source)
    assert manifest["visual_review_complete"] is False
    for item in manifest["outputs"]:
        assert item == module.record(item["path"])
    assert source.read_bytes() == before
    receipt = (output / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError):
        module.run(source, output, replay_policy="count-level")
    assert (output / "manifest.json").read_bytes() == receipt


def test_input_change_during_export_never_gets_success_manifest(tmp_path, monkeypatch):
    source = tmp_path / "results.json"
    source.write_text(json.dumps(report()))
    output = tmp_path / "export"
    original = plt.Figure.savefig
    def save_and_change(self, *args, **kwargs):
        source.write_text(source.read_text() + "\n")
        return original(self, *args, **kwargs)
    monkeypatch.setattr(plt.Figure, "savefig", save_and_change)
    figures = plt.get_fignums()
    with pytest.raises(ValueError, match="Input changed during rendering"):
        module.run(source, output, replay_policy="count-level")
    assert not (output / "manifest.json").exists()
    assert plt.get_fignums() == figures
