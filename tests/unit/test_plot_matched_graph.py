import json
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools import plot_matched_graph as module
from benchmark_tools.plot_matched_graph import plot


def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/matched_graph_scores_20260926/results.json"
    return json.loads(path.read_text())


def test_plot_preserves_all_conditions_and_aligned_axes():
    fig = plot(report())
    assert len(fig.axes) == 2
    assert len(fig.axes[0].get_yticklabels()) == 8
    assert fig.axes[0].get_ylim() == fig.axes[1].get_ylim()
    fig.canvas.draw()
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
