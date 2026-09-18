import copy
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools.plot_qfo_factorial import plot, validate


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_factorial_swiss_bootstrap_20260918.json"
    return json.loads(path.read_text())


def test_units_and_panel_inventory(report):
    scores = validate(report)
    assert scores.shape == (8, 3)
    assert scores[0, 0] == pytest.approx(67.75670237994498)
    figure = plot(report)
    assert len(figure.axes) == 4
    assert len(figure.axes[0].texts) == 24
    for ax in figure.axes[1:]:
        assert len(ax.get_yticks()) == 14
    interval = figure.axes[1].lines[1].get_xdata()
    assert np.allclose(interval, np.asarray(report["comparisons"][0]["metrics"]["F1"]["bonferroni_percentile_ci"]) * 100)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for text in figure.texts:
        bounds = text.get_window_extent(renderer)
        assert figure.bbox.contains(bounds.x0, bounds.y0)
        assert figure.bbox.contains(bounds.x1, bounds.y1)
    plt.close(figure)


@pytest.mark.parametrize("field,value", [("replicates", 20000), ("seed", 1),
    ("multiplicity_endpoints", 36), ("units", "percent"), ("alpha", .1)])
def test_changed_protocol(report, field, value):
    report[field] = value
    with pytest.raises(ValueError):
        validate(report)


@pytest.mark.parametrize("change", ["missing_cell", "missing_contrast", "identity", "effect", "interval", "nonfinite", "family", "f1"])
def test_reject_invalid_results(report, change):
    report = copy.deepcopy(report)
    value = report["comparisons"][0]["metrics"]["F1"]
    if change == "missing_cell":
        report["point_estimates"].pop("p0_c0_r0")
    elif change == "missing_contrast":
        report["comparisons"].pop()
    elif change == "identity":
        report["comparisons"][0]["weights"] = [0] * 8
    elif change == "effect":
        value["difference"] += .1
    elif change == "interval":
        value["bonferroni_percentile_ci"] = [1, -1]
    elif change == "nonfinite":
        value["paired_percentile_ci"][0] = float("nan")
    elif change == "family":
        value["family_wins"] = 19
    elif change == "f1":
        report["point_estimates"]["p0_c0_r0"]["F1"] += .1
    with pytest.raises(ValueError):
        validate(report)
