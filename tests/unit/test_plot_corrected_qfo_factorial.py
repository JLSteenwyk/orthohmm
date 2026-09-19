import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools.plot_qfo_factorial import METRICS, plot, validate


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[2] / (
        "benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json")
    return json.loads(path.read_text())


def test_corrected_panels_and_every_interval(report):
    scores = validate(report, "corrected")
    assert scores.shape == (8, 3)
    assert scores[-1, 0] == pytest.approx(83.35132180095781)
    figure = plot(report, "corrected")
    try:
        labels = "\n".join(t.get_text() for t in figure.texts)
        assert "Corrected-release inputs (984,137 genes)" in labels
        assert "4/14 adjusted F1 intervals exclude zero" in labels
        assert "Original-release" not in labels
        assert len(figure.axes[0].texts) == 24
        for ax, metric in zip(figure.axes[1:], METRICS):
            assert len(ax.get_yticks()) == 14
            for i, row in enumerate(report["comparisons"]):
                value = row["metrics"][metric]
                assert np.allclose(ax.lines[1 + i * 3].get_xdata(),
                                   np.asarray(value["bonferroni_percentile_ci"]) * 100)
                assert np.allclose(ax.lines[2 + i * 3].get_xdata(),
                                   np.asarray(value["paired_percentile_ci"]) * 100)
                assert np.allclose(ax.lines[3 + i * 3].get_xdata(), value["difference"] * 100)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        for text in figure.texts:
            bounds = text.get_window_extent(renderer)
            assert figure.bbox.contains(bounds.x0, bounds.y0)
            assert figure.bbox.contains(bounds.x1, bounds.y1)
    finally:
        plt.close(figure)


def test_cannot_silently_label_corrected_as_original(report):
    with pytest.raises(ValueError):
        validate(report)


def test_cannot_label_historical_as_corrected():
    path = Path(__file__).resolve().parents[2] / (
        "benchmark_tools/results/qfo_factorial_swiss_bootstrap_20260918.json")
    with pytest.raises(ValueError):
        validate(json.loads(path.read_text()), "corrected")


@pytest.mark.parametrize("change", ["status", "release", "protocol", "unknown", "missing", "effect"])
def test_corrected_identity_and_arithmetic(report, change):
    mode = "corrected"
    if change == "status":
        report["status"] = "paired_qfo_factorial_swiss_intervals"
    elif change == "release":
        report["input_release"] = "original"
    elif change == "protocol":
        report["corrected_protocol"]["sha256"] = "incorrect"
    elif change == "unknown":
        mode = "other"
    elif change == "missing":
        report["comparisons"].pop()
    else:
        report["comparisons"][0]["metrics"]["F1"]["difference"] += .1
    with pytest.raises(ValueError):
        validate(report, mode)
