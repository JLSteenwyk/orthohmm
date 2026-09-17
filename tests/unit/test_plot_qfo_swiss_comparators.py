import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools.plot_qfo_swiss_comparators import METRICS, POSITIONS, plot, validate


def report():
    return json.loads((Path(__file__).resolve().parents[2] /
                       "benchmark_tools/results/qfo_swiss_comparator_intervals_20260917.json").read_text())


def test_all_points_and_intervals_match_source():
    data = report()
    fig = plot(data)
    try:
        assert len(fig.axes) == 3
        for ax, metric in zip(fig.axes, METRICS):
            assert len(ax.lines) == 26
            for index, row in enumerate(data["comparisons"]):
                values = row["metrics"][metric]
                adjusted, nominal, point = ax.lines[2 + 3*index:5 + 3*index]
                assert adjusted.get_xdata() == pytest.approx(np.asarray(values["bonferroni_percentile_ci"]) * 100)
                assert nominal.get_xdata() == pytest.approx(np.asarray(values["paired_percentile_ci"]) * 100)
                assert point.get_xdata()[0] == values["difference"] * 100
                assert point.get_ydata()[0] == POSITIONS[index]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        texts = list(fig.texts)
        for ax in fig.axes:
            texts.extend([ax.title, ax.xaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()])
        for text in texts:
            if not text.get_text():
                continue
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0
            assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
    finally:
        plt.close(fig)


@pytest.mark.parametrize("mutation", ["missing", "nan", "reverse", "outside"])
def test_changed_or_unplottable_data_rejected(mutation):
    data = report()
    metric = data["comparisons"][0]["metrics"]["F1"]
    if mutation == "missing":
        data["comparisons"].pop()
    elif mutation == "nan":
        metric["difference"] = float("nan")
    elif mutation == "reverse":
        metric["bonferroni_percentile_ci"].reverse()
    else:
        metric["bonferroni_percentile_ci"][0] = -.9
    with pytest.raises(ValueError):
        validate(data)
