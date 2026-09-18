import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools.plot_swiss_domain_strata import METRICS, POSITIONS, plot, rows, validate


def report():
    return json.loads((Path(__file__).resolve().parents[2] /
                       "benchmark_tools/results/swiss_domain_strata_results_20260917.json").read_text())


def test_all_27_points_and_54_intervals_match():
    data = report()
    fig = plot(data)
    try:
        assert len(fig.axes) == 3
        for ax, metric in zip(fig.axes, METRICS):
            assert len(ax.lines) == 28
            for index, row in enumerate(rows(data)):
                adjusted, nominal, point = ax.lines[1 + 3*index:4 + 3*index]
                assert adjusted.get_xdata() == pytest.approx(np.asarray(row[metric]["bonferroni27"]) * 100)
                assert nominal.get_xdata() == pytest.approx(np.asarray(row[metric]["nominal95"]) * 100)
                assert point.get_xdata()[0] == row[metric]["difference"] * 100
                assert point.get_ydata()[0] == POSITIONS[index]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        texts = list(fig.texts)
        for ax in fig.axes:
            texts.extend([ax.title, ax.xaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()])
        for text in texts:
            if text.get_text():
                box = text.get_window_extent(renderer)
                assert box.x0 >= 0 and box.y0 >= 0
                assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
        assert all(row["interaction_higher_minus_lower"][metric]["bonferroni27"][0] <= 0 <=
                   row["interaction_higher_minus_lower"][metric]["bonferroni27"][1]
                   for row in data["contrasts"] for metric in METRICS)
    finally:
        plt.close(fig)


@pytest.mark.parametrize("mutation", ["missing", "nan", "reverse", "outside", "order", "bins"])
def test_invalid_data_rejected(mutation):
    data = report()
    values = data["contrasts"][0]["primary_strata"][0]["metrics"]["F1"]
    if mutation == "missing":
        data["contrasts"].pop()
    elif mutation == "nan":
        values["difference"] = float("nan")
    elif mutation == "reverse":
        values["bonferroni27"].reverse()
    elif mutation == "outside":
        values["bonferroni27"][0] = -.9
    elif mutation == "order":
        data["contrasts"][0]["primary_strata"].reverse()
    else:
        data["strata"]["median_pfam_types_below_two"].pop()
    with pytest.raises(ValueError):
        validate(data)
