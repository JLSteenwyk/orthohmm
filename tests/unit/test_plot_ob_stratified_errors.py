import copy
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_ob_stratified_errors as plotter


@pytest.fixture
def report():
    path = Path(plotter.__file__).parent / "results/ob_stratified_error_results_20260916.json"
    return json.loads(path.read_text())


def test_complete_endpoint_inventory(report):
    assert plotter.validate(report) == {"planned_endpoints": 84, "finite_points": 72,
                                       "interval_endpoints": 66, "nonestimable_endpoints": 12}


@pytest.mark.parametrize("change", ["missing", "multiplicity", "families", "point", "interval", "nonfinite", "small_ci", "empty_zero", "bootstrap"])
def test_changed_evidence_rejected(report, change):
    row = report["strata"]["size:large_gt_50"]
    value = row["comparisons"][plotter.COMPARATORS[0]]["metrics"]["f_score"]
    if change == "missing":
        report["strata"].pop("composition:missing")
    elif change == "multiplicity":
        report["multiplicity_endpoints"] = 66
    elif change == "families":
        row["families"] = row["families"][:-1]
    elif change == "point":
        row["point_estimates_percent"][plotter.BASELINE]["f_score"] += 1
    elif change == "interval":
        value["bonferroni_percentile_ci"] = [0., 0.]
        row["bootstrap"]["comparisons"] = copy.deepcopy(row["comparisons"])
    elif change == "nonfinite":
        value["difference_percentage_points"] = float("nan")
        row["bootstrap"]["comparisons"] = copy.deepcopy(row["comparisons"])
    elif change == "small_ci":
        report["strata"]["composition:concentrated"]["comparisons"][plotter.COMPARATORS[0]]["metrics"]["f_score"]["paired_percentile_ci"] = [-1, 1]
    elif change == "empty_zero":
        report["strata"]["composition:missing"]["comparisons"][plotter.COMPARATORS[0]]["metrics"]["f_score"]["difference_percentage_points"] = 0.
    elif change == "bootstrap":
        row["bootstrap"]["seed"] = 0
    with pytest.raises(ValueError):
        plotter.validate(report)


def test_all_strata_points_intervals_and_missing_rows_rendered(report):
    fig = plotter.plot(report)
    try:
        fig.canvas.draw()
        assert len(fig.axes) == 3
        for ax in fig.axes:
            points = [line for line in ax.lines if line.get_marker() in ("o", "s")]
            assert len(points) == 24
            assert sum(line.get_markerfacecolor() == "white" for line in points) == 2
            assert len(ax.lines) == 24 + 44 + 5
            assert [text.get_text() for text in ax.texts] == ["No families", "No families"]
            assert ax.get_xlim() == fig.axes[0].get_xlim()
        assert len(fig.axes[0].get_yticklabels()) == 14
        assert "0/22" in " ".join(text.get_text() for text in fig.texts)
        renderer = fig.canvas.get_renderer()
        for text in [*fig.texts, *fig.axes[0].get_yticklabels()]:
            bounds = text.get_window_extent(renderer)
            assert bounds.x0 >= 0 and bounds.y0 >= 0
            assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
    finally:
        plt.close(fig)
