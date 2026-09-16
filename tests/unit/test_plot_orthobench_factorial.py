import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_orthobench_factorial import plot, validate


@pytest.fixture
def report():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/orthobench_factorial_results_20260916.json").read_text())


def test_all_cells_and_all_endpoints_render(report):
    fig = plot(report)
    assert fig.axes[0].images[0].get_array().shape == (8, 3)
    assert len(fig.axes[0].texts) == 24
    for ax in fig.axes[1:]:
        assert len(ax.lines) == 39  # 12 interval/point triplets, zero and two separators.
        assert len(ax.get_yticks()) == 12
    plt.close(fig)


@pytest.mark.parametrize("change", ["multiplicity", "missing", "order", "interval", "point", "nan"])
def test_rejects_changed_or_inconsistent_results(report, change):
    if change == "multiplicity":
        report["multiplicity_endpoints"] = 12
    elif change == "missing":
        report["comparisons"].pop()
    elif change == "order":
        report["comparisons"].reverse()
    elif change == "interval":
        report["comparisons"][0]["metrics"]["f_score"]["bonferroni_percentile_ci"] = [0, 0]
    elif change == "point":
        report["comparisons"][0]["metrics"]["f_score"]["difference_percentage_points"] += 1
    else:
        report["point_estimates_percent"]["p0_c0_r0"]["precision"] = float("nan")
    with pytest.raises(ValueError):
        validate(report)
