import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_qfo_hit_coverage import plot


def fixture():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_hit_coverage_20260918/summary.json"
    return json.loads(path.read_text())


def test_bars_match_source_counts_and_zero_based_axes():
    report = fixture()
    fig = plot(report)
    try:
        counts, coverage = fig.axes
        assert [p.get_width() for p in counts.patches] == pytest.approx(
            [r["directed_hits"] / 1e6 for r in report["rows"]])
        assert [p.get_width() for p in coverage.patches] == pytest.approx(
            [100 * (r["genes"] - r["queries_without_cross_species_hits"]) / r["genes"]
             for r in report["rows"]])
        assert counts.get_xlim()[0] == 0
        assert coverage.get_xlim() == (0, 100)
        fig.canvas.draw()
    finally:
        plt.close(fig)


@pytest.mark.parametrize("problem", ["accuracy", "order", "genes", "count"])
def test_changed_identity_or_invalid_counts_rejected(problem):
    report = fixture()
    if problem == "accuracy":
        report["accuracy_evaluated"] = True
    elif problem == "order":
        report["rows"].reverse()
    elif problem == "genes":
        report["rows"][0]["genes"] -= 1
    else:
        report["rows"][0]["queries_without_cross_species_hits"] = 984138
    with pytest.raises(ValueError):
        plot(report)
