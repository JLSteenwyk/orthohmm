import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_ygob_validation import plot


@pytest.fixture
def report():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/ygob_frozen_results_20260916.json").read_text())


def test_all_methods_metrics_and_intervals_render(report):
    fig = plot(report)
    left, right = fig.axes
    assert len(left.patches) == 12
    assert len(right.lines) == 19  # Zero line plus six adjusted/nominal/point triplets.
    assert left.get_xlim() == (0, 100)
    assert len(right.get_yticklabels()) == 6
    plt.close(fig)


def test_unverified_results_rejected(report):
    report["frozen_evaluation_gates_verified"] = False
    with pytest.raises(ValueError, match="admitted"):
        plot(report)


def test_wrong_multiplicity_rejected(report):
    report["uncertainty"]["multiplicity_count"] = 1
    with pytest.raises(ValueError, match="uncertainty"):
        plot(report)


def test_nonfinite_score_rejected(report):
    report["scores"]["orthofinder_full"]["metrics"]["f1"] = float("nan")
    with pytest.raises(ValueError, match="score"):
        plot(report)
