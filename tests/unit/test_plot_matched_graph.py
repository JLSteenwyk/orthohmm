import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

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
