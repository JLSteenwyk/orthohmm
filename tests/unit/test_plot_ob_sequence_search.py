import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_ob_sequence_search import plot, validate


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/ob_sequence_search_results_20260916.json"
    return json.loads(path.read_text())


def test_all_effects_shown_with_shared_axis_and_no_clipping(report):
    figure = plot(report)
    figure.canvas.draw()
    assert len(figure.axes) == 4
    for ax in figure.axes[1:]:
        assert ax.get_xlim() == (-40., 40.)
        assert len(ax.lines) == 7
    renderer = figure.canvas.get_renderer()
    for text in figure.texts:
        bbox = text.get_window_extent(renderer)
        assert bbox.x0 >= 0 and bbox.y0 >= 0
        assert bbox.x1 <= figure.bbox.width and bbox.y1 <= figure.bbox.height
    plt.close(figure)


@pytest.mark.parametrize("problem", ["baseline", "seed", "panel", "score", "coverage", "effect", "interval", "families", "multiplicity"])
def test_changed_evidence_rejected(report, problem):
    if problem == "baseline":
        report["baseline"] = "other"
    elif problem == "seed":
        report["seed"] = 1
    elif problem == "panel":
        report["comparisons"].pop("top100")
    elif problem == "score":
        report["point_estimates_percent"]["all_hits"]["f_score"] += 1
    elif problem == "coverage":
        report["gene_coverage"]["top100"]["assigned_genes"] -= 1
    elif problem == "effect":
        report["comparisons"]["top100"]["metrics"]["recall"]["difference_percentage_points"] += 1
    elif problem == "interval":
        report["comparisons"]["top100"]["metrics"]["recall"]["bonferroni_percentile_ci"] = [30, -30]
    elif problem == "families":
        report["families"].pop()
    else:
        report["multiplicity"] = "uncorrected"
    with pytest.raises(ValueError):
        validate(report)
