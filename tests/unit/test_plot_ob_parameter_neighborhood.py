from copy import deepcopy

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_ob_parameter_neighborhood as plotting
from benchmark_tools.parameter_neighborhood_statistics import summarize


def fixture():
    records = [{"refog": f"RefOG{i:03d}", "genes": 2, "true_positive": 1,
                "false_positive": 0, "false_negative": 0} for i in range(70)]
    scores = {name: deepcopy(records) for name in (plotting.BASELINE, *plotting.VARIANTS)}
    report = summarize(scores, {})
    report["scores"] = {name: {"refog_records": values} for name, values in scores.items()}
    report["coverage_resources"] = {name: {"input_genes": 251378, "assigned_genes": 251378} for name in scores}
    return report


def test_complete_figure():
    figure = plotting.plot(fixture())
    figure.canvas.draw()
    assert len(figure.axes) == 4
    assert len(figure.axes[0].tables[0].get_celld()) == 32
    for ax in figure.axes[1:]:
        assert len(ax.lines) == 19
    plt.close(figure)


@pytest.mark.parametrize("change", ["scores", "ci", "seed", "coverage", "missing", "failed", "published"])
def test_changed_evidence_rejected(change):
    report = fixture()
    first = plotting.VARIANTS[0]
    if change == "scores":
        report["point_estimates_percent"][first]["f_score"] = 99
    elif change == "ci":
        report["comparisons"][first]["metrics"]["f_score"]["paired_percentile_ci"] = [-1, 1]
    elif change == "seed":
        report["seed"] = 1
    elif change == "coverage":
        report["coverage_resources"][first]["assigned_genes"] = 1
    elif change == "missing":
        report["scores"].pop(first)
    elif change == "failed":
        report["failed_variants"] = {first: "failed"}
    else:
        report["publication_ready"] = True
    with pytest.raises(ValueError):
        plotting.validate(report)
