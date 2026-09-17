import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_simulation_tree_robustness as plotter
from benchmark_tools.summarize_simulation_tree_panel import SEEDS, ARMS


def report(failed=False):
    score = dict(tp=5, fp=5, fn=5, input_genes=100, eligible_true_pairs=10,
                 f1=.5, precision=.5, recall=.5, undefined_ratios=[])
    rows = [dict(condition=c, seed=s, method=m, arm=a, truth_sha256="a" * 64,
                 **(dict(status="failed", reason="fixture") if failed else dict(status="complete", score=score.copy())))
            for c in plotter.CONDITIONS for s in SEEDS for m in plotter.METHODS for a in ARMS]
    return plotter.summarize(rows)


@pytest.mark.parametrize("problem", ["status", "missing", "point", "interval", "bootstrap", "record"])
def test_modified_evidence_rejected(problem):
    data = report()
    if problem == "status":
        data["status"] = "partial"
    elif problem == "missing":
        data["contrasts"].pop()
    elif problem == "point":
        data["contrasts"][0]["metrics"]["f1"]["difference_percentage_points"] = 1.
    elif problem == "interval":
        data["contrasts"][0]["metrics"]["f1"]["bonferroni_126_ci"] = [-1., 1.]
    elif problem == "bootstrap":
        data["bootstrap"]["multiplicity"] = 3
    else:
        data["records"].pop()
    with pytest.raises(ValueError):
        plotter.validate(data)


def test_all_126_endpoints_rendered():
    fig = plotter.plot(report())
    try:
        fig.canvas.draw()
        assert len(fig.axes) == 9
        assert all(len(ax.lines) == 43 for ax in fig.axes)
        assert "126 exploratory endpoints" in " ".join(t.get_text() for t in fig.texts)
    finally:
        plt.close(fig)


def test_unavailable_endpoints_never_draw_zero_effect():
    fig = plotter.plot(report(failed=True))
    try:
        fig.canvas.draw()
        assert all(len(ax.lines) == 1 and len(ax.texts) == 14 for ax in fig.axes)
        assert "560 failures retained" in " ".join(t.get_text() for t in fig.texts)
    finally:
        plt.close(fig)


def test_single_pair_effect_is_visible_without_interval():
    data = report()
    for row in data["records"]:
        if (row["condition"], row["method"], row["arm"]) == (plotter.CONDITIONS[0], plotter.METHODS[0], "generating"):
            if row["seed"] == SEEDS[0]:
                row["score"].update(tp=9, fp=1, fn=1, f1=.9, precision=.9, recall=.9)
            else:
                row.update(status="failed", reason="fixture")
                del row["score"]
    fig = plotter.plot(plotter.summarize(data["records"]))
    try:
        assert fig.axes[0].get_xlim()[1] > 40.
        assert len(fig.axes[0].lines) == 41
    finally:
        plt.close(fig)
