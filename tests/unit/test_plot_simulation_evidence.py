import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_simulation_evidence import plot, CONDITIONS, METHODS


def report(estimated):
    return {"panel": "variable_length_v2" if estimated else "fixed_length_v1", "conditions": {
        c: {"methods": {m: {"complete_seeds": list(range(10 if j < 2 else 0))} for j, m in enumerate(METHODS)},
            "contrasts": {m: {"included_seeds": [1, 2] if estimated else [],
                              "metrics": {"f1": {"difference_percentage_points": -2,
                                                  "paired_95_percent_ci": [-3, -1],
                                                  "bonferroni_14_ci": [-4, 0]}} if estimated else {}}
                          for m in METHODS[:2]}} for c in CONDITIONS}}


def test_intervals_are_paired_values_not_differences_of_table_means():
    fig = plot(report(True))
    try:
        ax = fig.axes[1]
        assert len(ax.lines) == 43
        assert list(ax.lines[1].get_xdata()) == [-4, 0]
        assert list(ax.lines[2].get_xdata()) == [-3, -1]
        assert list(ax.lines[3].get_xdata()) == [-2]
        assert all(t.get_text() == "n=2" for t in ax.get_yticklabels())
    finally:
        plt.close(fig)


def test_no_pairs_do_not_plot_zero_effect_estimates():
    fig = plot(report(False))
    try:
        assert len(fig.axes[1].lines) == 1
        assert "No admitted paired comparisons" in fig.axes[1].texts[0].get_text()
        assert fig.axes[0].images[0].get_array()[0].tolist() == [10, 10, 0, 0]
    finally:
        plt.close(fig)


def test_unknown_panel_rejected():
    with pytest.raises(ValueError):
        plot({"panel": "unknown"})
