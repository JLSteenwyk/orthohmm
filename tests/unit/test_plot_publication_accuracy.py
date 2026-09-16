import json
from pathlib import Path

import pytest

from benchmark_tools.plot_publication_accuracy import plotting_data, paired_differences, plt


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/publication_comparison_20260916.json"
    return json.loads(path.read_text())


def test_plot_values_and_pending_mask(report):
    rows = plotting_data(report)
    assert len(rows) == 8
    assert rows[0]["orthobench"]["f_score"] == pytest.approx(70.35899765495243)
    assert rows[0]["qfo"]["VGNC F"]["x"] == report["methods"][0]["qfo"]["metric_details"]["VGNC F"]["participant"]["metric_x"]
    report["methods"][0]["qfo"]["status"] = "awaiting_workflow_completion"
    assert plotting_data(report)[0]["qfo"] == {}


def test_wrong_axes_rejected(report):
    report["methods"][0]["qfo"]["metric_details"]["VGNC F"]["axes"]["x_axis"] = "PPV"
    with pytest.raises(ValueError, match="axes"):
        plotting_data(report)


def test_missing_method_rejected(report):
    report["methods"].pop()
    with pytest.raises(ValueError, match="panel"):
        plotting_data(report)


def test_invalid_metric_rejected(report):
    report["methods"][0]["orthobench"]["recall_percent"] = float("nan")
    with pytest.raises(ValueError, match="metric"):
        plotting_data(report)


def test_inconsistent_paired_estimate_rejected(report):
    report["methods"][0]["orthobench"]["uncertainty"]["metrics"]["f_score"]["difference_percentage_points"] = 100
    with pytest.raises(ValueError, match="disagrees"):
        plotting_data(report)


def test_paired_plot_has_correct_baseline_and_interval_endpoints(report):
    rows = plotting_data(report)
    figure = paired_differences(rows)
    ax = figure.axes[0]
    first = rows[0]["uncertainty"]["f_score"]
    segment = ax.collections[0].get_segments()[0]
    assert segment[:, 0] == pytest.approx(first["bonferroni_percentile_ci"])
    assert ax.lines[-1].get_xdata() == [0, 0]
    plt.close(figure)
