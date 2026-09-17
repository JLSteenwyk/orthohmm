import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_wgd_application import ENDPOINTS, METHODS, cases, plot


def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/biological_wgd_results_20260917.json"
    return json.loads(path.read_text())


def test_plot_uses_fixed_population_means_and_all_intervals():
    data = report()
    fig = plot(data)
    try:
        assert len(fig.axes) == 6
        for index, endpoint in enumerate(ENDPOINTS):
            for line, method in zip(fig.axes[index].lines, METHODS):
                assert line.get_xdata()[0] == data["methods"][method]["summary"]["endpoints"][endpoint]["mean"] * 100
            bottom = fig.axes[index + 3]
            assert len(bottom.lines) == 13  # Zero line and three marks for each contrast.
            contrasts = [c for c in data["uncertainty"]["comparisons"] if c["endpoint"] == endpoint]
            for j, contrast in enumerate(contrasts):
                assert list(bottom.lines[1 + 3 * j].get_xdata()) == contrast["bonferroni12_pp"]
                assert bottom.lines[3 + 3 * j].get_xdata()[0] == contrast["difference_pp"]
    finally:
        plt.close(fig)


def test_cases_preserve_all_six_and_reference_exclusion():
    data = report()
    text = cases(data)
    for example in data["prespecified_examples"]:
        assert " / ".join(example["orf_pair"]) in text
    assert text.count("| OrthoFinder MCL checkpoint (diagnostic) |") == 6
    assert "anchors_in_different_pillars" in text
    assert "| OrthoHMM phylogeny | separated | 0, 3 | 3/6 |" in text


def test_changed_population_fails_before_plotting():
    data = report()
    data["cohort_pairs"] = 239
    with pytest.raises(ValueError, match="Unexpected frozen population"):
        plot(data)
