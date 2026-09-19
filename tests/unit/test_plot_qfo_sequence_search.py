import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from benchmark_tools.plot_qfo_sequence_search import plot, validate, render


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_sequence_swiss_bootstrap_20260918.json"
    return json.loads(path.read_text())


def test_full_panel_and_percent_conversion(report):
    fig = plot(report)
    fig.canvas.draw()
    assert len(fig.axes) == 4
    for i, ax in enumerate(fig.axes[1:]):
        assert ax.get_xlim() == (-30, 30)
        assert len(ax.lines) == 7
        metric = ("F1", "PPV", "TPR")[i]
        for j, contrast in enumerate(report["comparisons"]):
            np.testing.assert_allclose(ax.lines[1 + j * 3].get_xdata(),
                100 * np.array(contrast["metrics"][metric]["bonferroni_percentile_ci"]))
    renderer = fig.canvas.get_renderer()
    for text in fig.texts:
        box = text.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0
        assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
    plt.close(fig)


@pytest.mark.parametrize("problem", ["status", "seed", "units", "families", "panel", "orientation",
    "point", "effect", "interval", "nan", "wins", "clip", "adjustment"])
def test_changed_evidence_rejected(report, problem):
    contrast = report["comparisons"][0]
    row = contrast["metrics"]["F1"]
    if problem == "status":
        report["uncertainty_admitted"] = False
    elif problem == "seed":
        report["seed"] += 1
    elif problem == "units":
        report["units"] = "percent"
    elif problem == "families":
        report["families"].pop()
    elif problem == "panel":
        report["comparisons"].reverse()
    elif problem == "orientation":
        contrast["reference"] = "top100"
    elif problem == "point":
        report["point_estimates"]["all_hits"]["F1"] += .01
    elif problem == "effect":
        row["difference"] += .01
    elif problem == "interval":
        row["bonferroni_percentile_ci"] = [1, -1]
    elif problem == "nan":
        row["paired_percentile_ci"][0] = float("nan")
    elif problem == "wins":
        row["family_wins"] += 1
    elif problem == "clip":
        row["bonferroni_percentile_ci"] = [-.5, .5]
    else:
        report["multiplicity_endpoints"] = 3
    with pytest.raises(ValueError):
        validate(report)


def test_render_refuses_changed_source_and_existing_output(tmp_path):
    results = tmp_path / "result.json"
    results.write_text("{}")
    output = tmp_path / "figure"
    with pytest.raises(ValueError):
        render(results, output)
    assert not output.exists()
    output.symlink_to(tmp_path / "absent")
    with pytest.raises(FileExistsError):
        render(results, output)
