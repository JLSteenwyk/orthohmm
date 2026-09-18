import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_dgx_descriptive_resources as module


def report():
    return json.loads((Path(__file__).resolve().parents[2] /
        "benchmark_tools/results/dgx_descriptive_resources_20260918.json").read_text())


def test_actual_panel_shows_all_repeats_without_fitted_lines():
    fig = module.plot(report())
    assert len(fig.axes) == 3
    for ax in fig.axes:
        assert len(ax.collections) == 9
        assert sum(len(c.get_offsets()) for c in ax.collections) == 27
        assert len(ax.lines) == 18
        assert ax.get_ylim()[0] == 0
    text = " ".join(t.get_text() for t in fig.texts)
    assert "not confidence intervals" in text
    assert "All host assessments remain inconclusive" in text
    assert "Sampled aggregate RSS is not plotted" in text
    plt.close(fig)


@pytest.mark.parametrize("problem", ["admitted", "missing", "failed", "host", "repeat", "index", "score", "summary", "genes"])
def test_changed_evidence_refused(problem):
    data = report()
    if problem == "admitted":
        data["scientific_timings_admitted"] = 27
    elif problem == "missing":
        data["runs"].pop()
    elif problem == "failed":
        data["runs"][0]["native_validation_status"] = "native_validation_failed"
    elif problem == "host":
        data["runs"][0]["host_status"] = "quiet"
    elif problem == "repeat":
        data["runs"][0]["repeat"] = 1
    elif problem == "index":
        data["runs"][0]["index"] = 1
    elif problem == "score":
        data["runs"][0]["native_elapsed_seconds"] = float("nan")
    elif problem == "summary":
        data["summaries"][0]["native_elapsed_seconds"]["median"] = 0
    else:
        data["runs"][0]["proteins"] = 1
    with pytest.raises(ValueError):
        module.validate(data)
