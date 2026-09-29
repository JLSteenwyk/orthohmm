from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.plot_ob_complete_strata import matrices, plot


def report():
    return json.loads((Path(__file__).resolve().parents[2] /
        "benchmark_tools/results/ob_complete_strata_20260928/report.json").read_text())


def test_complete_endpoint_inventory():
    values, counts = matrices(report())
    assert values.shape == (3, 14, 8)
    assert np.isfinite(values).sum() == 288
    assert np.isnan(values).sum() == 48
    assert sum(counts) == 350


def test_rendered_cells_and_text_bounds():
    import matplotlib.pyplot as plt
    data = report()
    values, _ = matrices(data)
    figure = plot(data)
    try:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        for k, axis in enumerate(figure.axes[:3]):
            expected = ["NA" if np.isnan(v) else f"{v:.1f}" for v in values[k].flat]
            assert [t.get_text() for t in axis.texts] == expected
            assert np.allclose(axis.images[0].get_array().filled(np.nan), values[k], equal_nan=True)
            for text in [*axis.texts, *axis.get_xticklabels(), *axis.get_yticklabels()]:
                if text.get_visible():
                    box = text.get_window_extent(renderer)
                    assert box.x0 >= 0 and box.y0 >= 0
                    assert box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
        for text in figure.texts:
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0
            assert box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
    finally:
        plt.close(figure)


@pytest.mark.parametrize("change", ["duplicate", "missing", "zero_empty", "nonfinite", "membership"])
def test_reject_invalid_display(change):
    data = deepcopy(report())
    rows = data["rows"]
    if change == "duplicate":
        rows[-1] = deepcopy(rows[0])
    elif change == "missing":
        rows.pop()
    elif change == "zero_empty":
        next(r for r in rows if not r["families"])["metrics_percent"] = dict(f_score=0, precision=0, recall=0)
    elif change == "nonfinite":
        next(r for r in rows if r["families"])["metrics_percent"]["f_score"] = float("nan")
    else:
        next(r for r in rows if r["families"])["families"] = ["wrong"]
    with pytest.raises(ValueError):
        matrices(data)
