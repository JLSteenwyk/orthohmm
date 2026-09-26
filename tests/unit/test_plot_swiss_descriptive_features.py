import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_swiss_descriptive_features import SOURCES, COMPLETE_SOURCES, validate, figure


def source(kind):
    root = Path(__file__).resolve().parents[2]
    with (root / "benchmark_tools/results" / SOURCES[kind][0]).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@pytest.mark.parametrize("kind,expected", [("identity", 48), ("fragment", 96)])
def test_actual_tables_and_plot_inventory(kind, expected):
    table = validate(source(kind), SOURCES[kind][2])
    fig, rows = figure(table, kind)
    try:
        fig.canvas.draw()
        assert len(rows) == expected
        assert all(r["difference"] is None for r in rows if r["method"] == "orthomcl_1_4")
        assert all(r["difference"] == 0 for r in rows if r["method"] == "orthofinder_3_1_5_full")
        assert len({r["method"] for r in rows}) == 8
    finally:
        plt.close(fig)


@pytest.mark.parametrize("problem", ["duplicate", "missing", "delta", "f1", "families", "nan", "fake_zero"])
def test_reject_changed_table(problem):
    rows = source("identity")
    if problem == "duplicate":
        rows.append(dict(rows[0]))
    elif problem == "missing":
        rows.pop()
    elif problem == "delta":
        rows[0]["delta_F1"] = "0"
    elif problem == "f1":
        rows[0]["F1"] = "0"
    elif problem == "families":
        rows[0]["families"] = "19"
    elif problem == "nan":
        rows[0]["PPV"] = "nan"
    elif problem == "fake_zero":
        next(r for r in rows if r["method"] == "orthomcl_1_4")["F1"] = "0"
    with pytest.raises(ValueError):
        validate(rows, SOURCES["identity"][2])


@pytest.mark.parametrize("kind,expected", [("identity", 48), ("fragment", 96)])
def test_complete_panel_requires_explicit_mode(kind, expected):
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    path, _, bins = COMPLETE_SOURCES[kind]
    with (root / path).open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with pytest.raises(ValueError):
        validate(rows, bins)
    with pytest.raises(ValueError):
        validate(source(kind), bins, complete=True)
    indexed = validate(rows, bins, complete=True)
    fig, endpoints = figure(indexed, kind, complete=True)
    try:
        fig.canvas.draw()
        assert len(endpoints) == expected
        assert all(r["difference"] is not None for r in endpoints)
        texts = [t.get_text() for ax in fig.axes for t in [*ax.texts, *ax.get_yticklabels()]]
        assert not any("unavailable" in t or "Not yet admitted" in t for t in texts)
    finally:
        plt.close(fig)


def test_complete_duplication_render(tmp_path):
    from benchmark_tools.plot_swiss_duplication_strata import render
    result = render(Path(__file__).resolve().parents[2], tmp_path / "complete", complete=True)
    assert result["endpoints"] == 48
    assert result["unavailable_endpoints"] == 0
    svg = (tmp_path / "complete/swiss_duplication_descriptive.svg").read_text()
    assert "Not yet admitted" not in svg
    assert "(unavailable)" not in svg
