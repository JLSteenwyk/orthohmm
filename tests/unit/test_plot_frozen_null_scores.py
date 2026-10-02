import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.plot_frozen_null_scores import export, load, plot, validate

BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
SCORES = BASE / "frozen_null_score_observations_20261002.json.gz"
RECEIPT = BASE / "frozen_null_score_receipt_20261002.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def report():
    return json.loads(gzip.decompress(SCORES.read_bytes()))


def test_every_planned_endpoint_and_zero_hit_glyph_is_preserved(report):
    rows = validate(report)
    assert len(rows) == 90
    assert len({(r["regime"], r["length"], r["band"], r["threshold"]) for r in rows}) == 90
    assert {r["trials"] for r in rows} == {10000}
    zeros = [r for r in rows if r["zero_hit_upper_limit"]]
    assert zeros and all(r["fraction"] == r["adjusted_low"] == 0 < r["adjusted_high"] for r in zeros)
    assert all(r["adjusted_low"] <= r["fraction"] <= r["adjusted_high"] for r in rows)


@pytest.mark.parametrize("field", ["hits", "fraction", "bonferroni_clopper_pearson"])
def test_changed_retained_statistic_is_rejected(report, field):
    changed = copy.deepcopy(report)
    endpoint = changed["summaries"][0]["bands"]["0"][0]
    endpoint[field] = [0, 1] if field == "bonferroni_clopper_pearson" else endpoint[field] + 1
    with pytest.raises(ValueError, match="disagrees"):
        validate(changed)


@pytest.mark.parametrize("change", ["duplicate", "missing", "negative", "calibration", "paired"])
def test_changed_inventory_or_scope_is_rejected(report, change):
    changed = copy.deepcopy(report)
    if change == "duplicate":
        changed["rows"][1] = changed["rows"][0]
    elif change == "missing":
        changed["summaries"].pop()
    elif change == "negative":
        changed["rows"][0]["scores"]["0"][0] = -1
    elif change == "calibration":
        changed["calibration_established"] = True
    else:
        changed["summaries"][0]["band_changed_scores"] += 1
    with pytest.raises(ValueError):
        validate(changed)


def test_external_hash_gates_are_not_replaced_by_local_metadata(tmp_path):
    with pytest.raises(ValueError):
        load(SCORES, "0" * 64, RECEIPT, sha(RECEIPT))
    with pytest.raises(ValueError):
        load(SCORES, sha(SCORES), RECEIPT, "0" * 64)
    changed = tmp_path / "changed.gz"
    changed.write_bytes(SCORES.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="pinned receipt"):
        load(changed, sha(changed), RECEIPT, sha(RECEIPT))


def test_all_panels_and_glyphs_render_within_figure(report):
    rows = validate(report)
    figure = plot(rows)
    try:
        figure.canvas.draw()
        assert len(figure.axes) == 9
        assert sum(len(ax.lines) for ax in figure.axes) == 99
        zero_glyphs = [line for ax in figure.axes for line in ax.lines if line.get_marker() == "v"]
        assert len(zero_glyphs) == sum(r["zero_hit_upper_limit"] for r in rows)
        renderer = figure.canvas.get_renderer()
        for ax in figure.axes:
            label = ax.texts[0].get_window_extent(renderer)
            for line in ax.lines:
                if line.get_marker() in ("o", "s", "v"):
                    assert not label.overlaps(line.get_window_extent(renderer))
        for text in [*figure.texts, *[ax.title for ax in figure.axes]]:
            assert figure.bbox.contains(*text.get_window_extent(renderer).get_points()[0])
            assert figure.bbox.contains(*text.get_window_extent(renderer).get_points()[1])
    finally:
        plt.close(figure)


def test_actual_export_recounts_all_rows_and_preserves_attempt(tmp_path):
    output = tmp_path / "figure"
    result = export(SCORES, sha(SCORES), RECEIPT, sha(RECEIPT), output)
    assert result["endpoints"] == 90 and result["independent_pairs"] == 90000
    assert result["source_score_counts_recomputed"] is True
    for key in ("native_scoring_rerun", "calibration_established", "publication_ready", "figure_visually_reviewed"):
        assert result[key] is False
    with (output / "plotted_values.tsv").open() as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 90
    for pin in result["outputs"]:
        path = Path(pin["path"])
        assert path.stat().st_size == pin["bytes"] and sha(path) == pin["sha256"]
    before = (output / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError):
        export(SCORES, sha(SCORES), RECEIPT, sha(RECEIPT), output)
    assert (output / "manifest.json").read_bytes() == before
