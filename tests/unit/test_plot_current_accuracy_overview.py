import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest
import matplotlib.pyplot as plt

from benchmark_tools.plot_current_accuracy_overview import SOURCES, export, plot, validate

BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


@pytest.fixture(scope="module")
def reports():
    return {key: json.loads((BASE / name).read_bytes()) for key, (name, _) in SOURCES.items()}


def test_current_points_preserve_all_methods_and_select_matched_sonic(reports):
    data, rows, changes = validate(reports)
    assert len(data) == 8 and len(rows) == 24
    assert len({(r["key"], r["benchmark"], r["metric"]) for r in rows}) == 24
    assert {r["units"] for r in rows} == {"percent"}
    sonic = next(r for r in data if r["key"] == "sonicparanoid_2_0_9")
    assert sonic["three_kingdoms_f1"] == 0.9912758996728462
    assert changes[0]["historical_f1"] == 0.9907944084555063
    assert changes[0]["difference_pp"] == pytest.approx(0.04814912173399799)
    assert changes[0]["current_input_status"] == "matched input verified"
    assert len(changes) == 1
    table = {r["key"]: r for r in reports["table"]["rows"]}
    for point in data:
        assert point["three_kingdoms_f1"] == table[point["key"]]["scores"]["ThreeKingdoms"]
        assert point["orthobench"] == {m: table[point["key"]]["orthobench_retained_evidence"][m + "_percent"]
                                      for m in ("precision", "recall", "f_score")}


@pytest.mark.parametrize("change", ["status", "units", "count", "inventory", "duplicate", "qfo_semantics",
    "qfo_axis", "qfo_coordinate", "history", "tk_semantics", "tk_input", "mean", "nonfinite"])
def test_inconsistent_retained_evidence_is_not_promoted(reports, change):
    changed = copy.deepcopy(reports)
    row = changed["table"]["rows"][0]
    if change == "status":
        changed["table"]["status"] = "historical"
    elif change == "units":
        row["scores"]["OrthoBench"] *= 100
    elif change == "count":
        changed["kingdoms"]["rows"][0]["counts"]["true_positive_gene_pairs"] += 1
    elif change == "inventory":
        changed["register"]["rows"].pop()
    elif change == "duplicate":
        changed["table"]["rows"][1] = row
    elif change == "qfo_semantics":
        row["qfo_prediction_semantics"] = "pre-clustering graph"
    elif change == "qfo_axis":
        changed["qfo"]["methods"][0]["details"]["GO"]["statistic"] = "F1"
    elif change == "qfo_coordinate":
        changed["qfo_figure"]["plotted_data"][0]["qfo"]["GO"]["x"] += 1
    elif change == "history":
        old = next(r for r in changed["historical"]["plotted_data"] if r["key"] == "sonicparanoid_2_0_9")
        old["three_kingdoms_f1"] = 0.9912758996728462
    elif change == "tk_semantics":
        changed["kingdoms"]["rows"][0]["semantics"] = "genome-wide orthology"
    elif change == "tk_input":
        row["three_kingdoms"]["input_status"] = "verified consumption"
    elif change == "mean":
        row["scores"]["QfO_secondary_mean"] = 0
    else:
        row["scores"]["GO"] = float("nan")
    with pytest.raises(ValueError):
        validate(changed)


def test_changed_source_hash_fails_without_output(tmp_path):
    name = SOURCES["table"][0]
    source = tmp_path / "inputs" / name
    source.parent.mkdir(parents=True)
    source.write_bytes((BASE / name).read_bytes() + b"\n")
    output = tmp_path / "figure"
    with pytest.raises(ValueError, match="Changed overview source"):
        export(tmp_path / "inputs", output)
    assert not output.exists()


def test_caption_boxes_avoid_axis_labels_and_legend(reports):
    figure = plot(validate(reports)[0])
    try:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        captions = [t.get_window_extent(renderer) for t in figure.texts[1:]]
        labels = [ax.xaxis.label.get_window_extent(renderer) for ax in figure.axes]
        legend = figure.legends[0].get_window_extent(renderer)
        assert not captions[0].overlaps(captions[1])
        for box in captions:
            assert not any(box.overlaps(label) for label in labels)
            assert not box.overlaps(legend)
            assert figure.bbox.contains(*box.get_points()[0])
            assert figure.bbox.contains(*box.get_points()[1])
    finally:
        plt.close(figure)


def test_actual_export_preserves_source_scores_and_historical_figure(tmp_path):
    historical = BASE / "figures_accuracy_orthomcl_complete_20260916/manifest.json"
    before = historical.read_bytes()
    output = tmp_path / "new"
    result = export(BASE, output)
    assert result["cross_dataset_score_cells_verified"] == 72
    assert result["corrected_qfo_points_verified"] == 48
    assert result["overview_plotted_values"] == 24
    for key in ("native_predictions_or_scores_changed", "uncertainty_recomputed", "native_inference_or_scoring_rerun",
                "historical_input_consumption_proven", "controlled_timing", "publication_ready",
                "existing_corrected_qfo_figure_rerendered", "figure_visually_reviewed"):
        assert result[key] is False
    with (output / "plotted_values.tsv").open() as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 24
    for pin in result["outputs"]:
        path = Path(pin["path"])
        assert path.stat().st_size == pin["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == pin["sha256"]
    assert historical.read_bytes() == before
    with pytest.raises(FileExistsError):
        export(BASE, output)


def test_prose_distinguishes_current_points_from_historical_figures():
    overview = " ".join((BASE / "CURRENT_ACCURACY_OVERVIEW_20261002.md").read_text().split())
    extended = " ".join((BASE / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text().split())
    assert "not the primary corrected-release results" in extended
    assert "On the original-input historical release" in extended
    assert "not a newly matched eight-tool experiment" in extended
    assert "0.9912758997" in overview and "+0.048149 percentage points" in overview
    assert "GO/EC/FAS are similarities, not F1" in overview
    assert "No new manuscript render or archive is generated" in overview
    assert "other historical consumption gaps remain" in overview
    for name in ("figures_current_accuracy_overview_20261002_v2/current_accuracy_overview.png",
                 "CURRENT_ACCURACY_OVERVIEW_20261002.md", "current_benchmark_scores_20260926_v2/scores.tsv"):
        assert name in extended and (BASE / name).is_file()
