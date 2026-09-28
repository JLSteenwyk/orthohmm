import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_ygob_overlap_strata as plotter


@pytest.fixture
def report():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/ygob_overlap_strata_20260928.json").read_text())


def test_all_values_recomputed_without_dropping_methods(report):
    rows = plotter.plotting_data(report)
    assert len(rows) == 24
    for row in rows:
        expected = report["strata"][row["stratum"]][row["method"]]
        assert row["score_percent"] == pytest.approx(100 * expected["metrics"][row["metric"]])
        assert row["fp"] == expected["counts"]["fp"]


@pytest.mark.parametrize("field", ["independent_confirmation", "confidence_intervals_computed"])
def test_new_inference_claims_rejected(report, field):
    report[field] = True
    with pytest.raises(ValueError):
        plotter.plotting_data(report)


@pytest.mark.parametrize("mutation", ["method", "stratum", "metric", "defined", "fp", "tp", "sizes"])
def test_inconsistent_records_rejected(report, mutation):
    item = report["strata"]["screen_negative"]["orthohmm_satellite_v2"]
    if mutation == "method":
        del report["strata"]["screen_positive"]["orthofinder_full"]
    elif mutation == "stratum":
        report["strata"]["extra"] = {}
    elif mutation == "metric":
        item["metrics"]["f1"] += .001
    elif mutation == "defined":
        item["defined"]["recall"] = False
    elif mutation == "fp":
        item["counts"]["fp"] = .3
    elif mutation == "tp":
        item["counts"]["tp"] = True
    elif mutation == "sizes":
        item["reference_genes"] += 1
    with pytest.raises(ValueError):
        plotter.plotting_data(report)


def test_six_panels_show_all_twenty_four_points(report):
    fig = plotter.plot(plotter.plotting_data(report))
    try:
        assert len(fig.axes) == 6
        assert all(ax.get_xlim() == (0, 100) for ax in fig.axes)
        assert sum(len(ax.collections) for ax in fig.axes) == 24
        assert all(len(ax.texts) == 4 for ax in fig.axes)
        notes = " ".join(t.get_text() for t in fig.texts)
        assert "2,893 singleton" in notes and "2,017 singleton" in notes
        assert "does not establish family independence" in notes
    finally:
        plt.close(fig)


def test_undefined_is_not_drawn_as_zero(report):
    for item in report["strata"]["screen_negative"].values():
        item["counts"] = dict(tp=0, fp=0, fn=0)
        item["truth_pairs"] = 0
        item["metrics"] = dict.fromkeys(plotter.METRICS, 0.0)
        item["defined"] = dict.fromkeys(plotter.METRICS, False)
    fig = plotter.plot(plotter.plotting_data(report))
    try:
        assert sum(len(ax.collections) for ax in fig.axes) == 12
        assert all(t.get_text() == "undefined" for ax in fig.axes[3:] for t in ax.texts)
    finally:
        plt.close(fig)


def test_render_manifest_and_refuse_overwrite(report, tmp_path):
    source = tmp_path / "source.json"
    source.write_text(json.dumps(report))
    sha = plotter.file_provenance(source)["sha256"]
    output = tmp_path / "figure"
    manifest = plotter.render(source, sha, output)
    assert manifest["rows"] == 24 and len(manifest["outputs"]) == 4
    for ref in manifest["outputs"]:
        assert plotter.file_provenance(Path(ref["path"])) == ref
    with (output / "plotted_values.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 24
    assert {row["label"] for row in rows} == set(plotter.NAMES)
    with pytest.raises(FileExistsError):
        plotter.render(source, sha, output)
