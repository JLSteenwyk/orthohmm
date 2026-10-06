"""Guard complete native diagnostic cohort joins and exact figure data."""

import copy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import plot_native_qfo_swiss_mechanism as plotter

RESULTS = Path(plotter.__file__).parent / "results"
FILES = {
    "search": ("native_qfo_swiss_search_support_20261006_v1.json", "native_qfo_swiss_search_support_code_readback_20261006.json"),
    "reconciliation": ("native_qfo_swiss_reconciliation_trace_20261006_v1.json", "native_qfo_swiss_reconciliation_newick_readback_20261006.json"),
    "strata": ("native_qfo_swiss_sequence_strata_20261006_v1/report.json", "native_qfo_swiss_sequence_strata_readback_20261006.json"),
}


@pytest.fixture
def reports():
    return [json.loads((RESULTS / FILES[k][0]).read_text()) for k in ("search", "reconciliation", "strata")]


def test_actual_figure_data_matches_complete_reports(reports):
    stage, scores, differences = plotter.figure_data(*reports)
    assert len(stage) == len(scores) == len(differences) == 10
    assert [(r["count"], r["total"]) for r in stage if r["panel"] == "A"] == [
        (60, 334), (6, 334), (268, 334), (479, 1689), (20, 1689), (1190, 1689)]
    assert [(r["count"], r["total"]) for r in stage if r["panel"] == "B"] == [
        (280, 334), (54, 334), (1689, 1689), (0, 1689)]
    for panel in ("A", "B"):
        for label in ("TP", "FP"):
            assert sum(r["percentage"] for r in stage if r["panel"] == panel and r["label"] == label) == pytest.approx(100)
    original = {(r["cell"], r["stratum"]): r for r in reports[2]["rows"]}
    assert all(all(r[k] == original[r["cell"], r["stratum"]][k] for k in ("families", "F1", "PPV", "TPR")) for r in scores)
    assert any(r["difference_pp"] < 0 for r in differences) and any(r["difference_pp"] > 0 for r in differences)


@pytest.mark.parametrize("fault", ("cohort", "duplicate", "truth", "root_type", "root_summary", "equality", "support",
                                 "search_score", "hit_row", "search_summary", "strata_duplicate", "strata_missing",
                                 "strata_size", "strata_nan", "strata_f1", "membership", "delta", "delta_duplicate"))
def test_corrupted_cohort_or_statistics_refuse(reports, fault):
    search, recon, strata = reports
    row = search["cases"][0]
    if fault == "cohort":
        search["changed_pairs"] += 1
    elif fault == "duplicate":
        search["cases"][1] = copy.deepcopy(row)
    elif fault == "truth":
        row["after"] = "TP"
    elif fault == "root_type":
        row["same_root_hog"] = "True"
    elif fault == "root_summary":
        recon["summary"]["FP_same_root"] += 1
    elif fault == "equality":
        row["selected_directed_score_multisets_identical"] = False
    elif fault == "support":
        row["direct_search"][plotter.CELLS[0]]["support"] = "wrong"
    elif fault in ("search_score", "hit_row"):
        hit = next(h for case in search["cases"] for h in case["direct_search"][plotter.CELLS[0]]["gene_a_to_b"])
        hit["score" if fault == "search_score" else "row"] = float("nan") if fault == "search_score" else -1
    elif fault == "search_summary":
        search["summary"][0]["pairs"] += 1
    elif fault == "strata_duplicate":
        strata["rows"][1] = strata["rows"][0]
    elif fault == "strata_missing":
        strata["rows"][0]["stratum"] = "wrong"
    elif fault == "strata_size":
        strata["rows"][0]["families"] += 1
    elif fault == "strata_nan":
        strata["rows"][0]["F1"] = float("nan")
    elif fault == "strata_f1":
        strata["rows"][0]["F1"] += .001
    elif fault == "membership":
        next(r for r in strata["rows"] if r["cell"] == plotter.CELLS[1] and r["stratum"] == "all")["family_members"][0] = "wrong"
    elif fault == "delta":
        strata["differences"][0]["PPV"] += .001
    elif fault == "delta_duplicate":
        strata["differences"][1] = strata["differences"][0]
    with pytest.raises(ValueError):
        plotter.figure_data(search, recon, strata)


@pytest.mark.parametrize("kind", plotter.KINDS)
def test_actual_report_readback_source_bindings(kind):
    refs = [plotter.record(RESULTS / path) for path in FILES[kind]]
    evidence = []
    assert plotter.paired_reports(kind, *refs, evidence)["schema"] == plotter.KINDS[kind][0]
    assert len(evidence) == 4


@pytest.mark.parametrize("fault", ("report", "source", "scope", "schema", "new_inference"))
def test_changed_readback_binding_refuses(tmp_path, fault):
    refs = [plotter.record(RESULTS / path) for path in FILES["strata"]]
    readback = json.loads(Path(refs[1]["path"]).read_text())
    if fault == "report":
        readback["report"]["sha256"] = "0" * 64
    elif fault == "source":
        readback["source"] = plotter.record(__file__)
    elif fault == "scope":
        readback["publication_ready"] = True
    elif fault == "schema":
        readback["schema"] = "wrong"
    elif fault == "new_inference":
        readback["new_uncertainty"] = True
    path = tmp_path / "readback.json"
    path.write_text(json.dumps(readback))
    with pytest.raises(ValueError):
        plotter.paired_reports("strata", refs[0], plotter.record(path), [])


def test_real_join_and_generated_tables_with_render_stub(tmp_path, monkeypatch):
    inputs = {kind: [plotter.record(RESULTS / path) for path in files] for kind, files in FILES.items()}
    calls = []
    monkeypatch.setattr(plotter, "render", lambda *args: calls.append(args))
    result = plotter.run(inputs, tmp_path / "figure")
    assert len(calls) == 1 and result["panels"] == 4 and result["raw_search_rescanned"] is False
    assert result["scientific_evidence_replayed"] is result["publication_ready"] is False
    for name in ("stages.tsv", "strata.tsv", "differences.tsv"):
        with (tmp_path / "figure" / name).open(newline="") as stream:
            assert len(list(csv.DictReader(stream, delimiter="\t"))) == 10
    assert all(plotter.record(ref["path"]) == ref for ref in result["outputs"])
    with pytest.raises(ValueError, match="occupied"):
        plotter.run(inputs, tmp_path / "figure")


def test_actual_render_legends_do_not_cover_plot_areas(reports, tmp_path, monkeypatch):
    held = []
    original_close = plotter.plt.close
    monkeypatch.setattr(plotter.plt, "close", lambda fig: held.append(fig))
    plotter.render(*plotter.figure_data(*reports), tmp_path)
    assert len(held) == 1
    fig = held[0]
    fig.canvas.draw()
    try:
        assert len(fig.axes) == 4
        assert all(not axis.get_legend().get_window_extent().overlaps(axis.bbox) for axis in fig.axes)
        assert all(not axis.get_legend().get_window_extent().overlaps(axis._left_title.get_window_extent())
                   for axis in fig.axes)
        assert [len(axis.patches) for axis in fig.axes] == [6, 4, 0, 10]
        assert all((tmp_path / ("native_swiss_mechanism." + suffix)).is_file() for suffix in ("png", "pdf", "svg"))
    finally:
        original_close(fig)
