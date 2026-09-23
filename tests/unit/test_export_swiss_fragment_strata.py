import json
from pathlib import Path

import pytest

from benchmark_tools.export_swiss_fragment_strata import build_rows, export, REFERENCE
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_swiss_historical_fragments import family_bins


def inputs():
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    counts = json.loads((results / "qfo_fastoma_swiss_uncertainty_22098.json").read_text())
    families = {f: [f] for f in counts["families"]}
    annotations = {f: dict(fragment_flag=False, incomplete_sequence_features=[],
        selection_class="baseline_release") for f in families}
    first, second = list(families)[:2]
    annotations[first].update(fragment_flag=True, selection_class="later_sequence_version")
    annotations[second] = None
    admission = dict(status="historical_annotation_panel_checked_with_explicit_missingness",
        annotation_panel_admitted=True, prediction_statistics_evaluated=False, publication_ready=False,
        families=families, annotations=annotations, matched=17, missing=1,
        strata=family_bins(families, annotations), baseline_only_strata=family_bins(families, annotations, True))
    return counts, admission


def test_all_methods_bins_and_missing_are_preserved():
    counts, admission = inputs()
    rows = build_rows(counts, admission)
    assert len(rows) == 56
    for row in rows:
        if row["method"] == "orthomcl_1_4" or row["stratum"] == "baseline_only_annotation_positive":
            assert row["F1"] is row["delta_F1"] is None
        elif row["method"] == REFERENCE:
            assert row["delta_F1"] == row["delta_PPV"] == row["delta_TPR"] == 0
        if row["stratum"] == "all" and row["F1"] is not None:
            assert row["F1"] == pytest.approx(counts["point_estimates"][row["method"]]["F1"], abs=1e-12)


@pytest.mark.parametrize("problem", ["status", "scored", "admitted", "bins", "baseline", "coverage"])
def test_invalid_annotations_rejected(problem):
    counts, admission = inputs()
    if problem == "status":
        admission["status"] = "collecting"
    elif problem == "scored":
        admission["prediction_statistics_evaluated"] = True
    elif problem == "admitted":
        admission["annotation_panel_admitted"] = False
    elif problem in ("bins", "baseline"):
        admission["strata" if problem == "bins" else "baseline_only_strata"]["all_matched_unflagged"].pop()
    else:
        admission["missing"] = 0
    with pytest.raises(ValueError):
        build_rows(counts, admission)


def test_export_files_and_no_overwrite(tmp_path):
    _, admission = inputs()
    admission["records"] = []
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(admission))
    counts = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json"
    output = tmp_path / "export"
    report = export(counts, path, record(path)["sha256"], output)
    assert len(report["rows"]) == 56
    assert report["new_inferential_claims"] is False
    for item in report["outputs"]:
        check(item)
    text = (output / "scores.md").read_text()
    assert "Unflagged is not proven complete" in text
    assert "NA" in text
    assert len((output / "scores.tsv").read_text().splitlines()) == 57
    with pytest.raises(FileExistsError):
        export(counts, path, record(path)["sha256"], output)
    with pytest.raises(ValueError):
        export(counts, path, "0"*64, tmp_path / "bad_export")
    assert not (tmp_path / "bad_export").exists()
