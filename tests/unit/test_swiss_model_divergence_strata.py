"""Invented counts/features only: preserve frozen bins and the real statistic."""

import copy
import csv
import inspect
import math

import pytest

from benchmark_tools import export_swiss_model_divergence_strata as exporter
from benchmark_tools import readback_swiss_model_divergence_strata as reader


@pytest.fixture
def documents():
    names = [f"F{i:02d}" for i in range(18)]
    memberships = {f: [f"{f}_g{j:03d}" for j in range(31 if i < 17 else 36)]
                   for i, f in enumerate(names)}
    flags = {k: False for k in ("independent_confirmation", "publication_ready",
             "new_accuracy_or_resource_admission", "new_uncertainty", "scientific_timings_admitted")}
    flags["new_bootstrap_draws"] = 0
    features = dict(schema="swiss_model_divergence_features_v1", status="features_constructed_unverified",
                    failed_families=[], model="WAG+G4", seed=20261007,
                    unit="model_estimated_expected_amino_acid_substitutions_per_site",
                    prediction_statistics_evaluated=False, memberships=memberships,
                    median_family_distance=8.5,
                    strata=dict(all=names, lower_or_equal_median=names[:9], higher_than_median=names[9:]), **flags)
    admission = dict(schema="swiss_model_divergence_edge_readback_v1", status="features_verified",
                     families_checked=18, proteins_checked=563, scheduler=dict(state="COMPLETED", exit_code="0:0"),
                     strata=copy.deepcopy(features["strata"]),
                     features={f: dict(median_pair_distance=i) for i, f in enumerate(names)}, **flags)
    rows = []
    for cell_index, cell in enumerate(exporter.CELLS):
        for i, family in enumerate(names):
            positive, negative = 20 + 7 * i, 150 + 3 * i
            fn, fp = 5 + i % 3 + cell_index, 10 + 5 * i - cell_index
            counts = dict(TP=positive - fn, FN=fn, FP=fp, TN=negative - fp)
            rows.append(dict(cell=cell, family=family, counts_without_prior=counts,
                             **{m: float(v) for m, v in reader.point(counts).items()}))
    overall = [dict(suite="sequence", stratum="all", cell=c,
                    **{m: float(v) for m, v in reader.macro(
                       [reader.point(r["counts_without_prior"]) for r in rows if r["cell"] == c]).items()})
               for c in exporter.CELLS]
    counts = dict(schema="native_qfo_three_cell_strata_v1", memberships=memberships,
                  family_rows=rows, rows=overall,
                  cells=[dict(cell=c, **(dict(timing_eligible=False, timing_admitted=False)
                         if c == exporter.CELLS[1] else {})) for c in exporter.CELLS], **flags)
    count_reader = dict(family_rows_checked=54, score_rows_checked=60, differences_checked=40,
                        inherited_score_rows_reproduced=40, **flags)
    return features, admission, counts, count_reader


def test_complete_invented_projection_matches_independent_rational_arithmetic(documents):
    features, _, counts, _ = documents
    report = exporter.project(*documents)
    rows, diffs = reader.numerical_readback(report, counts, features)
    assert len(rows) == 9 and len(diffs) == 6 and len(report["family_rows"]) == 54
    assert [r["families"] for r in rows] == [18] * 3 + [9] * 6
    assert report["cells"][1]["timing_eligible"] is False
    assert report["new_bootstrap_draws"] == 0 and report["new_uncertainty"] is False
    for diff in report["differences"]:
        index = {(r["cell"], r["stratum"]): r for r in report["rows"]}
        for m in exporter.METRICS:
            assert diff[m + "_pp"] == pytest.approx(100 * (
                index[diff["candidate"], diff["stratum"]][m] - index[diff["reference"], diff["stratum"]][m]))


def test_macro_statistic_is_not_pooled_pairs_or_mean_family_f1(documents):
    report = exporter.project(*documents)
    families = [r for r in report["family_rows"] if r["cell"] == exporter.CELLS[0]]
    mean_f1 = sum(r["F1"] for r in families) / len(families)
    pooled = exporter.family_statistics({k: sum(r["counts_without_prior"][k] for r in families)
                                         for k in ("TP", "FP", "FN", "TN")})
    assert not math.isclose(report["rows"][0]["F1"], mean_f1, abs_tol=1e-5)
    assert not math.isclose(report["rows"][0]["F1"], pooled["F1"], abs_tol=1e-5)


@pytest.mark.parametrize("alteration", ["failed", "unverified", "mismatch", "bad_counts", "inflated", "partial"])
def test_invalid_scope_partial_verification_and_counts_rejected(documents, alteration):
    features, admitted, counts, counted = copy.deepcopy(documents)
    if alteration == "failed":
        features["failed_families"] = ["F00"]
    elif alteration == "unverified":
        admitted["status"] = "unverified"
    elif alteration == "mismatch":
        counts["memberships"]["F00"] = []
    elif alteration == "bad_counts":
        counts["family_rows"][0]["counts_without_prior"]["TP"] = -1
    elif alteration == "inflated":
        counted["publication_ready"] = True
    elif alteration == "partial":
        counted["differences_checked"] = 39
    with pytest.raises(ValueError):
        exporter.project(features, admitted, counts, counted)


def test_empty_bin_has_na_not_zero_and_no_rebalancing(documents):
    features, admission, counts, counted = copy.deepcopy(documents)
    features["strata"]["lower_or_equal_median"] = features["strata"]["all"]
    features["strata"]["higher_than_median"] = []
    admission["strata"] = copy.deepcopy(features["strata"])
    report = exporter.project(features, admission, counts, counted)
    reader.numerical_readback(report, counts, features)
    assert all(r["status"] == "empty" and r["F1"] is None for r in report["rows"][-3:])
    assert all(r["F1_pp"] is None for r in report["differences"][-2:])
    assert "| higher_than_median | 0 | p0_c0_r0 | NA | NA | NA |" in exporter.table(report)


@pytest.mark.parametrize("alteration", ["number", "members", "missing", "unit", "semantics", "count"])
def test_independent_reader_detects_table_corruption(documents, alteration):
    features, _, counts, _ = documents
    report = exporter.project(*documents)
    if alteration == "number":
        report["rows"][0]["F1"] += .01
    elif alteration == "members":
        report["rows"][3]["family_members"] = []
    elif alteration == "missing":
        report["differences"].pop()
    elif alteration == "unit":
        report["differences"][0]["F1_pp"] /= 100
    elif alteration == "semantics":
        report["rows"][1]["prediction_semantics"] = "group_clique"
    elif alteration == "count":
        report["family_rows"] = copy.deepcopy(report["family_rows"])
        report["family_rows"][0]["counts_without_prior"]["TN"] += 1
    with pytest.raises(ValueError):
        reader.numerical_readback(report, counts, features)


def test_tsv_and_human_tables_preserve_all_metrics_labels_and_units(documents, tmp_path):
    features, _, counts, _ = documents
    report = exporter.project(*documents)
    expected, _ = reader.numerical_readback(report, counts, features)
    path = tmp_path / "scores.tsv"
    fields = ("stratum", "cell", "families", "status", *exporter.METRICS, "prediction_semantics")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        writer.writerows(report["rows"])
    reader.verify_tsv(path, expected, fields, exporter.METRICS)
    text = exporter.table(report)
    assert text.count("| p0_c0_r0 |") == 3
    assert "Scores are percentages; differences are percentage points" in text
    assert "No subgroup confidence intervals" in text
    assert "failed timing remains ineligible" in text
    path.write_text(path.read_text().replace("descriptive", "confirmed", 1))
    with pytest.raises(ValueError):
        reader.verify_tsv(path, expected, fields, exporter.METRICS)


def test_projection_refuses_existing_directory_before_reading_any_inputs(tmp_path):
    with pytest.raises(ValueError, match="Existing projection"):
        exporter.export(tmp_path, "absent", "sha", "absent", "sha", tmp_path, "a" * 40)


def test_reader_is_independent_of_float_exporter():
    source = inspect.getsource(reader)
    assert "import export_swiss" not in source and "from benchmark_tools" not in source
    assert "Fraction" in source and ".distance(" not in source
