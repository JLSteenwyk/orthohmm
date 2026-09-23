import json
import csv
from fractions import Fraction
from pathlib import Path

import pytest

from benchmark_tools.export_swiss_duplication_strata import build_rows, export, REFERENCE

RESULTS = Path(__file__).resolve().parents[1] / "benchmark_tools/results"
COUNTS = RESULTS / "qfo_fastoma_swiss_uncertainty_22098.json"
FEATURES = RESULTS / "swiss_duplication_features_v2_20260923.json"


def inputs():
    return json.loads(COUNTS.read_text()), json.loads(FEATURES.read_text())


def test_all_methods_missing_and_full_statistic():
    counts, features = inputs()
    rows = build_rows(counts, features)
    assert len(rows) == 32
    for row in rows:
        if row["method"] == "orthomcl_1_4" or row["stratum"] == "missing_duplication_fraction":
            assert row["F1"] is row["delta_F1"] is None
        elif row["method"] == REFERENCE:
            assert row["delta_F1"] == row["delta_PPV"] == row["delta_TPR"] == 0
        if row["stratum"] == "all" and row["F1"] is not None:
            assert row["F1"] == pytest.approx(counts["point_estimates"][row["method"]]["F1"], abs=1e-12)


@pytest.mark.parametrize("problem", ["native", "scored", "bins", "median", "fraction", "count", "family"])
def test_changed_features_rejected(problem):
    counts, features = inputs()
    if problem == "native":
        features["independent_native_counts_checked"] = False
    elif problem == "scored":
        features["prediction_statistics_evaluated"] = True
    elif problem == "bins":
        features["primary_strata"]["lower_duplication_fraction"].pop()
    elif problem == "median":
        features["median_fraction"] = "0"
    elif problem == "fraction":
        features["families"]["APP"]["duplication_fraction"] = "0"
    elif problem == "count":
        features["families"]["APP"]["informative_nodes"] += 1
    else:
        del features["families"]["APP"]
    with pytest.raises(ValueError):
        build_rows(counts, features)


def test_export_and_no_overwrite(tmp_path):
    report = export(COUNTS, FEATURES, tmp_path / "table")
    assert len(report["rows"]) == 32 and report["new_inferential_claims"] is False
    assert len((tmp_path / "table/scores.tsv").read_text().splitlines()) == 33
    assert "not evolutionary duplication rates" in (tmp_path / "table/scores.md").read_text()
    with pytest.raises(FileExistsError):
        export(COUNTS, FEATURES, tmp_path / "table")


def test_independent_rational_reproduction_of_export(tmp_path):
    export(COUNTS, FEATURES, tmp_path / "table")
    counts, features = inputs()
    groups = {"all": counts["families"], **features["primary_strata"]}
    expected = {}
    for method in counts["reconstructed_counts"]["methods"]:
        for group, names in groups.items():
            values = []
            if method["status"] == "counts_verified":
                for family in method["families"]:
                    if family["family"] in names:
                        c = family["counts_without_prior"]
                        values.append((Fraction(c["TP"]+2, c["TP"]+c["FP"]+4),
                                       Fraction(c["TP"]+2, c["TP"]+c["FN"]+4)))
            if values:
                p = sum(v[0] for v in values)/len(values)
                r = sum(v[1] for v in values)/len(values)
                expected[method["method"], group] = dict(F1=2*p*r/(p+r), PPV=p, TPR=r)
            else:
                expected[method["method"], group] = dict.fromkeys(("F1", "PPV", "TPR"))
    with (tmp_path / "table/scores.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            for metric, value in expected[row["method"], row["stratum"]].items():
                ref = expected[REFERENCE, row["stratum"]][metric]
                delta = None if value is None or ref is None else value-ref
                for key, exact in ((metric, value), ("delta_"+metric, delta)):
                    if exact is None:
                        assert row[key] == "NA"
                    else:
                        assert float(row[key]) == pytest.approx(float(exact), rel=0, abs=1e-12)
