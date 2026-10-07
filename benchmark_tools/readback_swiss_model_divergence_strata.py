"""Independent rational verification of the fixed-model distance score tables."""

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
FALSE_FLAGS = ("new_uncertainty", "independent_confirmation", "new_accuracy_or_resource_admission",
               "publication_ready", "scientific_timings_admitted", "raw_scorer_repeated", "alignment_repeated")
PINS = dict(counts="55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5",
            counts_reader="6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969",
            protocol="a596371c443471585b9a277cbdd1466bb3f9870200f3187d20286461d3699bfd")


def assert_that(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def point(counts):
    assert_that(set(counts) == {"TP", "FP", "FN", "TN"}
                and all(type(v) is int and v >= 0 for v in counts.values())
                and sum(counts.values()) > 0, "Invalid retained integer counts")
    tp, fp, fn = (Fraction(counts[k] + 2, 2) for k in ("TP", "FP", "FN"))
    p, r = tp / (tp + fp), tp / (tp + fn)
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def macro(points):
    if not points:
        return dict.fromkeys(METRICS)
    p, r = (sum(row[m] for row in points) / len(points) for m in ("PPV", "TPR"))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def verify_value(actual, expected):
    assert_that(actual is None if expected is None else
                type(actual) in (float, int) and math.isfinite(actual)
                and math.isclose(actual, float(expected), rel_tol=0, abs_tol=1e-12),
                "Incorrect numeric or empty-bin statistic")


def numerical_readback(report, counts, features):
    assert_that(report["family_rows"] == counts["family_rows"], "Retained family rows changed")
    memberships = features["memberships"]
    assert_that(report["memberships"] == counts["memberships"] == memberships,
                "Member universes differ")
    assert_that([(r["cell"], r["family"]) for r in report["family_rows"]] ==
                [(c, f) for c in CELLS for f in sorted(memberships)], "Incomplete family count rows")
    index = {}
    for row in report["family_rows"]:
        index[row["cell"], row["family"]] = point(row["counts_without_prior"])
        for m in METRICS:
            verify_value(row[m], index[row["cell"], row["family"]][m])
    assert_that(report["bins"] == features["strata"] and report["median_family_distance"] ==
                features["median_family_distance"], "Fixed bins or cutoff changed")
    expected_rows, expected_diffs = [], []
    for name in ("all", "lower_or_equal_median", "higher_than_median"):
        families = report["bins"][name]
        values = {c: macro([index[c, f] for f in families]) for c in CELLS}
        for cell in CELLS:
            expected_rows.append(dict(cell=cell, stratum=name, families=len(families), family_members=families,
                                      status="descriptive" if families else "empty",
                                      prediction_semantics="native_pair" if cell == CELLS[1] else "group_clique",
                                      **values[cell]))
        for contrast, cell in (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2])):
            expected_diffs.append(dict(contrast=contrast, candidate=cell, reference=CELLS[0],
                                      stratum=name, families=len(families), family_members=families,
                                      status="descriptive" if families else "empty",
                                      **{m + "_pp": 100 * (values[cell][m] - values[CELLS[0]][m])
                                         if families else None for m in METRICS}))
    for key, expected, numeric in (("rows", expected_rows, METRICS),
                                   ("differences", expected_diffs, tuple(m + "_pp" for m in METRICS))):
        assert_that(len(report[key]) == len(expected), "Incomplete score/difference rows")
        for actual, computed in zip(report[key], expected):
            assert_that(set(actual) == set(computed), "Changed row schema")
            for field, value in computed.items():
                if field in numeric:
                    verify_value(actual[field], value)
                else:
                    assert_that(actual[field] == value, "Incorrect row membership/semantics/status")
    return expected_rows, expected_diffs


def verify_tsv(path, expected, fields, numeric):
    with Path(path).open(newline="", encoding="ascii") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert_that(reader.fieldnames == list(fields), "Wrong TSV fields")
        rows = list(reader)
    assert_that(len(rows) == len(expected), "Incomplete TSV rows")
    for actual, computed in zip(rows, expected):
        for field in fields:
            if field in numeric:
                value = None if actual[field] == "" else float(actual[field])
                verify_value(value, computed[field])
            else:
                assert_that(actual[field] == str(computed[field]), "Wrong TSV label/status/size")


def verify(path, repo, source_commit):
    path, repo = Path(path).resolve(), Path(repo).resolve()
    checked = [record(path)]
    report = json.loads(path.read_text())

    def checked_path(ref):
        assert_that(record(ref["path"]) == ref, "Changed direct projection binding")
        if ref not in checked:
            checked.append(ref)
        return Path(ref["path"])

    assert_that(report["schema"] == "native_qfo_swiss_model_divergence_strata_v1"
                and report["model"] == "WAG+G4" and report["seed"] == 20261007
                and report["unit"] == "model_estimated_expected_amino_acid_substitutions_per_site"
                and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
                and all(report[k] is False for k in FALSE_FLAGS), "Inflated projection scope")
    assert_that(set(report["inputs"]) == {*PINS, "features", "feature_reader"}, "Unexpected input set")
    documents = {}
    for key, ref in report["inputs"].items():
        if key in PINS:
            assert_that(ref["sha256"] == PINS[key], "Changed prospective pinned input")
        bound_path = checked_path(ref)
        if bound_path.suffix == ".json":
            documents[key] = json.loads(bound_path.read_text())
    feature_reader, count_reader = documents["feature_reader"], documents["counts_reader"]
    assert_that(feature_reader["report"] == report["inputs"]["features"]
                and feature_reader["status"] == "features_verified"
                and feature_reader["families_checked"] == 18 and feature_reader["proteins_checked"] == 563
                and feature_reader["scheduler"]["state"] == "COMPLETED"
                and feature_reader["scheduler"]["exit_code"] == "0:0"
                and count_reader["report"] == report["inputs"]["counts"]
                and count_reader["family_rows_checked"] == 54, "Incomplete inherited verification")
    assert_that(report["features"] == feature_reader["features"]
                and report["cells"] == documents["counts"]["cells"], "Descriptors/native identities changed")
    assert_that(report["cells"][1]["timing_eligible"] is False
                and report["cells"][1]["timing_admitted"] is False, "Failed timing promoted")
    for ref in report["inherited_sources"]:
        checked_path(ref)
    source = checked_path(report["source"])
    assert_that(source.read_bytes() == subprocess.check_output(["git", "show", report["source_commit"]
                + ":benchmark_tools/export_swiss_model_divergence_strata.py"], cwd=repo), "Unfrozen exporter")
    rows, diffs = numerical_readback(report, documents["counts"], documents["features"])
    verify_tsv(checked_path(report["outputs"]["rows"]), rows,
               ("stratum", "cell", "families", "status", *METRICS, "prediction_semantics"), METRICS)
    verify_tsv(checked_path(report["outputs"]["differences"]), diffs,
               ("stratum", "contrast", "families", "status", *(m + "_pp" for m in METRICS)),
               tuple(m + "_pp" for m in METRICS))
    text = checked_path(report["outputs"]["table"]).read_text()
    for row in rows:
        values = ["NA" if row[m] is None else f"{100 * float(row[m]):.3f}" for m in METRICS]
        line = "| " + " | ".join([row["stratum"], str(row["families"]), row["cell"], *values]) + " |"
        assert_that(text.count(line) == 1, "Missing/incorrect human score row")
    for row in diffs:
        values = ["NA" if row[m + "_pp"] is None else f"{float(row[m + '_pp']):+.3f}" for m in METRICS]
        line = "| " + " | ".join([row["stratum"], str(row["families"]), row["contrast"], *values]) + " |"
        assert_that(text.count(line) == 1, "Missing/incorrect human difference row")
    assert_that(Path(__file__).read_bytes() == subprocess.check_output(["git", "show", source_commit
                + ":benchmark_tools/readback_swiss_model_divergence_strata.py"], cwd=repo), "Unfrozen reader")
    return dict(schema="swiss_model_divergence_strata_rational_readback_v1",
                status="projection_verified", report=record(path), source=record(__file__),
                source_commit=source_commit, checked_inputs=checked, families_checked=18, proteins_checked=563,
                family_rows_checked=54, score_rows_checked=9, differences_checked=6,
                rational_macro_points=[{k: (float(v) if k in METRICS and v is not None else v)
                                        for k, v in row.items()} for row in rows],
                limitations=["Independent exact arithmetic; inherited raw counts and distance verification.",
                             "No new intervals, admission, true history, causal inference or generalization."],
                **{k: False for k in FALSE_FLAGS}, new_bootstrap_draws=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    assert_that(not args.output.exists() and not args.output.is_symlink(), "Existing readback; never overwrite")
    result = verify(args.report, args.repo, args.source_commit)
    with args.output.open("x", encoding="ascii") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "family_rows_checked", "score_rows_checked", "differences_checked")}))
