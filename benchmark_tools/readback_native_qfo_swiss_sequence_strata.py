"""Stdlib-only independent native SwissTrees frozen-bin statistic readback."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

CELLS = ("p0_c0_r0", "p0_c0_r1")
METRICS = ("F1", "PPV", "TPR")
STRATA_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
HEADER = "# Dataset<tab>Protein ID 1<tab>Protein ID 2<tab>Correctness (TP:True positive, FP: False positive, TN: True Negative. FN: False Negative)"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def raw_counts(path, memberships):
    counts = {f: Counter({k: 0 for k in ("TP", "FP", "FN", "TN")}) for f in memberships}
    relations, members = {}, {f: set() for f in memberships}
    with gzip.open(path, "rt", newline="") as stream:
        require(stream.readline().rstrip("\r\n") == HEADER, "Wrong raw header")
        for row in csv.reader(stream, delimiter="\t"):
            require(len(row) == 4, "Wrong raw width")
            family, a, b, label = row
            require(family in counts and a and b and a != b and label in counts[family], "Invalid raw relation")
            key = (family, *sorted((a, b)))
            require(key not in relations, "Duplicate raw relation")
            relations[key] = label in ("TP", "FN")
            counts[family][label] += 1
            members[family].update((a, b))
    require(all(sorted(members[f]) == sorted(memberships[f]) and sum(counts[f].values()) > 0
                for f in memberships), "Wrong raw family membership/coverage")
    return counts, relations


def stats(counts, members):
    if not members:
        return dict.fromkeys(METRICS)
    precision = math.fsum((counts[f]["TP"] + 2) / (counts[f]["TP"] + counts[f]["FP"] + 4)
                          for f in members) / len(members)
    recall = math.fsum((counts[f]["TP"] + 2) / (counts[f]["TP"] + counts[f]["FN"] + 4)
                       for f in members) / len(members)
    return dict(F1=2 / (1 / precision + 1 / recall), PPV=precision, TPR=recall)


def equal(actual, expected):
    return (actual is None if expected is None else type(actual) in (int, float)
            and math.isfinite(actual) and math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12))


def verify(path, digest):
    report_ref = record(path)
    require(report_ref["sha256"] == digest, "Changed projection report")
    report = json.loads(Path(path).read_text())
    require(report["schema"] == "native_qfo_swiss_sequence_strata_v1"
            and report["source"] == record(Path(__file__).with_name("export_native_qfo_swiss_sequence_strata.py"))
            and report["strata"]["sha256"] == STRATA_SHA
            and all(report[k] is False for k in ("new_uncertainty", "new_accuracy_or_resource_admission",
                                                "independent_confirmation", "publication_ready")),
            "Wrong projection source/scope")
    checked = [report_ref, report["source"], report["binding"], report["strata"], *report["checked_inputs"], *report["outputs"]]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed directly supplied evidence")
    feature = json.loads(Path(report["strata"]["path"]).read_text())
    memberships = feature["family_memberships"]
    bins = {"all": sorted(memberships), **feature["primary_strata"], **feature["secondary_strata"]}
    require(len(memberships) == 18 and len(bins) == 11, "Wrong feature universe")
    binding = json.loads(Path(report["binding"]["path"]).read_text())
    expected_keys = [(cell, name) for cell in CELLS for name in bins]
    require([(r["cell"], r["stratum"]) for r in report["rows"]] == expected_keys, "Incomplete projection rows")
    counts_by_cell, truth_anchor, raw_rows = {}, None, 0
    for cell in CELLS:
        ref = binding["bound_cells"][cell]["count_audit"]
        require(ref in checked, "Missing original count audit binding")
        audit = json.loads(Path(ref["path"]).read_text())
        selected = [r for r in audit["cells"] if r["cell"] == cell]
        require(len(selected) == 1 and selected[0]["raw_file"] in checked, "Missing native raw binding")
        counts, truth = raw_counts(selected[0]["raw_file"]["path"], memberships)
        require(truth_anchor is None or truth_anchor == truth, "Different raw reference truth")
        truth_anchor = truth
        raw_rows += len(truth)
        counts_by_cell[cell] = counts
        expected = stats(counts, list(memberships))
        require(all(equal(selected[0]["aggregate"][k], expected[k]) for k in METRICS), "Native aggregate differs")
    expected_family_keys = [(cell, family) for cell in CELLS for family in sorted(memberships)]
    require([(r["cell"], r["family"]) for r in report["family_rows"]] == expected_family_keys, "Incomplete family rows")
    for row in report["family_rows"]:
        counts = counts_by_cell[row["cell"]]
        require(row["counts_without_prior"] == counts[row["family"]]
                and all(equal(row[k], value) for k, value in stats(counts, [row["family"]]).items()),
                "Native family count/statistic differs")
    for row in report["rows"]:
        members = bins[row["stratum"]]
        require(row["family_members"] == members and row["families"] == len(members)
                and row["status"] == ("descriptive" if members else "empty_bin")
                and row["prediction_semantics"] == ("group_clique" if row["cell"] == CELLS[0] else "resolved_native_pairs")
                and all(equal(row[k], value) for k, value in stats(counts_by_cell[row["cell"]], members).items()),
                "Native stratum statistic/membership/semantics differs")
    require([r["stratum"] for r in report["differences"]] == list(bins), "Incomplete descriptive differences")
    for row in report["differences"]:
        members = bins[row["stratum"]]
        left, right = (stats(counts_by_cell[cell], members) for cell in CELLS)
        require(row["families"] == len(members) and row["family_members"] == members
                and row["status"] == ("descriptive" if members else "empty_bin")
                and all(equal(row[k], None if left[k] is None else right[k] - left[k]) for k in METRICS),
                "Wrong descriptive difference")
    table_ref = next(ref for ref in report["outputs"] if Path(ref["path"]).name == "scores.tsv")
    with open(table_ref["path"], newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == ["cell", "stratum", "families", "status", *METRICS, "prediction_semantics"],
                "Wrong TSV columns")
        table = list(reader)
    require(len(table) == len(report["rows"]), "Incomplete TSV")
    for row, observed in zip(report["rows"], table):
        require(observed == {k: "NA" if row[k] is None else str(row[k]) for k in observed}, "TSV differs")
    require(raw_rows == report["raw_relations_checked"], "Wrong raw relation count")
    for ref in checked:
        require(record(ref["path"]) == ref, "Evidence changed during independent readback")
    return dict(schema="native_qfo_swiss_sequence_strata_readback_v1", report=report_ref, source=record(__file__),
        checked_inputs=checked, raw_rows_checked=raw_rows, families_checked=len(memberships),
        family_rows_checked=len(report["family_rows"]), projection_rows_checked=len(report["rows"]),
        differences_checked=len(report["differences"]), new_uncertainty=False,
        new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False,
        limitations=["Independent stdlib raw-row counts and alternative arithmetic on the same original data.",
                     "Direct diagnostic bindings, not repeated transitive admission or new biological replication.",
                     "No new uncertainty, causal subgroup claim, inference timing or total-HMM control."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("report", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("raw_rows_checked", "family_rows_checked", "projection_rows_checked", "differences_checked")}))
