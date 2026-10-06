"""Independent integer-bin and raw-count readback of native duplication strata."""

import argparse
import csv
from functools import cmp_to_key
import gzip
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_sequence_strata import (
    CELLS, METRICS, equal, raw_counts, record, require, stats,
)

NAMES = ("all", "lower_duplication_fraction", "upper_duplication_fraction", "missing_duplication_fraction")
FEATURE_SHA = "97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec"
IDENTIFIERS_SHA = "1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d"


def check_memberships(memberships, mapping, identifiers):
    require(set(mapping["families"]) == set(memberships), "Wrong mapped family inventory")
    for name, genes in memberships.items():
        original = mapping["families"][name]
        expected = set()
        for value in original["mapped_labels"].values():
            require(type(value) is int and value > 0, "Invalid retained entry ID")
            expected.add(value)
        observed = set()
        for gene in genes:
            number = identifiers.get(gene)
            require(type(number) is int and number > 0 and number not in observed, "Unknown, invalid or duplicate entry ID")
            observed.add(number)
        require(observed == expected and len(observed) == original["mapped_members"]
                and original["exact_match"] is True, "Mapped-entry inventory disagrees")


def rational_text(n, d):
    divisor = math.gcd(n, d)
    n, d = n // divisor, d // divisor
    return str(n) if d == 1 else f"{n}/{d}"


def independent_bins(memberships, feature):
    require(memberships and set(memberships) == set(feature["families"]), "Changed family universe")
    genes = [g for members in memberships.values() for g in members]
    require(all(members for members in memberships.values()) and len(set(genes)) == len(genes), "Changed gene universe")
    ratios = {}
    for family, value in feature["families"].items():
        keys = ("explicit_duplication_nodes", "explicit_speciation_nodes", "default_speciation_nodes",
                "informative_nodes", "child_overlap_nodes", "mapped_members")
        require(all(type(value[k]) is int and value[k] >= 0 for k in keys), "Invalid annotation counts")
        numerator, denominator = value[keys[0]], value[keys[3]]
        require(numerator + value[keys[1]] + value[keys[2]] == denominator
                and value[keys[4]] <= denominator and value[keys[5]] == len(memberships[family]), "Invalid annotation totals")
        require(value["duplication_fraction"] == (rational_text(numerator, denominator) if denominator else None),
                "Invalid stored fraction")
        ratios[family] = (numerator, denominator) if denominator else None
    # Compare integer cross-products, rather than the primary Fraction/median code.
    def compare(a, b):
        difference = a[0] * b[1] - b[0] * a[1]
        return (difference > 0) - (difference < 0)
    values = sorted([v for v in ratios.values() if v is not None], key=cmp_to_key(compare))
    middle = None
    if values:
        left, right = values[(len(values) - 1) // 2], values[len(values) // 2]
        middle = (left[0] * right[1] + right[0] * left[1], 2 * left[1] * right[1])
    require(feature["median_fraction"] == (rational_text(*middle) if middle else None), "Wrong median")
    bins = {name: [] for name in NAMES}
    for family, value in sorted(ratios.items()):
        bins[NAMES[0]].append(family)
        target = NAMES[3] if value is None else NAMES[2 if compare(value, middle) > 0 else 1]
        bins[target].append(family)
    require({k: bins[k] for k in NAMES[1:]} == feature["primary_strata"], "Wrong original bins")
    return bins


def verify_projection(report, counts, feature):
    bins = independent_bins(report["memberships"], feature)
    require(report["bins"] == bins, "Wrong reported bins")
    require([(r["cell"], r["stratum"]) for r in report["rows"]] ==
            [(cell, name) for cell in CELLS for name in NAMES], "Incomplete projected rows")
    require([(r["cell"], r["family"]) for r in report["family_rows"]] ==
            [(cell, family) for cell in CELLS for family in sorted(report["memberships"])], "Incomplete family rows")
    for row in report["family_rows"]:
        require(row["counts_without_prior"] == counts[row["cell"]][row["family"]]
                and all(equal(row[k], v) for k, v in stats(counts[row["cell"]], [row["family"]]).items()),
                "Wrong family counts/statistics")
    for row in report["rows"]:
        members = bins[row["stratum"]]
        require(row["family_members"] == members and row["families"] == len(members)
                and row["status"] == ("descriptive" if members else "empty_bin")
                and row["prediction_semantics"] == ("group_clique" if row["cell"] == CELLS[0] else "resolved_native_pairs")
                and all(equal(row[k], v) for k, v in stats(counts[row["cell"]], members).items()), "Wrong projected statistic")
    require([r["stratum"] for r in report["differences"]] == list(NAMES), "Incomplete differences")
    for row in report["differences"]:
        members = bins[row["stratum"]]
        left, right = (stats(counts[cell], members) for cell in CELLS)
        require(row["family_members"] == members and row["families"] == len(members)
                and row["status"] == ("descriptive" if members else "empty_bin")
                and all(equal(row[k], None if not members else right[k] - left[k]) for k in METRICS), "Wrong difference")


def verify(path, digest):
    report_ref = record(path)
    require(report_ref["sha256"] == digest, "Changed report")
    report = json.loads(Path(path).read_text())
    require(report["schema"] == "native_qfo_swiss_duplication_strata_v1"
            and report["source"] == record(Path(__file__).with_name("export_native_qfo_swiss_duplication_strata.py"))
            and report["features"]["sha256"] == FEATURE_SHA
            and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
            and all(report[k] is False for k in ("original_tree_traversal_repeated", "new_uncertainty",
                "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")), "Wrong report scope")
    checked = [report_ref, report["source"], *report["checked_inputs"], *report["outputs"],
               record(Path(__file__).with_name("readback_native_qfo_swiss_sequence_strata.py"))]
    require(report["identifiers"]["sha256"] == IDENTIFIERS_SHA
            and all(ref in checked for ref in (report["features"], report["mapping"], report["native"], report["identifiers"],
                report["native_readback"], report["original_protocol"])), "Missing supplied identity")
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed supplied evidence")
    native = json.loads(Path(report["native"]["path"]).read_text())
    feature = json.loads(Path(report["features"]["path"]).read_text())
    mapping = json.loads(Path(report["mapping"]["path"]).read_text())
    require(report["family_rows"] == native["family_rows"] and len(report["memberships"]) == 18
            and sum(map(len, report["memberships"].values())) == 563,
            "Changed native/mapping inventory")
    with gzip.open(report["identifiers"]["path"], "rt") as stream:
        identifiers = json.load(stream)["mapping"]
    check_memberships(report["memberships"], mapping, identifiers)
    del identifiers
    require(report["inherited_feature_inputs"] == [r for r in feature["checked_inputs"] if r not in report["checked_inputs"]],
            "Wrong inherited evidence inventory")
    binding = json.loads(Path(native["binding"]["path"]).read_text())
    require([c["cell"] for c in report["cells"]] == list(CELLS), "Wrong cell inventory")
    counts, truth_anchor, raw_rows = {}, None, 0
    for cell in report["cells"]:
        name = cell["cell"]
        audit_ref = binding["bound_cells"][name]["count_audit"]
        require(audit_ref in checked and cell["raw_file"] in checked, "Missing native raw/audit binding")
        selected = [r for r in json.loads(Path(audit_ref["path"]).read_text())["cells"] if r["cell"] == name]
        require(len(selected) == 1 and all(cell[k] == selected[0][k] for k in cell), "Changed native cell binding")
        if name == CELLS[1]:
            require(cell["timing_eligible"] is False and cell["timing_admitted"] is False
                    and selected[0]["resources"] is None, "Failed timing relabeled")
        counts[name], truth = raw_counts(cell["raw_file"]["path"], report["memberships"])
        require(len(truth) == 10765 and (truth_anchor is None or truth_anchor == truth), "Different raw truth/coverage")
        truth_anchor = truth
        raw_rows += len(truth)
    verify_projection(report, counts, feature)
    require(raw_rows == report["raw_relations_checked"] == 21530, "Wrong raw row count")
    table_ref = next(r for r in report["outputs"] if Path(r["path"]).name == "scores.tsv")
    with open(table_ref["path"], newline="") as stream:
        table = csv.DictReader(stream, delimiter="\t")
        fields = ["cell", "stratum", "families", "status", *METRICS, "prediction_semantics"]
        require(table.fieldnames == fields, "Wrong TSV columns")
        observed = list(table)
    require(observed == [{k: "NA" if row[k] is None else str(row[k]) for k in fields} for row in report["rows"]],
            "Changed TSV rows")
    for ref in checked:
        require(record(ref["path"]) == ref, "Evidence changed during readback")
    return dict(schema="native_qfo_swiss_duplication_strata_readback_v1", report=report_ref, source=record(__file__),
                checked_inputs=checked, bins=report["bins"], raw_rows_checked=raw_rows,
                family_rows_checked=len(report["family_rows"]), projection_rows_checked=len(report["rows"]),
                differences_checked=len(report["differences"]), original_tree_traversal_repeated=False,
                new_uncertainty=False, new_accuracy_or_resource_admission=False,
                independent_confirmation=False, publication_ready=False,
                limitations=["Integer cross-product bin decisions and alternative raw/statistic arithmetic on the same original data.",
                    "Original mapping/acquisition/Darwin traversal inherited, not fresh biological or transitive validation.",
                    "No new intervals, causal annotation mechanism, total-HMM control or timing repair."])


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
