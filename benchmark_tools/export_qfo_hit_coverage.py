"""Export a compact, source-bound search diagnostic, not orthology accuracy."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

REPORT_SHA = "0e2b668fc5b2829f0a73175245e6c491986b0f2b1b89b2a4111776917d521e0f"
VARIANTS = ("hmm", "all_hits", "top100")
FIELDS = ("genes", "directed_hits", "self_hits", "nonself_hits", "cross_species_hits",
          "queries_without_hits", "queries_without_nonself_hits",
          "queries_without_cross_species_hits", "targets_without_hits")


def validate(report):
    if (report["status"] != "corrected_qfo_label_free_hit_comparison_complete"
            or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False
            or set(report["searches"]) != set(VARIANTS)):
        raise ValueError("Require complete, unscored search diagnostic")
    searches = report["searches"]
    labels = report["species_labels"]
    if len(labels) != len(set(labels)) or not labels:
        raise ValueError("Invalid species labels")
    n = searches["hmm"]["genes"]
    for row in searches.values():
        if any(type(row[k]) is not int or row[k] < 0 for k in FIELDS):
            raise ValueError("Invalid search counts")
        if (row["genes"] != n or n < 1
                or row["directed_hits"] != row["self_hits"] + row["nonself_hits"]
                or row["cross_species_hits"] > row["nonself_hits"]
                or row["self_hits"] > n
                or not 0 <= row["queries_without_hits"] <= row["queries_without_nonself_hits"]
                    <= row["queries_without_cross_species_hits"] <= n
                or row["targets_without_hits"] > n):
            raise ValueError("Inconsistent search counts")
        directions = row["species_directions"]
        expected = {(a, b) for a in range(len(labels)) for b in range(len(labels))}
        if (len(directions) != len(expected)
                or {(r["query_species"], r["target_species"]) for r in directions} != expected
                or sum(r["directed_hits"] for r in directions) != row["directed_hits"]
                or sum(r["directed_hits"] for r in directions
                       if r["query_species"] != r["target_species"]) != row["cross_species_hits"]
                or sum(row["normalized_score_histogram"]["counts"]) != row["directed_hits"]):
            raise ValueError("Direction or score histogram totals differ")
    pairs = (("hmm", "all_hits"), ("hmm", "top100"), ("all_hits", "top100"))
    if set(report["overlaps"]) != {a + "_vs_" + b for a, b in pairs}:
        raise ValueError("Wrong overlap comparisons")
    for a, b in pairs:
        for label, field in (("all", "directed_hits"), ("nonself", "nonself_hits")):
            row = report["overlaps"][a + "_vs_" + b][label]
            common, left, right = row["intersection"], searches[a][field], searches[b][field]
            if type(common) is not int or not 0 <= common <= min(left, right):
                raise ValueError("Invalid intersection")
            union = left + right - common
            expected = {"first_only": left - common, "second_only": right - common,
                        "union": union, "jaccard": common / union if union else None,
                        "fraction_of_first_recovered": common / left if left else None,
                        "fraction_of_second_recovered": common / right if right else None}
            if any(row[key] != value for key, value in expected.items()):
                raise ValueError("Overlap arithmetic differs")
    if report["overlaps"]["all_hits_vs_top100"]["all"]["second_only"]:
        raise ValueError("Top100 is not a subset of all hits")


def export(source, output):
    if output.exists():
        raise FileExistsError(output)
    report = read_frozen(source, REPORT_SHA)
    validate(report)
    evidence = [record(source), report["source"], *report["helpers"],
                report["canonicalization"], *report["checkpoints"].values()]
    for item in evidence:
        check(item)
    rows = [{"variant": name, **{k: report["searches"][name][k] for k in FIELDS}}
            for name in VARIANTS]
    limits = ["Directed search hits are not ortholog predictions or accuracy estimates.",
              "Equal E-values do not match sensitivity, calibration or computational effort.",
              "Top100 is a post-search per-query/per-target-species reporting cap, not an HMM prefilter.",
              "Only report/source/helper/manifest hashes and summary arithmetic are rechecked here; the full hit arrays are not rescored.",
              "The source job checks full checkpoint integrity; this export is not an independent replication of its hit-set intersections.",
              "Missing hits do not distinguish prefilter rejection from scoring rejection."]
    output.mkdir(parents=True)
    with (output / "coverage.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, ["variant", *FIELDS], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = {"status": "source_bound_label_free_coverage_summary", "accuracy_evaluated": False,
              "publication_ready": False, "rows": rows, "overlaps": report["overlaps"],
              "checked_records": evidence, "source": record(__file__), "limitations": limits,
              "table": record(output / "coverage.tsv")}
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.source.resolve(), args.output.resolve())
