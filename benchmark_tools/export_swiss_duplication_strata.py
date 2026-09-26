"""Join prespecified mapped-tree features to admitted descriptive method scores."""

import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path

from benchmark_tools.export_swiss_descriptive_strata import build, COUNTS_SHA, METRICS
from benchmark_tools.extract_swiss_duplication_features import bins
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

FEATURE_SHA = "97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec"
HELPER_SHA = "5c91edf93b0ab0d563a06a7c90826a8a9c144994bff3dffed755930762a1d121"
FEATURE_HELPER_SHA = "96355e41e5733dbf049e113e348bd6819cc215448e936b7d0e99070908d4b73e"
REFERENCE = "orthofinder_3_1_5_full"


def build_rows(counts, features):
    if (features["status"] != "mapped_tree_duplication_features_unscored"
            or features["independent_native_counts_checked"] is not True
            or features["prediction_statistics_evaluated"] is not False
            or features["publication_ready"] is not False):
        raise ValueError("Require unscored native-checked features")
    for row in features["families"].values():
        keys = ("explicit_duplication_nodes", "explicit_speciation_nodes", "default_speciation_nodes",
                "child_overlap_nodes", "informative_nodes", "mapped_members")
        if any(type(row[k]) is not int or row[k] < 0 for k in keys):
            raise ValueError("Invalid feature counts")
        n = row["informative_nodes"]
        if sum(row[k] for k in keys[:3]) != n:
            raise ValueError("Feature counts do not sum")
        expected = None if n == 0 else str(Fraction(row["explicit_duplication_nodes"], n))
        if row["duplication_fraction"] != expected:
            raise ValueError("Feature fraction differs from counts")
    median, groups = bins(features["families"])
    if median != features["median_fraction"] or groups != features["primary_strata"]:
        raise ValueError("Changed prespecified strata")
    rows = build(counts, dict(family_memberships=features["families"],
                             primary_strata=groups, secondary_strata={}))
    references = {row["stratum"]: row for row in rows if row["method"] == REFERENCE}
    for row in rows:
        for metric in METRICS:
            ref = references[row["stratum"]][metric]
            row["delta_" + metric] = None if row[metric] is None or ref is None else row[metric] - ref
    return rows


def export(counts, features_path, output, *, counts_sha=COUNTS_SHA):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    features = read_frozen(features_path, FEATURE_SHA)
    helpers = [record(Path(__file__).with_name(name)) for name in
               ("export_swiss_descriptive_strata.py", "extract_swiss_duplication_features.py")]
    if [r["sha256"] for r in helpers] != [HELPER_SHA, FEATURE_HELPER_SHA]:
        raise ValueError("Changed feature or statistic implementation")
    inputs = [record(counts), record(features_path), record(__file__), *helpers,
              *features["checked_inputs"]]
    for item in inputs:
        check(item)
    rows = build_rows(read_frozen(counts, counts_sha), features)
    output.mkdir(parents=True)
    fields = ["method", "stratum", "families", "status", *METRICS,
              *["delta_" + k for k in METRICS], "prediction_semantics"]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    lines = ["# Descriptive Duplication-Annotation Strata", "",
        "Scores are percentages; differences are percentage points versus full OrthoFinder.",
        "F1 is the harmonic mean of macro precision and recall, not mean family F1.",
        "No new intervals or significance tests. Missing is not zero.", ""]
    for name in dict.fromkeys(row["stratum"] for row in rows):
        selected = [row for row in rows if row["stratum"] == name]
        lines += [f"## {name} ({selected[0]['families']} families)", "",
                  "| Method | F1 | Precision | Recall | F1 difference | Status |",
                  "|---|---:|---:|---:|---:|---|"]
        for row in selected:
            values = ["NA" if row[k] is None else f"{100*row[k]:.3f}" for k in (*METRICS, "delta_F1")]
            lines.append("| " + " | ".join([row["method"], *values, row["status"]]) + " |")
        lines.append("")
    lines.append("Reference-derived, development-exposed annotation fractions are not evolutionary duplication rates or causal explanations. Default-S nodes are not explicit speciation observations. FastOMA uses a supplied tree; sequence-only OrthoFinder is a group-clique diagnostic.")
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    for item in inputs:
        check(item)
    manifest = dict(status="descriptive_swiss_duplication_strata_exported", inputs=inputs,
        rows=rows, outputs=[record(p) for p in sorted(output.iterdir())],
        new_inferential_claims=False, publication_ready=False)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "features", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--counts-sha256", default=COUNTS_SHA)
    args = parser.parse_args()
    export(args.counts.resolve(), args.features.resolve(), args.output.absolute(), counts_sha=args.counts_sha256)
