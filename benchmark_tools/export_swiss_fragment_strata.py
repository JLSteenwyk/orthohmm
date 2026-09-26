"""Export descriptive historical fragment strata only after annotation admission."""

import argparse
import csv
import json
from pathlib import Path

from benchmark_tools.export_swiss_descriptive_strata import build, COUNTS_SHA, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_swiss_historical_fragments import family_bins

HELPER_SHA = "5c91edf93b0ab0d563a06a7c90826a8a9c144994bff3dffed755930762a1d121"
VERIFIER_SHA = "b316148f54992552a187ee9aeea8d5dbae41abf05665a650ced62f0d1f386a50"
REFERENCE = "orthofinder_3_1_5_full"


def build_rows(counts, admission):
    if (admission["status"] != "historical_annotation_panel_checked_with_explicit_missingness"
            or admission["annotation_panel_admitted"] is not True
            or admission["prediction_statistics_evaluated"] is not False
            or admission["publication_ready"] is not False):
        raise ValueError("Require admitted unscored annotations")
    families, annotations = admission["families"], admission["annotations"]
    matched = sum(r is not None for r in annotations.values())
    if (matched, len(annotations)-matched) != (admission["matched"], admission["missing"]):
        raise ValueError("Annotation coverage mismatch")
    combined = {}
    for baseline, key, prefix in ((False, "strata", "historical"), (True, "baseline_only_strata", "baseline_only")):
        bins = family_bins(families, annotations, baseline)
        if bins != admission[key]:
            raise ValueError("Changed frozen annotation bins")
        combined.update({prefix + "_" + name: members for name, members in bins.items()})
    rows = build(counts, dict(family_memberships=families, primary_strata=combined, secondary_strata={}))
    references = {r["stratum"]: r for r in rows if r["method"] == REFERENCE}
    for row in rows:
        for metric in METRICS:
            reference = references[row["stratum"]][metric]
            row["delta_" + metric] = None if row[metric] is None or reference is None else row[metric] - reference
    return rows


def export(counts, admission, admission_sha, output, *, counts_sha=COUNTS_SHA):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    helpers = [record(Path(__file__).with_name(name)) for name in
               ("export_swiss_descriptive_strata.py", "verify_swiss_historical_fragments.py")]
    if [r["sha256"] for r in helpers] != [HELPER_SHA, VERIFIER_SHA]:
        raise ValueError("Changed statistic or annotation-bin implementation")
    features = read_frozen(admission, admission_sha)
    inputs = [record(counts), record(admission), record(__file__), *helpers, *features["records"]]
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
    lines = ["# Historical Fragment Annotation Strata", "",
        "Descriptive percentages; differences are percentage points versus full OrthoFinder.",
        "F1 is the harmonic mean of macro precision and recall. No new intervals or significance tests.",
        f"Annotation coverage: {features['matched']} matched; {features['missing']} missing.", ""]
    for name in dict.fromkeys(r["stratum"] for r in rows):
        selected = [r for r in rows if r["stratum"] == name]
        lines.extend([f"## {name} ({selected[0]['families']} families)", "",
            "| Method | F1 | Precision | Recall | F1 difference | Status |",
            "|---|---:|---:|---:|---:|---|"])
        for row in selected:
            values = ["NA" if row[k] is None else f"{100*row[k]:.3f}" for k in (*METRICS, "delta_F1")]
            lines.append("| " + " | ".join([row["method"], *values, row["status"]]) + " |")
        lines.append("")
    lines.append("Unflagged is not proven complete. Missing is not zero. Baseline-only bins treat later-version annotations as missing. These development-exposed associations are not causal or confirmatory. FastOMA uses a supplied tree; sequence-only OrthoFinder is a group-clique diagnostic.")
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    for item in inputs:
        check(item)
    result = dict(status="descriptive_swiss_fragment_strata_exported", inputs=inputs, rows=rows,
        annotation_coverage=dict(matched=features["matched"], missing=features["missing"]),
        outputs=[record(p) for p in sorted(output.iterdir())],
        new_inferential_claims=False, publication_ready=False)
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha", required=True)
    parser.add_argument("--counts-sha256", default=COUNTS_SHA)
    args = parser.parse_args()
    export(args.counts.resolve(), args.admission.resolve(), args.admission_sha, args.output.absolute(), counts_sha=args.counts_sha256)
