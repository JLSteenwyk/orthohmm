"""Export prespecified descriptive SwissTrees identity strata, without new inference."""

import argparse
import csv
import json
from pathlib import Path

from benchmark_tools.export_swiss_descriptive_strata import build, COUNTS_SHA, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

ADMISSION_SHA = "cb04162af62fbd58fcfe8f02cbb78bc53bcf49ac20a487aae89911bc3a4e7b2d"
HELPER_SHA = "65ddc574fb312968b4266bb62476a8f312185104a8d85c75bcf423894e11f728"
REFERENCE = "orthofinder_3_1_5_full"


def rows_with_differences(counts, admission):
    if (admission["status"] != "corrected_swiss_identity_features_verified"
            or admission["prediction_statistics_evaluated"] is not False
            or admission["publication_ready"] is not False):
        raise ValueError("Require admitted unscored identity features")
    strata = admission["strata"]
    if set(strata) != {"lower_identity", "higher_identity", "missing_identity"}:
        raise ValueError("Wrong identity bins")
    flattened = [f for members in strata.values() for f in members]
    if len(set(flattened)) != len(flattened) or set(flattened) != set(admission["family_memberships"]):
        raise ValueError("Identity bins do not partition reference families")
    rows = build(counts, dict(family_memberships=admission["family_memberships"],
                            primary_strata=strata, secondary_strata={}))
    reference = {r["stratum"]: r for r in rows if r["method"] == REFERENCE}
    for row in rows:
        ref = reference[row["stratum"]]
        for metric in METRICS:
            row["delta_" + metric] = None if row[metric] is None or ref[metric] is None else row[metric] - ref[metric]
    return rows


def export(counts, admission, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    helper = record(Path(__file__).with_name("export_swiss_descriptive_strata.py"))
    if helper["sha256"] != HELPER_SHA:
        raise ValueError("Descriptive statistic implementation changed")
    inputs = [record(counts), record(admission), helper, record(__file__)]
    features = read_frozen(admission, ADMISSION_SHA)
    rows = rows_with_differences(read_frozen(counts, COUNTS_SHA), features)
    output.mkdir(parents=True)
    fields = ["method", "stratum", "families", "status", *METRICS,
              *["delta_" + k for k in METRICS], "prediction_semantics"]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in r.items()} for r in rows)
    lines = ["# Corrected SwissTrees Identity Strata", "",
        "Descriptive percentages; differences are percentage points versus full OrthoFinder.",
        "F1 is the harmonic mean of macro precision and recall. No additional intervals or significance claims.", ""]
    for name in dict.fromkeys(r["stratum"] for r in rows):
        selected = [r for r in rows if r["stratum"] == name]
        lines.extend([f"## {name} ({selected[0]['families']} families)", "",
            "| Method | F1 | Precision | Recall | F1 difference | Status |",
            "|---|---:|---:|---:|---:|---|"])
        for row in selected:
            values = ["NA" if row[k] is None else f"{100*row[k]:.3f}" for k in (*METRICS, "delta_F1")]
            lines.append("| " + " | ".join([row["method"], *values, row["status"]]) + " |")
        lines.append("")
    lines.append("Development-exposed, alignment-dependent identity bins are not calibrated evolutionary distances. Missing is not zero. FastOMA uses a supplied tree; sequence-only OrthoFinder is a group-clique diagnostic.")
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    for item in inputs:
        check(item)
    result = dict(status="descriptive_swiss_identity_strata_exported", inputs=inputs,
        median_family_identity=features["median_family_identity"], rows=rows,
        outputs=[record(p) for p in sorted(output.iterdir())],
        new_inferential_claims=False, publication_ready=False)
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    export(args.counts.resolve(), args.admission.resolve(), args.output.absolute())
