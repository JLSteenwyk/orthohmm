"""Descriptive all-method SwissTrees strata, using frozen corrected counts and bins."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

COUNTS_SHA = "121cc8adf3cd63878f19006ec5500a13f879d042eccd565e1ae5f5a863f43fb1"
STRATA_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
METRICS = ("F1", "PPV", "TPR")


def scores(values):
    p = math.fsum(v[0] for v in values) / len(values)
    r = math.fsum(v[1] for v in values) / len(values)
    return dict(F1=2*p*r/(p+r), PPV=p, TPR=r)


def build(report, strata):
    if (report["status"] != "corrected_swiss_comparison_intervals_audited"
            or report["scientific_inputs_admitted"] is not True):
        raise ValueError("Require admitted corrected counts")
    families = report["families"]
    if len(families) != 18 or len(set(families)) != 18 or set(strata["family_memberships"]) != set(families):
        raise ValueError("Family universe differs")
    bins = {"all": families, **strata["primary_strata"], **strata["secondary_strata"]}
    for members in bins.values():
        if len(set(members)) != len(members) or not set(members) <= set(families):
            raise ValueError("Invalid bin membership")
    rows = []
    methods = report["reconstructed_counts"]["methods"]
    if len(methods) != 8 or {m["method"] for m in methods} != set(report["point_estimates"]):
        raise ValueError("Method inventory differs")
    for method in methods:
        name, status = method["method"], method["status"]
        values = {}
        if status == "counts_verified":
            if [f["family"] for f in method["families"]] != families:
                raise ValueError("Method family order/coverage differs")
            for family in method["families"]:
                raw = family["counts_without_prior"]
                if set(raw) != {"TP", "FP", "FN", "TN"} or any(type(v) is not int or v < 0 for v in raw.values()):
                    raise ValueError("Invalid family counts")
                tp, fp, fn = (raw[k]/2 + 1 for k in ("TP", "FP", "FN"))
                value = (tp/(tp+fp), tp/(tp+fn))
                values[family["family"]] = value
                if any(not math.isclose(scores([value])[k], family["statistics_with_prior"][k], rel_tol=0, abs_tol=1e-12) for k in METRICS):
                    raise ValueError("Family statistic does not reproduce")
            aggregate = scores(list(values.values()))
            if any(not math.isclose(aggregate[k], report["point_estimates"][name][k], rel_tol=0, abs_tol=1e-12) for k in METRICS):
                raise ValueError("Full-family statistic does not reproduce")
        elif status != "not_admitted" or report["point_estimates"][name] is not None:
            raise ValueError("Unknown status or imputed missing method")
        for bin_name, members in bins.items():
            available = bool(members) and bool(values)
            rows.append(dict(method=name, stratum=bin_name, families=len(members),
                family_members=members, status="descriptive" if available else
                "method_not_admitted" if not values else "empty_bin",
                prediction_semantics=method.get("prediction_semantics"),
                **(scores([values[f] for f in members]) if available else dict.fromkeys(METRICS))))
    return rows


def export(counts, strata, output, *, counts_sha=COUNTS_SHA):
    if output.exists():
        raise FileExistsError(output)
    inputs = [record(counts), record(strata)]
    rows = build(read_frozen(counts, counts_sha), read_frozen(strata, STRATA_SHA))
    output.mkdir(parents=True)
    with (output / "scores.tsv").open("x") as handle:
        fields = ["method", "stratum", "families", "status", *METRICS, "prediction_semantics"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    with (output / "scores.md").open("x") as handle:
        handle.write("# Descriptive Corrected SwissTrees Strata\n\n")
        handle.write("F1 is the harmonic mean of macro precision and recall, not mean family F1. Scores are percentages. Missing is not zero.\n\n")
        for bin_name in dict.fromkeys(r["stratum"] for r in rows):
            selected = [r for r in rows if r["stratum"] == bin_name]
            handle.write(f"## {bin_name} ({selected[0]['families']} families)\n\n| Method | F1 | Precision | Recall | Status |\n|---|---:|---:|---:|---|\n")
            for row in selected:
                values = ["NA" if row[k] is None else f"{100*row[k]:.3f}" for k in METRICS]
                handle.write("| " + " | ".join([row["method"], *values, row["status"]]) + " |\n")
            handle.write("\n")
        handle.write("Development-exposed descriptive bins; no additional intervals or significance claims. Relative shortness does not establish fragmentation. FastOMA uses a supplied OrthoFinder tree; sequence-only OrthoFinder is a group-clique diagnostic.\n")
    for item in inputs:
        check(item)
    result = dict(status="descriptive_swiss_strata_exported", inputs=inputs, source=record(__file__),
        rows=rows, outputs=[record(p) for p in sorted(output.iterdir())],
        new_inferential_claims=False, publication_ready=False)
    with (output / "manifest.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "strata", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--counts-sha256", default=COUNTS_SHA)
    args = parser.parse_args()
    export(args.counts.resolve(), args.strata.resolve(), args.output.absolute(), counts_sha=args.counts_sha256)
