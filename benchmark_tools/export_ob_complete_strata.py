"""Describe all eight OrthoBench methods in the unchanged fourteen feature bins."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from benchmark_tools.assemble_ob_stratified_errors import validate_strata, STRATA_SHA
from benchmark_tools.bootstrap_orthobench import weighted_records, statistics, METRICS
from benchmark_tools.prepare_ob_error_strata import CATEGORIES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

PINS = {
    "ob_error_strata_prepared_20260916.json": STRATA_SHA,
    "orthobench_paired_uncertainty_20260916.json": "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9",
    "retained_ob_comparator_readback_20260926.json": "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55",
    "ob_stratified_error_results_20260916.json": "a8a072e9a3917fa7dabf25dba4affb8c3fd6b30c02286dda3efe882baceb3cac",
}


def summarize(scores, strata):
    indexed, names, sizes = {}, None, None
    for method, score in scores.items():
        current_names, current_sizes, _ = weighted_records(score["refog_records"])
        if names is not None and (names != current_names or not np.array_equal(sizes, current_sizes)):
            raise ValueError("Method reference names or sizes differ")
        names, sizes = current_names, current_sizes
        indexed[method] = {r["refog"]: r for r in score["refog_records"]}
    if not indexed or set(strata) != {d + ":" + c for d, cs in CATEGORIES.items() for c in cs}:
        raise ValueError("Require scores and all frozen strata")
    for dimension, categories in CATEGORIES.items():
        if sorted(n for category in categories for n in strata[dimension + ":" + category]) != names:
            raise ValueError("Dimension must partition every family once")
    rows = []
    for label, selected in sorted(strata.items()):
        for method in scores:
            counts = weighted_records([indexed[method][n] for n in sorted(selected)])[2].sum(axis=0) if selected else None
            rows.append(dict(stratum=label, method=method, families=sorted(selected), family_count=len(selected),
                metrics_percent=dict(zip(METRICS, statistics(counts).tolist())) if selected else None,
                weighted_counts=dict(zip(("tp", "fp", "fn"), counts.tolist())) if selected else None,
                status="descriptive" if selected else "empty_nonestimable"))
    return rows


def run(root, output, protocol_sha):
    if output.exists():
        raise FileExistsError(output)
    base = root / "benchmark_tools/results"
    protocol = record(base / "OB_COMPLETE_STRATA_PROTOCOL_20260928.md")
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Protocol changed")
    inputs, refs = {}, [protocol, record(Path(__file__))]
    for name, sha in PINS.items():
        ref = record(base / name)
        if ref["sha256"] != sha:
            raise ValueError("Changed frozen source: " + name)
        inputs[name] = json.loads((base / name).read_text())
        refs.append(ref)
    frozen, main, others, prior = [inputs[name] for name in PINS]
    validate_strata(frozen)
    if not others["all_scores_agree"]:
        raise ValueError("Comparator readback not admitted")
    refs.extend(others["checked_records"])
    refs.extend([*main["inputs"]["predictions"].values(), *main["inputs"]["references"], *main["inputs"]["uncertain"]])
    refs.extend(record(Path(__file__).with_name(name)) for name in (
        "bootstrap_orthobench.py", "assemble_ob_stratified_errors.py", "prepare_ob_error_strata.py"))
    for ref in refs:
        check(ref)
    scores = dict(main["scores"])
    for row in others["rows"]:
        if row["key"] in scores or row["agrees_with_retained_score"] is not True:
            raise ValueError("Duplicate or unverified comparator")
        scores[row["key"]] = row["score"]
    if len(scores) != 8:
        raise ValueError("Require all eight retained methods")
    rows = summarize(scores, frozen["strata"])
    matches = 0
    for row in rows:
        if row["method"] in main["scores"]:
            previous = prior["strata"][row["stratum"]]["point_estimates_percent"][row["method"]]
            current = row["metrics_percent"]
            if (previous is None) != (current is None) or (current is not None and any(abs(current[k] - previous[k]) > 1e-10 for k in METRICS)):
                raise ValueError("Original three-method point estimate changed")
            matches += 1
    for ref in refs:
        check(ref)
    output.mkdir(parents=True)
    result = dict(status="all_method_orthobench_descriptive_strata", checked_records=refs, rows=rows,
        original_method_rows_reproduced=matches, new_intervals_calculated=False, publication_ready=False,
        limitations=["Development-exposed descriptive extension; no new independent confirmation or causal inference.",
                     "Existing three-method confidence intervals and 84-endpoint adjustment are not extended to new methods.",
                     "Full-reference sufficient statistics and low-certainty conventions precede stratum restriction.",
                     "No new native inference, official-scorer invocation or proof of historical input consumption.",
                     "Family descriptors are not validated fragmentation, domain or duplication-history annotations."])
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with (output / "scores.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stratum", "method", "families", *METRICS, "status"])
        for row in rows:
            writer.writerow([row["stratum"], row["method"], row["family_count"],
                *[row["metrics_percent"][k] if row["metrics_percent"] is not None else "NA" for k in METRICS], row["status"]])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.output.absolute(), args.protocol_sha256)
    print(json.dumps(dict(rows=len(result["rows"]), original_rows_reproduced=result["original_method_rows_reproduced"])))
