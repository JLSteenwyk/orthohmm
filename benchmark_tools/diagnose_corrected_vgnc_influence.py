"""Exhaustive single-block deletion sensitivity for corrected VGNC scored tables."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

from benchmark_tools.audit_corrected_vgnc_components import SOURCE_SHA, read_cells, read_reference
from benchmark_tools.diagnose_vgnc_block_influence import CATEGORIES, deleted_scores, metrics
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

BASELINE = "orthofinder_3_1_5_full"


def removals(path, reference):
    _, totals = read_cells(path, reference)
    removed = {block: Counter() for block in reference}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts = {c: int(row[c]) for c in CATEGORIES}
            for block in {row["block_left"], row["block_right"]}:
                removed[block].update(counts)
    return totals, {b: {c: value[c] for c in CATEGORIES} for b, value in removed.items()}


def contrast(candidate, baseline, full_difference):
    if set(candidate) != set(baseline) or not candidate:
        raise ValueError("Require matched nonempty block inventories")
    values = {b: candidate[b] - baseline[b] for b in candidate}
    low = min(values, key=lambda b: (values[b], b))
    high = min(values, key=lambda b: (-values[b], b))
    return dict(full_difference=full_difference, minimum_deleted_difference=values[low],
                maximum_deleted_difference=values[high], minimum_block=low, maximum_block=high,
                positive=sum(v > 0 for v in values.values()), zero=sum(v == 0 for v in values.values()),
                negative=sum(v < 0 for v in values.values()))


def run(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source_ref = record(source)
    if source_ref["sha256"] != SOURCE_SHA:
        raise ValueError("Corrected mapping checksum differs")
    data = json.loads(source.read_text())
    refs = [source_ref, data["reference_table"], *data["checked_records"],
            *[m["table"] for m in data["methods"]], record(Path(__file__)),
            record(Path(__file__).with_name("audit_corrected_vgnc_components.py")),
            record(Path(__file__).with_name("diagnose_vgnc_block_influence.py"))]
    for ref in refs:
        check(ref)
    reference = read_reference(Path(data["reference_table"]["path"]))
    keys = [m["key"] for m in data["methods"]]
    if len(reference) != 16844 or len(keys) != 8 or len(set(keys)) != 8 or BASELINE not in keys:
        raise ValueError("Wrong corrected panel")
    output.mkdir(parents=True, exist_ok=False)
    table = output / "all_block_deletions.tsv.gz"
    methods, scores = {}, {}
    with gzip.open(table, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "block", *["removed_" + c for c in CATEGORIES],
                         "precision", "recall", "f1", "f1_change"])
        for method in data["methods"]:
            key = method["key"]
            totals, removed = removals(Path(method["table"]["path"]), reference)
            if totals != method["counts"]:
                raise ValueError("Counts differ from retained mapping")
            full = metrics(totals)
            scores[key] = {}
            changes = []
            for block in sorted(reference):
                _, value = deleted_scores(totals, removed[block])
                if any(v is None for v in value.values()):
                    raise ValueError("Undefined deletion score")
                scores[key][block] = value["f1"]
                delta = value["f1"] - full["f1"]
                changes.append(dict(block=block, f1_change=delta, removed=removed[block]))
                writer.writerow([key, block, *[removed[block][c] for c in CATEGORIES],
                                 value["precision"], value["recall"], value["f1"], delta])
            methods[key] = dict(label=method["label"], counts=totals, full_metrics=full,
                blocks=len(reference), minimum_f1_change=min(r["f1_change"] for r in changes),
                maximum_f1_change=max(r["f1_change"] for r in changes),
                largest_absolute_changes=sorted(changes, key=lambda r: (-abs(r["f1_change"]), r["block"]))[:10])
    comparisons = {key: contrast(scores[key], scores[BASELINE],
                    methods[key]["full_metrics"]["f1"] - methods[BASELINE]["full_metrics"]["f1"])
                   for key in keys if key != BASELINE}
    for ref in refs:
        check(ref)
    result = dict(status="corrected_vgnc_single_block_sensitivity_described", methods=methods,
                  baseline=BASELINE, comparisons=comparisons, checked_records=refs, table=record(table),
                  rows=len(keys)*len(reference), uncertainty_admitted=False, publication_ready=False,
                  limitations=[
                      "Exploratory fixed-table deletion sensitivity after previous score inspection, not a confidence interval.",
                      "Deletes each block and incident scored pairs once; does not rerun native eligibility, inference or scoring.",
                      "Separate deletions share rows and are dependent; sign stability does not establish joint-deletion or population robustness.",
                      "Original observed counts and scores remain unchanged; no method or endpoint selection."])
    with (output / "report.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.output)
