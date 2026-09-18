"""Exploratory deletion influence on fixed scored VGNC tables, without CIs."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_vgnc_family_dependencies import audit, blocks, reference_data
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

CATEGORIES = ("TP", "FP", "FN")
CONTRASTS = ((1, 0), (3, 2), (2, 0), (3, 1))


def metrics(counts):
    tp, fp, fn = (counts[c] for c in CATEGORIES)
    if any(type(v) is not int or v < 0 for v in (tp, fp, fn)):
        raise ValueError("Invalid confusion counts")
    return dict(precision=tp / (tp + fp) if tp + fp else None,
                recall=tp / (tp + fn) if tp + fn else None,
                f1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None)


def incidents(rows, mapping):
    totals = Counter()
    removed = {block: Counter() for block in set(mapping.values())}
    seen = set()
    for row in rows:
        a, b, category, left, right, _, _ = row
        if category not in CATEGORIES or left not in mapping or right not in mapping or a == b:
            raise ValueError("Invalid scored row")
        pair = tuple(sorted((a, b)))
        if pair in seen:
            raise ValueError("Duplicate or conflicting scored pair")
        seen.add(pair)
        totals[category] += 1
        # Distinct endpoints each receive one incident row, not two within a block.
        for block in {mapping[left], mapping[right]}:
            removed[block][category] += 1
    return ({c: totals[c] for c in CATEGORIES},
            {block: {c: counts[c] for c in CATEGORIES} for block, counts in removed.items()})


def deleted_scores(total, removed):
    remaining = {c: total[c] - removed[c] for c in CATEGORIES}
    return remaining, metrics(remaining)


def run(root, output):
    if output.exists():
        raise FileExistsError(output)
    dependency = audit(root)
    reference = next(r for r in dependency["checked_inputs"] if Path(r["path"]).name == "vgnc-orthologs.txt.gz")
    truth, _ = reference_data(Path(reference["path"]))
    mapping, _ = blocks(truth)
    output.mkdir(parents=True)
    stages, results = [], []
    table = output / "all_block_deletions.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "block", "removed_tp", "removed_fp", "removed_fn", "remaining_tp", "remaining_fp", "remaining_fn", "precision", "recall", "f1", "f1_change"])
        for stage in dependency["stages"]:
            with gzip.open(stage["raw"]["path"], "rt", newline="") as raw:
                total, removed = incidents(csv.reader(raw, delimiter="\t"), mapping)
            if total != stage["counts"]:
                raise ValueError("Scored inventory differs")
            full = metrics(total)
            scores, changes = {}, []
            for block in sorted(removed):
                remaining, value = deleted_scores(total, removed[block])
                if any(v is None for v in value.values()):
                    raise ValueError("Deletion produces undefined native ratios")
                delta = value["f1"] - full["f1"]
                scores[block] = value
                changes.append(dict(block=block, f1_change=delta, removed=removed[block]))
                writer.writerow([stage["stage_index"], block, *[removed[block][c] for c in CATEGORIES],
                                 *[remaining[c] for c in CATEGORIES], value["precision"], value["recall"], value["f1"], delta])
            stages.append(scores)
            results.append(dict(stage_index=stage["stage_index"], full_counts=total, full_metrics=full,
                blocks=len(scores), minimum_f1_change=min(r["f1_change"] for r in changes),
                maximum_f1_change=max(r["f1_change"] for r in changes),
                largest_absolute_f1_changes=sorted(changes, key=lambda r: (-abs(r["f1_change"]), r["block"]))[:10]))
    contrasts = []
    for on, off in CONTRASTS:
        values = [stages[on][block]["f1"] - stages[off][block]["f1"] for block in sorted(stages[on])]
        contrasts.append(dict(on=on, off=off, full_f1_difference=results[on]["full_metrics"]["f1"] - results[off]["full_metrics"]["f1"],
            minimum_deleted_f1_difference=min(values), maximum_deleted_f1_difference=max(values),
            positive=sum(v > 0 for v in values), zero=sum(v == 0 for v in values), negative=sum(v < 0 for v in values)))
    for item in dependency["checked_inputs"]:
        check(item)
    result = dict(status="fixed_scored_table_block_influence_described", source=record(__file__),
        helper=record(Path(__file__).with_name("audit_vgnc_family_dependencies.py")),
        inputs=dependency["checked_inputs"], reference=dependency["reference"], table=record(table),
        stages=results, paired_contrasts=contrasts, uncertainty_admitted=False, publication_ready=False,
        limitations=["Exploratory influence after prior score inspection; not confirmatory inference.",
            "Deletes incident scored rows only; does not rerun prediction, reference construction or native eligibility.",
            "Separate deletions share rows and are dependent; ranges are not confidence intervals.",
            "Single-block sign stability does not imply robustness to joint deletions or new datasets.",
            "Historical four-stage evidence only, not corrected-release or competitor results."])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute())
