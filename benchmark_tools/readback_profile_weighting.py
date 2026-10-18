"""Independently recount saved synthetic weighting scores and generated truth."""

import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


MODES = ("none", "henikoff", "henikoff_cluster", "species_henikoff_cluster")


def recount(raw_path, generated):
    clusters, queries = generated["clusters"], generated["queries"]
    n_train = len(generated["names"])
    if sorted(gene for group in clusters for gene in group) != list(range(n_train)):
        raise ValueError("Generated training ownership differs")
    with np.load(raw_path, allow_pickle=False) as raw:
        scores, thresholds, evalues = (raw[name] for name in ("scores", "thresholds", "evalues"))
        if (scores.shape != (len(clusters), n_train + len(queries))
                or evalues.shape != (len(clusters), len(queries))
                or thresholds.shape != (len(clusters),)
                or not np.isfinite(scores).all() or np.isnan(evalues).any()):
            raise ValueError("Raw score dimensions or values differ")
        observed = np.array([min(scores[index, group]) for index, group in enumerate(clusters)])
        if not np.array_equal(observed, thresholds):
            raise ValueError("Raw thresholds differ from actual training minima")
        positive_count = negative_count = true_assigned = assigned = negative_rejected = top1 = 0
        correct_total, margins, categories = 0, [], {}
        for column, query in enumerate(queries):
            values = scores[:, n_train + column]
            eligible = [index for index in range(len(clusters))
                        if values[index] >= thresholds[index] and evalues[index, column] < 1e-4]
            prediction = max(eligible, key=lambda index: (values[index], -index)) if eligible else None
            predicted = prediction is not None
            truth = query["family_id"]
            positive = truth is not None
            correct = prediction == truth if positive else not predicted
            assigned += int(predicted)
            correct_total += int(correct)
            if positive:
                positive_count += 1
                true_assigned += int(prediction == truth)
                top1 += int(max(range(len(clusters)), key=lambda index: (values[index], -index)) == truth)
                margins.append(float(values[truth] - thresholds[truth]))
            else:
                negative_count += 1
                negative_rejected += int(not predicted)
            category = categories.setdefault(query["category"], {"correct": 0, "assigned": 0, "queries": 0})
            category["correct"] += int(correct)
            category["assigned"] += int(predicted)
            category["queries"] += 1
        for category in categories.values():
            category["accuracy"] = category["correct"] / category["queries"]
        precision = true_assigned / assigned if assigned else 0.
        recall = true_assigned / positive_count if positive_count else 0.
        return {
            "top1_positive_family_correct": top1,
            "strict_assignment": {
                "true_family_assignments": true_assigned, "positive_queries": positive_count,
                "negative_queries": negative_count, "predictions_made": assigned,
                "precision": precision, "recall": recall,
                "f_score": 2 * precision * recall / (precision + recall) if precision + recall else 0.,
                "negative_rejection_rate": negative_rejected / negative_count if negative_count else 0.,
                "overall_accuracy": correct_total / len(queries) if queries else 0.,
                "by_category": categories,
            },
            "own_family_threshold_margin": {
                "minimum": min(margins), "median": float(np.median(margins)), "maximum": max(margins),
            },
        }


def readback(summary):
    report = json.loads(summary.read_text())
    if report["status"] != "completed_synthetic_diagnostic_not_biological_validation":
        raise ValueError("Require completed producer report")
    checked = [record(summary), *report["sources"], *report["runtime"]["files"]]
    rows, cases = [], set()
    for result in report["results"]:
        check(result["generated"])
        checked.append(result["generated"])
        generated = json.loads(Path(result["generated"]["path"]).read_text())
        if [row["mode"] for row in result["rows"]] != list(MODES):
            raise ValueError("Incomplete or reordered weighting arms")
        first = result["rows"][0]
        key = first["seed"], first["extra_copies_per_family"]
        if key in cases:
            raise ValueError("Duplicate diagnostic case")
        cases.add(key)
        if len(generated["clusters"]) != 36 or len(generated["queries"]) != 180:
            raise ValueError("Wrong prespecified family/query inventory")
        for row in result["rows"]:
            if (row["seed"], row["extra_copies_per_family"]) != key:
                raise ValueError("Mixed case identity")
            check(row["raw"])
            checked.append(row["raw"])
            observed = recount(row["raw"]["path"], generated)
            if any(row[field] != value for field, value in observed.items()):
                raise ValueError("Saved synthetic result differs from raw recount")
            rows.append({"seed": row["seed"], "extra_copies": row["extra_copies_per_family"],
                "mode": row["mode"], **observed,
                "effective_count_min": min(row["effective_sequence_counts"]),
                "effective_count_max": max(row["effective_sequence_counts"])})
    if cases != {(seed, copies) for seed in (18102026, 19102026) for copies in (0, 4)}:
        raise ValueError("Wrong prespecified seeds or expansion cases")
    for pin in checked:
        check(pin)
    return {"status": "raw_synthetic_scores_thresholds_and_assignments_independently_recounted",
            "producer": record(summary), "checked_files": len(checked), "rows": rows,
            "reader_source": record(Path(__file__).resolve()),
            "reader_revision": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1], text=True).strip(),
            "effective_counts_scope": "Producer-reported effective counts; raw score/threshold/assignment statistics are independently recounted, not a second implementation of profile construction.",
            "scope": report["scope"], "biological_accuracy_evaluated": False,
            "next_action": report["next_action"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError("Preserve previous readbacks")
    report = readback(args.summary)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(record(args.output)))


if __name__ == "__main__":
    main()
