#!/usr/bin/env python3
"""Score predictions on a predefined OrthoBench RefOG partition."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from typing import Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    read_clusters,
)


def read_named_groups(directory: Path, names: Sequence[str]) -> dict[str, set[str]]:
    groups = {}
    for name in names:
        path = directory / name
        if not path.is_file():
            raise ValueError(f"RefOG file not found: {path}")
        groups[name] = {
            line.strip() for line in path.read_text().splitlines() if line.strip()
        }
    return groups


def score_partition(
    predictions: Sequence[set[str]],
    references: Mapping[str, set[str]],
    uncertain: Mapping[str, set[str]],
) -> dict:
    """Reproduce the official equal-RefOG pair score for selected groups."""

    total_true_positive = 0.0
    total_false_positive = 0.0
    total_false_negative = 0.0
    exact = 0
    records = []
    for name, original_reference in references.items():
        excluded = uncertain.get(name, set())
        reference = original_reference - excluded
        false_positive = 0.0
        false_negative = 0.0
        true_positive = 0.0
        not_present = set(reference)
        splits = 0
        for prediction in predictions:
            overlap = len(reference & prediction)
            if not overlap:
                continue
            filtered_prediction = prediction - excluded
            overlap = len(reference & filtered_prediction)
            if not overlap:
                continue
            splits += 1
            not_present.difference_update(filtered_prediction)
            true_positive += overlap * (overlap - 1) / 2
            false_positive += overlap * (len(filtered_prediction) - overlap)
            false_negative += (len(reference) - overlap) * overlap
        false_negative += len(not_present) * (len(reference) - 1)
        false_negative /= 2
        if false_negative == 0 and false_positive == 0:
            exact += 1
        denominator = len(reference) - 1
        if denominator <= 0:
            raise ValueError(f"RefOG {name} has fewer than two scored genes")
        total_true_positive += true_positive / denominator
        total_false_positive += false_positive / denominator
        total_false_negative += false_negative / denominator
        records.append({
            "refog": name,
            "genes": len(reference),
            "true_positive": true_positive,
            "false_positive": false_positive,
            "false_negative": false_negative,
            "splits": splits,
            "exact": false_negative == 0 and false_positive == 0,
        })
    precision = total_true_positive / (
        total_true_positive + total_false_positive
    )
    recall = total_true_positive / (
        total_true_positive + total_false_negative
    )
    f_score = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "refogs": len(references),
        "f_score": 100 * f_score,
        "precision": 100 * precision,
        "recall": 100 * recall,
        "exact_refogs": exact,
        "weighted_true_positive": total_true_positive,
        "weighted_false_positive": total_false_positive,
        "weighted_false_negative": total_false_negative,
        "refog_records": records,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--partition-json", required=True, type=Path)
    parser.add_argument(
        "--partition", choices=("development", "validation", "all"), required=True
    )
    parser.add_argument("--json", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = json.loads(args.partition_json.read_text())
    if args.partition == "all":
        names = sorted(manifest["development"] + manifest["validation"])
    else:
        names = manifest[args.partition]
    references = read_named_groups(args.refogs, names)
    low_certainty = args.refogs / "low_certainty_assignments"
    uncertain = {
        name: (
            {
                line.strip()
                for line in (low_certainty / name).read_text().splitlines()
                if line.strip()
            }
            if (low_certainty / name).is_file()
            else set()
        )
        for name in names
    }
    predictions = read_clusters(args.predictions)
    result = score_partition(predictions, references, uncertain)
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": "Official-formula OrthoBench partition score",
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "inputs": {
            "predictions": file_provenance(args.predictions),
            "partition_manifest": file_provenance(args.partition_json),
            "refog_directory": str(args.refogs.resolve()),
        },
        "partition": args.partition,
        "refog_names": names,
        "score": result,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.json)
    print(
        f"{args.partition}: F={result['f_score']:.6f} "
        f"P={result['precision']:.6f} R={result['recall']:.6f} "
        f"exact={result['exact_refogs']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
