#!/usr/bin/env python3
"""Analyze label-blind satellite merge evidence after inference is frozen."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import sys
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    load_refogs,
)


FEATURES = (
    "iteration",
    "source_size",
    "target_size",
    "target_seed_families",
    "support",
    "margin",
    "species_overlap_fraction",
    "forward_hits",
    "reverse_hits",
    "forward_average",
    "reverse_average",
    "forward_maximum",
    "reverse_maximum",
    "forward_coverage",
    "reverse_coverage",
    "forward_normalized_support",
    "reverse_normalized_support",
)


def _percentile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def summarize_values(values: Sequence[float]) -> dict[str, float | int]:
    if not values:
        return {"count": 0}
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return {"count": len(values), "finite_count": 0}
    return {
        "count": len(values),
        "finite_count": len(finite_values),
        "minimum": min(finite_values),
        "p10": _percentile(finite_values, 0.10),
        "p25": _percentile(finite_values, 0.25),
        "median": statistics.median(finite_values),
        "mean": statistics.fmean(finite_values),
        "p75": _percentile(finite_values, 0.75),
        "p90": _percentile(finite_values, 0.90),
        "maximum": max(finite_values),
    }


def annotate_merge(
    merge: dict,
    gene_refogs: dict[str, set[int]],
) -> dict[str, int | str]:
    source = [str(gene) for gene in merge["source_genes"]]
    target = [str(gene) for gene in merge["target_genes"]]
    reference_genes = set(gene_refogs)
    same_refog_pairs = 0
    cross_refog_pairs = 0
    for left, right in itertools.product(source, target):
        left_refogs = gene_refogs.get(left, set())
        right_refogs = gene_refogs.get(right, set())
        if not left_refogs or not right_refogs:
            continue
        if left_refogs & right_refogs:
            same_refog_pairs += 1
        else:
            cross_refog_pairs += 1
    source_reference = len(set(source) & reference_genes)
    target_reference = len(set(target) & reference_genes)
    reference_nonreference_pairs = (
        source_reference * (len(target) - target_reference)
        + target_reference * (len(source) - source_reference)
    )
    if cross_refog_pairs:
        category = "cross_refog"
    elif same_refog_pairs:
        category = "same_refog"
    elif reference_nonreference_pairs:
        category = "reference_nonreference"
    else:
        category = "unlabelled"
    return {
        "category": category,
        "same_refog_pairs": same_refog_pairs,
        "cross_refog_pairs": cross_refog_pairs,
        "reference_nonreference_pairs": reference_nonreference_pairs,
        "source_reference_genes": source_reference,
        "target_reference_genes": target_reference,
    }


def analyze_merges(merges: Sequence[dict], refogs: Sequence[set[str]]) -> dict:
    gene_refogs: dict[str, set[int]] = {}
    for refog_id, refog in enumerate(refogs):
        for gene in refog:
            gene_refogs.setdefault(gene, set()).add(refog_id)

    categories: dict[str, list[dict]] = {}
    contributions = {
        "same_refog_pairs": 0,
        "cross_refog_pairs": 0,
        "reference_nonreference_pairs": 0,
    }
    for merge in merges:
        annotation = annotate_merge(merge, gene_refogs)
        categories.setdefault(str(annotation["category"]), []).append(merge)
        for field in contributions:
            contributions[field] += int(annotation[field])

    category_summaries = {}
    for category, records in sorted(categories.items()):
        category_summaries[category] = {
            "merges": len(records),
            "features": {
                feature: summarize_values(
                    [float(record[feature]) for record in records]
                )
                for feature in FEATURES
            },
        }
    return {
        "merges": len(merges),
        "reference_genes": len(gene_refogs),
        "refogs": len(refogs),
        "pair_contributions": contributions,
        "categories": category_summaries,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-json", required=True, type=Path)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    trace_payload = json.loads(args.trace_json.read_text())
    merges = trace_payload.get("merges")
    if not isinstance(merges, list):
        raise ValueError("trace JSON must contain a merges list")
    analysis = analyze_merges(merges, load_refogs(args.refogs))
    payload = {
        "schema_version": 1,
        "description": "Post-inference candidate merge evidence analysis",
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "inputs": {
            "trace": file_provenance(args.trace_json),
            "refog_directory": str(args.refogs.resolve()),
        },
        "analysis": analysis,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)
    print(json.dumps({
        "categories": {
            category: summary["merges"]
            for category, summary in analysis["categories"].items()
        },
        "pair_contributions": analysis["pair_contributions"],
    }, sort_keys=True))
    print(args.json.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
