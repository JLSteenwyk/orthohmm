#!/usr/bin/env python3
"""Validate calibrated profile-profile repair on held-out split families."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np

from benchmark_tools.reciprocal_profile_merge_holdout import (
    build_split_evidence,
)
from orthohmm.search.matrices import get_background_freqs
from orthohmm.search.profile_profile import build_profile_profile_edges


MATRIX_NAME = "BLOSUM62"
LOCKED_SIMILARITY_THRESHOLD = 0.4
LOCKED_MAX_COMBINED_GENES = 80


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    """Evaluate one deterministic split without final-benchmark labels."""
    evidence = build_split_evidence(seed, families, cpu)
    result = build_profile_profile_edges(
        evidence.profiles,
        evidence.clusters,
        evidence.gene_names,
        get_background_freqs(MATRIX_NAME),
        cpu,
        similarity_threshold=LOCKED_SIMILARITY_THRESHOLD,
        max_combined_genes=LOCKED_MAX_COMBINED_GENES,
    )
    cluster_by_gene = np.full(len(evidence.gene_names), -1, dtype=np.int32)
    for cluster_id, cluster in enumerate(evidence.clusters):
        cluster_by_gene[list(cluster)] = cluster_id
    predicted_pairs = {
        tuple(sorted((
            int(cluster_by_gene[int(source)]),
            int(cluster_by_gene[int(target)]),
        )))
        for source, target in zip(result.edges.sources, result.edges.targets)
    }
    true_positives = len(predicted_pairs & evidence.true_pairs)
    false_positives = len(predicted_pairs - evidence.true_pairs)
    false_negatives = len(evidence.true_pairs - predicted_pairs)
    precision = true_positives / len(predicted_pairs) if predicted_pairs else 0.0
    recall = true_positives / len(evidence.true_pairs)
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "seed": seed,
        "families": families,
        "profiles": len(evidence.profiles),
        "candidate_pairs": result.candidate_pairs,
        "predicted_pairs": len(predicted_pairs),
        "reported_mutual_pairs": result.reciprocal_pairs,
        "profile_edges": len(result.edges),
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "passed": precision == 1.0 and recall == 1.0,
    }


def _git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args], check=False, capture_output=True, text=True
    ).stdout.strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-seed", type=int, default=1103)
    parser.add_argument("--holdout-seed", type=int, default=2909)
    parser.add_argument("--families", type=int, default=36)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    development = evaluate_split(
        args.development_seed, args.families, args.cpu
    )
    holdout = evaluate_split(args.holdout_seed, args.families, args.cpu)
    passed = development["passed"] and holdout["passed"]
    payload = {
        "schema_version": 1,
        "description": (
            "Synthetic development/holdout test of self-normalized, "
            "mutual-nearest profile-profile alignment for split families."
        ),
        "label_policy": (
            "Uses generated family labels only; no OrthoBench, QfO, Three "
            "Kingdoms, species identities, or competing-tool output is read."
        ),
        "selection_policy": (
            "The 0.40 threshold was locked below the minimum true-family "
            "development similarity (0.6676) and above the maximum close-"
            "paralog similarity (0.3735), maximizing sensitivity within the "
            "perfect-separation interval. An exact-4-mer prefilter capacity "
            "of 30 was the smallest tested development setting with complete "
            "true-pair recall and was then applied unchanged to holdout."
        ),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git": {
            "commit": _git_output("rev-parse", "HEAD"),
            "dirty": bool(_git_output("status", "--short")),
        },
        "configuration": {
            "matrix": MATRIX_NAME,
            "similarity_threshold": LOCKED_SIMILARITY_THRESHOLD,
            "max_combined_genes": LOCKED_MAX_COMBINED_GENES,
            "development_seed": args.development_seed,
            "holdout_seed": args.holdout_seed,
            "families": args.families,
            "cpu": args.cpu,
            "prefilter_k": 4,
            "max_candidates_per_profile": 30,
        },
        "development": development,
        "holdout": holdout,
        "passed": passed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "passed": passed,
        "development_pairs": development["predicted_pairs"],
        "holdout_pairs": holdout["predicted_pairs"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
