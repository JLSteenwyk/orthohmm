#!/usr/bin/env python3
"""Validate direct HMM fallback without final benchmark labels."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np

from benchmark_tools.hmm_generalization_holdout import (
    _score_pairs,
    classification_metrics,
    generate_synthetic_dataset,
    species_sequences,
)
from orthohmm.search.evalue import batch_estimate_evalues
from orthohmm.search.matrices import (
    get_background_freqs,
    get_ka_params,
    get_matrix,
)
from orthohmm.search.msa_profile import build_msa_profile
from orthohmm.search.profile_expansion import (
    ProfileHits,
    build_direct_profile_fallback_edges,
    compute_profile_self_thresholds,
)


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    dataset = generate_synthetic_dataset(seed, families=families)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names: list[str] = []
    sequences: list[str] = []
    family_genes: dict[int, list[int]] = {}

    for family_id, members in dataset.training_sequences.items():
        family_genes[family_id] = []
        for member_index, sequence in enumerate(members):
            family_genes[family_id].append(len(gene_names))
            gene_names.append(f"family_{family_id}_member_{member_index}")
            sequences.append(sequence)

    query_start = len(gene_names)
    query_gene_ids = []
    for query in dataset.queries:
        query_gene_ids.append(len(gene_names))
        gene_names.append(query.query_id)
        sequences.append(query.sequence)
    database = species_sequences(gene_names, sequences, "direct_fallback")

    profiles = {}
    clusters = []
    training_gene_to_family = {}
    for family_id in range(families):
        members = family_genes[family_id]
        clusters.append(members)
        for gene in members:
            training_gene_to_family[gene] = family_id
        profile = build_msa_profile(
            [sequences[gene] for gene in members],
            [gene_names[gene] for gene in members],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build family profile {family_id}")
        profiles[family_id] = profile
    clusters.extend([[gene] for gene in query_gene_ids])

    thresholds = compute_profile_self_thresholds(
        profiles,
        gene_names,
        database,
        cpu,
        matrix_name=MATRIX_NAME,
        calibrate_weakest_member=True,
        calibration_profile_ids=set(profiles),
    )
    pairs = np.asarray([
        (profile_id, gene)
        for profile_id in range(families)
        for gene in query_gene_ids
    ], dtype=np.int32)
    _profile_ids, profile_lengths, scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        scores,
        np.repeat(profile_lengths, len(query_gene_ids)),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    )
    significant = evalues < EVALUE_THRESHOLD
    significant_pairs = pairs[significant]
    profile_hits = ProfileHits(
        profile_cluster_ids=significant_pairs[:, 0],
        gene_ids=significant_pairs[:, 1],
        scores=scores[significant],
        evalues=evalues[significant],
        candidate_count=len(pairs),
    )
    empty_i = np.asarray([], dtype=np.int32)
    empty_f = np.asarray([], dtype=np.float64)
    edges = build_direct_profile_fallback_edges(
        gene_names,
        clusters,
        profile_hits,
        thresholds,
        empty_i,
        empty_i,
        empty_f,
        set(profiles),
    )

    predictions: list[int | None] = [None] * len(query_gene_ids)
    for source, target in zip(edges.sources, edges.targets):
        source = int(source)
        target = int(target)
        if source >= query_start:
            query_gene, member_gene = source, target
        else:
            query_gene, member_gene = target, source
        predictions[query_gene - query_start] = training_gene_to_family[
            member_gene
        ]
    metrics = classification_metrics(dataset.queries, predictions)
    weights = np.asarray(edges.weights, dtype=np.float64)
    passed = bool(
        metrics["precision"] == 1.0
        and metrics["negative_rejection_rate"] == 1.0
        and metrics["recall"] > 0.0
        and len(edges) == metrics["predictions_made"]
        and (not len(weights) or (weights.min() >= 1.0 and weights.max() <= 5.0))
    )
    return {
        "seed": seed,
        "families": families,
        "queries": len(dataset.queries),
        "significant_profile_hits": int(significant.sum()),
        "historical_pairwise_unanchored_edges": 0,
        "direct_profile_fallback_edges": len(edges),
        "edge_weight_minimum": float(weights.min()) if len(weights) else None,
        "edge_weight_maximum": float(weights.max()) if len(weights) else None,
        "metrics": metrics,
        "passed": passed,
    }


def _git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args], check=False, capture_output=True, text=True
    ).stdout.strip()


def _sha256(path: Path) -> str:
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

    development = evaluate_split(args.development_seed, args.families, args.cpu)
    holdout = evaluate_split(args.holdout_seed, args.families, args.cpu)
    passed = development["passed"] and holdout["passed"]
    payload = {
        "schema_version": 1,
        "description": (
            "Synthetic development/holdout evaluation of direct profile-HMM "
            "edges for calibrated winners without pairwise anchors."
        ),
        "label_policy": (
            "Uses generated family labels only; no OrthoBench, QfO, Three "
            "Kingdoms, or competing-tool output is read."
        ),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git": {
            "commit": _git_output("rev-parse", "HEAD"),
            "dirty": bool(_git_output("status", "--short")),
        },
        "configuration": {
            "matrix": MATRIX_NAME,
            "evalue_threshold": EVALUE_THRESHOLD,
            "development_seed": args.development_seed,
            "holdout_seed": args.holdout_seed,
            "families": args.families,
            "cpu": args.cpu,
        },
        "development": development,
        "holdout": holdout,
        "passed": passed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "sha256": _sha256(args.output),
        "passed": passed,
        "development_edges": development["direct_profile_fallback_edges"],
        "holdout_edges": holdout["direct_profile_fallback_edges"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
