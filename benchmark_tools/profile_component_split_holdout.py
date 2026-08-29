#!/usr/bin/env python3
"""Validate HMM-supported component splitting without benchmark labels."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Sequence

import numpy as np

from benchmark_tools.hmm_generalization_holdout import (
    _score_pairs,
    generate_synthetic_dataset,
    species_sequences,
)
from orthohmm.refinement import refine_cluster_indices
from orthohmm.search.evalue import batch_estimate_evalues
from orthohmm.search.matrices import (
    get_background_freqs,
    get_ka_params,
    get_matrix,
)
from orthohmm.search.msa_profile import build_msa_profile
from orthohmm.search.profile_expansion import compute_profile_self_thresholds


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4
COMPONENT_EDGE_THRESHOLD = 1.5


def _pairwise_metrics(labels: Sequence[int], clusters) -> dict[str, float | int]:
    truth = {
        (left, right)
        for left, right in combinations(range(len(labels)), 2)
        if labels[left] == labels[right]
    }
    predicted = set()
    for cluster in clusters:
        predicted.update(combinations(sorted(int(gene) for gene in cluster), 2))
    true_positive = len(truth & predicted)
    precision = true_positive / len(predicted) if predicted else 0.0
    recall = true_positive / len(truth) if truth else 0.0
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "true_pairs": len(truth),
        "predicted_pairs": len(predicted),
        "true_positive_pairs": true_positive,
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "clusters": len(clusters),
    }


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    dataset = generate_synthetic_dataset(seed, families=families)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names: list[str] = []
    sequences: list[str] = []
    labels: list[int] = []
    family_genes: dict[int, list[int]] = {}

    for family_id, members in dataset.training_sequences.items():
        family_genes[family_id] = []
        for member_index, sequence in enumerate(members):
            family_genes[family_id].append(len(gene_names))
            gene_names.append(f"family_{family_id}_member_{member_index}")
            sequences.append(sequence)
            labels.append(family_id)

    query_gene_ids = []
    next_negative_label = families
    for query in dataset.queries:
        query_gene_ids.append(len(gene_names))
        gene_names.append(query.query_id)
        sequences.append(query.sequence)
        if query.family_id is None:
            labels.append(next_negative_label)
            next_negative_label += 1
        else:
            labels.append(int(query.family_id))

    database = species_sequences(gene_names, sequences, "component_split")
    profiles = {}
    for family_id in range(families):
        members = family_genes[family_id]
        profile = build_msa_profile(
            [sequences[gene] for gene in members],
            [gene_names[gene] for gene in members],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build family profile {family_id}")
        profiles[family_id] = profile

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
    score_matrix = scores.reshape(families, len(query_gene_ids))
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        scores,
        np.repeat(profile_lengths, len(query_gene_ids)),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    ).reshape(score_matrix.shape)

    predictions: list[int | None] = []
    for query_index in range(len(query_gene_ids)):
        eligible = [
            profile_id
            for profile_id in range(families)
            if evalues[profile_id, query_index] < EVALUE_THRESHOLD
            and score_matrix[profile_id, query_index] >= thresholds[profile_id]
        ]
        predictions.append(
            max(
                eligible,
                key=lambda profile_id: (
                    score_matrix[profile_id, query_index], -profile_id
                ),
            )
            if eligible
            else None
        )

    parent_clusters = []
    gene_to_species = np.empty(len(gene_names), dtype=np.int32)
    query_cursor = 0
    for pair_index in range(families // 2):
        left_family = pair_index * 2
        right_family = left_family + 1
        query_count = 10
        query_members = query_gene_ids[query_cursor:query_cursor + query_count]
        query_cursor += query_count
        cluster = (
            family_genes[left_family]
            + family_genes[right_family]
            + query_members
        )
        parent_clusters.append(cluster)
        gene_to_species[np.asarray(cluster, dtype=np.int32)] = pair_index

    edge_sources = []
    edge_targets = []
    edge_weights = []
    for family_id, members in family_genes.items():
        for source, target in combinations(members, 2):
            edge_sources.append(source)
            edge_targets.append(target)
            edge_weights.append(2.0)
        paired_family = family_id ^ 1
        if family_id < paired_family:
            for source, target in zip(members, family_genes[paired_family]):
                edge_sources.append(source)
                edge_targets.append(target)
                edge_weights.append(1.0)

    profile_anchor_edges = 0
    negative_profile_edges = 0
    correct_profile_edges = 0
    for query, query_gene, prediction in zip(
        dataset.queries, query_gene_ids, predictions
    ):
        if prediction is None:
            continue
        edge_sources.append(query_gene)
        edge_targets.append(family_genes[prediction][0])
        edge_weights.append(2.0)
        profile_anchor_edges += 1
        if query.family_id is None:
            negative_profile_edges += 1
        elif prediction == query.family_id:
            correct_profile_edges += 1

    empty_i = np.asarray([], dtype=np.int32)
    empty_f = np.asarray([], dtype=np.float32)
    edge_sources_array = np.asarray(edge_sources, dtype=np.int32)
    edge_targets_array = np.asarray(edge_targets, dtype=np.int32)
    edge_weights_array = np.asarray(edge_weights, dtype=np.float32)
    common = dict(
        clusters=parent_clusters,
        hit_queries=empty_i,
        hit_targets=empty_i,
        hit_scores=empty_f,
        rbnh_queries=edge_sources_array,
        rbnh_targets=edge_targets_array,
        gene_to_species=gene_to_species,
        rbnh_scores=edge_weights_array,
    )
    baseline = refine_cluster_indices(**common)
    candidate = refine_cluster_indices(
        **common,
        component_split_high_duplication=True,
    )
    baseline_metrics = _pairwise_metrics(labels, baseline)
    candidate_metrics = _pairwise_metrics(labels, candidate)
    passed = (
        negative_profile_edges == 0
        and correct_profile_edges == profile_anchor_edges
        and candidate_metrics["precision"] == 1.0
        and candidate_metrics["f_score"] > baseline_metrics["f_score"]
        and candidate_metrics["recall"] > baseline_metrics["recall"]
    )
    return {
        "seed": seed,
        "families": families,
        "genes": len(gene_names),
        "parent_clusters": len(parent_clusters),
        "profile_anchor_edges": profile_anchor_edges,
        "correct_profile_edges": correct_profile_edges,
        "negative_profile_edges": negative_profile_edges,
        "weak_cross_family_edges": families // 2 * 5,
        "baseline_degree_split": baseline_metrics,
        "hmm_component_split": candidate_metrics,
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
            "Synthetic profile-HMM anchor and paralog-component split "
            "development/holdout evaluation."
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
            "component_edge_threshold": COMPONENT_EDGE_THRESHOLD,
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
        "development_f_score": development["hmm_component_split"]["f_score"],
        "holdout_f_score": holdout["hmm_component_split"]["f_score"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
