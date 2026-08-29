#!/usr/bin/env python3
"""Validate reciprocal profile-HMM repair of synthetic split families."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np

from benchmark_tools.hmm_generalization_holdout import (
    _mutate_sequence,
    _score_pairs,
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
    build_reciprocal_profile_edges,
)


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4
LOCKED_THRESHOLD_RATIO = 0.7
LOCKED_MIN_SUPPORT = 2
DEVELOPMENT_RATIOS = (0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8)


@dataclass(frozen=True)
class SplitEvidence:
    gene_names: tuple[str, ...]
    clusters: tuple[tuple[int, ...], ...]
    gene_to_species: np.ndarray
    profile_hits: ProfileHits
    self_thresholds: dict[int, float]
    true_pairs: frozenset[tuple[int, int]]
    families: int


def build_split_evidence(seed: int, families: int, cpu: int) -> SplitEvidence:
    """Build two three-species profiles for every synthetic family."""
    dataset = generate_synthetic_dataset(seed, families=families)
    rng = np.random.default_rng(seed + 99173)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names: list[str] = []
    sequences: list[str] = []
    gene_to_species: list[int] = []
    clusters: list[tuple[int, ...]] = []
    profiles = {}

    for family_id, original_members in dataset.training_sequences.items():
        extra_member = _mutate_sequence(
            original_members[1],
            np.ones(len(original_members[1]), dtype=np.float64),
            0.18,
            rng,
        )
        members = original_members + (extra_member,)
        for member_indices in ((0, 1, 5), (2, 3, 4)):
            cluster_id = len(clusters)
            cluster = []
            profile_sequences = []
            profile_names = []
            for member_index in member_indices:
                gene_id = len(gene_names)
                gene_name = f"family_{family_id}_species_{member_index}"
                cluster.append(gene_id)
                gene_names.append(gene_name)
                sequences.append(members[member_index])
                gene_to_species.append(member_index)
                profile_sequences.append(members[member_index])
                profile_names.append(gene_name)
            clusters.append(tuple(cluster))
            profile = build_msa_profile(
                profile_sequences,
                profile_names,
                matrix,
                background,
            )
            if profile is None:
                raise RuntimeError(f"could not build split profile {cluster_id}")
            profiles[cluster_id] = profile

    database = species_sequences(gene_names, sequences, "split_families")
    pairs = np.asarray(
        [
            (profile_id, gene_id)
            for profile_id in range(len(profiles))
            for gene_id in range(len(gene_names))
        ],
        dtype=np.int32,
    )
    profile_ids, profile_lengths, scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    score_matrix = scores.reshape(len(profiles), len(gene_names))
    self_thresholds = {
        profile_id: float(
            np.min(score_matrix[profile_id, list(clusters[profile_id])])
        )
        for profile_id in range(len(profiles))
    }
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        scores,
        np.repeat(profile_lengths, len(gene_names)),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    )
    significant = evalues < EVALUE_THRESHOLD
    profile_hits = ProfileHits(
        profile_cluster_ids=profile_ids[pairs[significant, 0]],
        gene_ids=pairs[significant, 1].copy(),
        scores=scores[significant].astype(np.float64, copy=False),
        evalues=evalues[significant],
        candidate_count=len(pairs),
    )
    true_pairs = frozenset(
        (family_id * 2, family_id * 2 + 1)
        for family_id in range(families)
    )
    return SplitEvidence(
        tuple(gene_names),
        tuple(clusters),
        np.asarray(gene_to_species, dtype=np.int32),
        profile_hits,
        self_thresholds,
        true_pairs,
        families,
    )


def evaluate_evidence(
    evidence: SplitEvidence,
    threshold_ratio: float,
    min_support: int,
) -> dict:
    """Evaluate inferred cluster pairs without exposing final benchmarks."""
    edges, reported_pairs = build_reciprocal_profile_edges(
        evidence.gene_names,
        evidence.clusters,
        evidence.profile_hits,
        evidence.self_thresholds,
        evidence.gene_to_species,
        threshold_ratio=threshold_ratio,
        min_support=min_support,
    )
    cluster_by_gene = np.full(len(evidence.gene_names), -1, dtype=np.int32)
    for cluster_id, cluster in enumerate(evidence.clusters):
        cluster_by_gene[list(cluster)] = cluster_id
    predicted_pairs = {
        tuple(sorted((
            int(cluster_by_gene[int(source)]),
            int(cluster_by_gene[int(target)]),
        )))
        for source, target in zip(edges.sources, edges.targets)
    }
    if len(predicted_pairs) != reported_pairs:
        raise RuntimeError("edge-derived and reported profile pairs disagree")
    true_positives = len(predicted_pairs & evidence.true_pairs)
    false_positives = len(predicted_pairs - evidence.true_pairs)
    false_negatives = len(evidence.true_pairs - predicted_pairs)
    precision = (
        true_positives / len(predicted_pairs) if predicted_pairs else 0.0
    )
    recall = true_positives / len(evidence.true_pairs)
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    close_paralog_false_pairs = sum(
        1
        for left, right in predicted_pairs - evidence.true_pairs
        if (left // 2) // 2 == (right // 2) // 2
    )
    return {
        "threshold_ratio": threshold_ratio,
        "minimum_support_per_direction": min_support,
        "true_pairs": len(evidence.true_pairs),
        "predicted_pairs": len(predicted_pairs),
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "close_paralog_false_pairs": close_paralog_false_pairs,
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "profile_edges": len(edges),
    }


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    evidence = build_split_evidence(seed, families, cpu)
    return evaluate_evidence(
        evidence, LOCKED_THRESHOLD_RATIO, LOCKED_MIN_SUPPORT
    )


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

    development_evidence = build_split_evidence(
        args.development_seed, args.families, args.cpu
    )
    development_screen = [
        evaluate_evidence(
            development_evidence, threshold_ratio, LOCKED_MIN_SUPPORT
        )
        for threshold_ratio in DEVELOPMENT_RATIOS
    ]
    development = evaluate_evidence(
        development_evidence,
        LOCKED_THRESHOLD_RATIO,
        LOCKED_MIN_SUPPORT,
    )
    holdout = evaluate_split(args.holdout_seed, args.families, args.cpu)
    passed = all(
        result["precision"] == 1.0 and result["recall"] == 1.0
        for result in (development, holdout)
    )
    payload = {
        "schema_version": 1,
        "description": (
            "Synthetic development/holdout test of reciprocal profile-HMM "
            "edges for species-complementary split protein families."
        ),
        "label_policy": (
            "Uses generated family and species labels only; no OrthoBench, "
            "QfO, Three Kingdoms, OrthoFinder, or other tool output is read."
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
            "training_members_per_family": 6,
            "split_members_per_profile": 3,
            "cpu": args.cpu,
            "locked_threshold_ratio": LOCKED_THRESHOLD_RATIO,
            "locked_minimum_support_per_direction": LOCKED_MIN_SUPPORT,
            "require_disjoint_species": True,
        },
        "selection_policy": (
            "The 0.70 ratio was locked at the conservative end of the "
            "perfect-precision/perfect-recall development plateau before "
            "evaluating the holdout configuration."
        ),
        "development_parameter_screen": development_screen,
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
        "development_f_score": development["f_score"],
        "holdout_f_score": holdout["f_score"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
