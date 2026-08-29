#!/usr/bin/env python3
"""Validate single-copy gating for jackknife profile-HMM calibration."""

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
    compute_profile_self_thresholds,
    select_single_copy_profile_ids,
)


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4


def _predictions(
    query_scores: np.ndarray,
    evalues: np.ndarray,
    thresholds: dict[int, float],
) -> list[int | None]:
    predictions = []
    for query_index in range(query_scores.shape[1]):
        eligible = [
            profile_id
            for profile_id in range(query_scores.shape[0])
            if evalues[profile_id, query_index] < EVALUE_THRESHOLD
            and query_scores[profile_id, query_index]
            >= thresholds[profile_id]
        ]
        predictions.append(
            max(
                eligible,
                key=lambda profile_id: (
                    query_scores[profile_id, query_index], -profile_id
                ),
            )
            if eligible
            else None
        )
    return predictions


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    """Compare global and single-copy jackknife calibration on one split."""
    dataset = generate_synthetic_dataset(seed, families=families)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names: list[str] = []
    sequences: list[str] = []
    gene_to_species: list[int] = []
    family_genes: dict[int, list[int]] = {}

    for family_id, members in dataset.training_sequences.items():
        family_genes[family_id] = []
        for species_id, sequence in enumerate(members):
            family_genes[family_id].append(len(gene_names))
            gene_names.append(f"family_{family_id}_species_{species_id}")
            sequences.append(sequence)
            gene_to_species.append(species_id)

    query_start = len(gene_names)
    for query in dataset.queries:
        gene_names.append(query.query_id)
        sequences.append(query.sequence)
        gene_to_species.append(-1)
    database = species_sequences(gene_names, sequences, "jackknife_gate")

    profiles = {}
    clusters: list[list[int]] = []
    for family_id in range(families):
        cluster = family_genes[family_id]
        clusters.append(cluster)
        profile = build_msa_profile(
            [sequences[gene] for gene in cluster],
            [gene_names[gene] for gene in cluster],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build family profile {family_id}")
        profiles[family_id] = profile

    for pair_id in range(families // 2):
        cluster = (
            family_genes[pair_id * 2][:3]
            + family_genes[pair_id * 2 + 1][:3]
        )
        profile_id = len(profiles)
        clusters.append(cluster)
        profile = build_msa_profile(
            [sequences[gene] for gene in cluster],
            [gene_names[gene] for gene in cluster],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build mixed profile {profile_id}")
        profiles[profile_id] = profile

    selected_profiles = select_single_copy_profile_ids(
        profiles, clusters, gene_to_species
    )
    strict_thresholds = compute_profile_self_thresholds(
        profiles, gene_names, database, cpu
    )
    global_thresholds = compute_profile_self_thresholds(
        profiles,
        gene_names,
        database,
        cpu,
        matrix_name=MATRIX_NAME,
        calibrate_weakest_member=True,
    )
    selective_thresholds = compute_profile_self_thresholds(
        profiles,
        gene_names,
        database,
        cpu,
        matrix_name=MATRIX_NAME,
        calibrate_weakest_member=True,
        calibration_profile_ids=selected_profiles,
    )

    query_count = len(dataset.queries)
    pairs = np.asarray(
        [
            (profile_id, query_start + query_id)
            for profile_id in range(len(profiles))
            for query_id in range(query_count)
        ],
        dtype=np.int32,
    )
    _profile_ids, profile_lengths, scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    query_scores = scores.reshape(len(profiles), query_count)
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        scores,
        np.repeat(profile_lengths, query_count),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    ).reshape(query_scores.shape)

    evaluations = {}
    for name, thresholds in (
        ("strict", strict_thresholds),
        ("global_jackknife", global_thresholds),
        ("single_copy_jackknife", selective_thresholds),
    ):
        predictions = _predictions(query_scores, evalues, thresholds)
        metrics = classification_metrics(dataset.queries, predictions)
        metrics["contaminated_profile_wins"] = sum(
            prediction is not None and prediction >= families
            for prediction in predictions
        )
        evaluations[name] = metrics

    contaminated_ids = set(profiles) - set(range(families))
    contaminated_kept_strict = sum(
        selective_thresholds[profile_id] == strict_thresholds[profile_id]
        for profile_id in contaminated_ids
    )
    passed = (
        selected_profiles == set(range(families))
        and contaminated_kept_strict == len(contaminated_ids)
        and evaluations["single_copy_jackknife"]["precision"] == 1.0
        and evaluations["single_copy_jackknife"]["negative_rejection_rate"]
        == 1.0
        and evaluations["single_copy_jackknife"]["recall"]
        == evaluations["global_jackknife"]["recall"]
    )
    return {
        "seed": seed,
        "families": families,
        "single_copy_profiles": families,
        "contaminated_profiles": len(contaminated_ids),
        "selected_profiles": len(selected_profiles),
        "globally_relaxed_contaminated_profiles": sum(
            global_thresholds[profile_id] < strict_thresholds[profile_id]
            for profile_id in contaminated_ids
        ),
        "selectively_relaxed_contaminated_profiles": sum(
            selective_thresholds[profile_id] < strict_thresholds[profile_id]
            for profile_id in contaminated_ids
        ),
        "contaminated_profiles_kept_strict": contaminated_kept_strict,
        "evaluations": evaluations,
        "passed": passed,
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
            "Synthetic development/holdout comparison of strict, global "
            "jackknife, and single-copy-gated jackknife profile thresholds."
        ),
        "label_policy": (
            "Uses generated family, contamination, and species labels only; "
            "no final benchmark labels or competing-tool output is read."
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
        "sha256": sha256_file(args.output),
        "passed": passed,
        "development_recall": development["evaluations"]
        ["single_copy_jackknife"]["recall"],
        "holdout_recall": holdout["evaluations"]
        ["single_copy_jackknife"]["recall"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
