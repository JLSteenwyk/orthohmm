#!/usr/bin/env python3
"""Calibrate two-member HMM profiles without final benchmark labels."""

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
from orthohmm.search.profile_expansion import compute_profile_self_thresholds


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4
DEVELOPMENT_RATIOS = (0.6, 0.7, 0.8)
LOCKED_RATIO = 0.7


def _predictions(
    scores: np.ndarray,
    evalues: np.ndarray,
    thresholds: dict[int, float],
) -> list[int | None]:
    predictions = []
    for query_index in range(scores.shape[1]):
        eligible = [
            profile_id
            for profile_id in range(scores.shape[0])
            if evalues[profile_id, query_index] < EVALUE_THRESHOLD
            and scores[profile_id, query_index] >= thresholds[profile_id]
        ]
        predictions.append(
            max(
                eligible,
                key=lambda profile_id: (
                    scores[profile_id, query_index], -profile_id
                ),
            )
            if eligible
            else None
        )
    return predictions


def evaluate_split(
    seed: int,
    families: int,
    cpu: int,
    ratios: Sequence[float],
) -> dict:
    dataset = generate_synthetic_dataset(seed, families=families)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names = []
    sequences = []
    profiles = {}

    for family_id, members in dataset.training_sequences.items():
        pair = members[:2]
        member_ids = [
            f"family_{family_id}_member_{member_index}"
            for member_index in range(2)
        ]
        profile = build_msa_profile(
            list(pair), member_ids, matrix, background
        )
        if profile is None:
            raise RuntimeError(f"could not build pair profile {family_id}")
        profiles[family_id] = profile
        gene_names.extend(member_ids)
        sequences.extend(pair)

    query_start = len(gene_names)
    gene_names.extend(query.query_id for query in dataset.queries)
    sequences.extend(query.sequence for query in dataset.queries)
    database = species_sequences(gene_names, sequences, "pair_profiles")
    query_count = len(dataset.queries)
    pairs = np.asarray([
        (profile_id, query_start + query_index)
        for profile_id in range(families)
        for query_index in range(query_count)
    ], dtype=np.int32)
    _profile_ids, profile_lengths, flat_scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    scores = flat_scores.reshape(families, query_count)
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        flat_scores,
        np.repeat(profile_lengths, query_count),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    ).reshape(scores.shape)

    evaluations = {}
    for ratio in ratios:
        thresholds = compute_profile_self_thresholds(
            profiles,
            gene_names,
            database,
            cpu,
            matrix_name=MATRIX_NAME,
            calibrate_weakest_member=True,
            calibration_profile_ids=set(profiles),
            pair_profile_threshold_ratio=ratio,
        )
        metrics = classification_metrics(
            dataset.queries,
            _predictions(scores, evalues, thresholds),
        )
        evaluations[str(ratio)] = {
            **metrics,
            "threshold_minimum": min(thresholds.values()),
            "threshold_median": float(np.median(list(thresholds.values()))),
            "threshold_maximum": max(thresholds.values()),
        }
    return {
        "seed": seed,
        "families": families,
        "training_members_per_profile": 2,
        "queries": query_count,
        "evaluations": evaluations,
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

    development = evaluate_split(
        args.development_seed,
        args.families,
        args.cpu,
        DEVELOPMENT_RATIOS,
    )
    valid_ratios = [
        ratio
        for ratio in DEVELOPMENT_RATIOS
        if development["evaluations"][str(ratio)]["precision"] == 1.0
        and development["evaluations"][str(ratio)][
            "negative_rejection_rate"
        ] == 1.0
    ]
    selected_ratio = min(valid_ratios) if valid_ratios else None
    holdout = evaluate_split(
        args.holdout_seed,
        args.families,
        args.cpu,
        (LOCKED_RATIO,),
    )
    holdout_metrics = holdout["evaluations"][str(LOCKED_RATIO)]
    passed = (
        selected_ratio == LOCKED_RATIO
        and holdout_metrics["precision"] == 1.0
        and holdout_metrics["negative_rejection_rate"] == 1.0
        and holdout_metrics["recall"] >= 0.8
    )
    payload = {
        "schema_version": 1,
        "description": (
            "Synthetic development/holdout calibration of two-member "
            "profile-HMM thresholds."
        ),
        "label_policy": (
            "Uses generated family labels only; no OrthoBench, QfO, Three "
            "Kingdoms, or competing-tool output is read."
        ),
        "selection_policy": (
            "Choose the smallest development ratio in the prespecified grid "
            "with precision and negative rejection both equal to 1.0, then "
            "apply it unchanged to holdout."
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
            "development_ratios": list(DEVELOPMENT_RATIOS),
            "locked_ratio": LOCKED_RATIO,
        },
        "selected_development_ratio": selected_ratio,
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
        "selected_development_ratio": selected_ratio,
        "holdout_f_score": holdout_metrics["f_score"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
