#!/usr/bin/env python3
"""Benchmark HMM profile decisions without using final benchmark labels.

The synthetic families in this module are deterministic but deliberately
heterogeneous: model lengths vary, pairs of families share a superfamily
ancestor, sites evolve at different rates, and held-out queries include
divergent orthologs, fragments, close paralogs, and unrelated proteins.

This is a development/holdout diagnostic, not a biological benchmark.  Its
purpose is to reject profile-search changes that only improve a named final
benchmark and to expose in-sample threshold optimism in a reproducible way.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from orthohmm.search.evalue import batch_estimate_evalues
from orthohmm.search.matrices import (
    ALPHABET,
    get_background_freqs,
    get_ka_params,
    get_matrix,
)
from orthohmm.search.msa_profile import MSAProfile, build_msa_profile
from orthohmm.search.prefilter import prefilter_candidates
from orthohmm.search.profile_expansion import _pack_profiles, _score_profile_pairs
from orthohmm.search.sequences import SpeciesSequences, encode_sequence


MATRIX_NAME = "BLOSUM62"
TRAIN_DISTANCES = (0.08, 0.12, 0.18, 0.24, 0.31)
POSITIVE_CATEGORIES = ("close", "standard", "distant", "fragment")
NEGATIVE_CATEGORIES = ("close_paralog", "unrelated")


@dataclass(frozen=True)
class SyntheticQuery:
    query_id: str
    sequence: str
    family_id: int | None
    category: str


@dataclass(frozen=True)
class SyntheticDataset:
    seed: int
    training_sequences: Mapping[int, tuple[str, ...]]
    queries: tuple[SyntheticQuery, ...]


def species_sequences(
    ids: Sequence[str],
    sequences: Sequence[str],
    name: str,
) -> SpeciesSequences:
    """Pack sequence strings in the same representation used by OrthoHMM."""
    encoded = [encode_sequence(sequence) for sequence in sequences]
    lengths = np.asarray([len(sequence) for sequence in encoded], dtype=np.int32)
    offsets = np.zeros(len(encoded), dtype=np.int64)
    if len(encoded) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.int64)
    flat = (
        np.concatenate(encoded).astype(np.uint8, copy=False)
        if encoded
        else np.empty(0, dtype=np.uint8)
    )
    return SpeciesSequences(name, list(ids), flat, offsets, lengths)


def _random_sequence(length: int, rng: np.random.Generator) -> str:
    background = get_background_freqs(MATRIX_NAME)
    return "".join(rng.choice(list(ALPHABET), size=length, p=background))


def _mutate_sequence(
    sequence: str,
    site_rates: np.ndarray,
    distance: float,
    rng: np.random.Generator,
) -> str:
    """Evolve a sequence with BLOSUM-weighted replacements."""
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    mutated = list(sequence)
    for position, amino_acid in enumerate(sequence):
        probability = min(0.95, distance * float(site_rates[position]))
        if rng.random() >= probability:
            continue
        source = ALPHABET.index(amino_acid)
        weights = np.exp2(matrix[source].astype(np.float64) / 2.0) * background
        weights[source] = 0.0
        weights /= weights.sum()
        mutated[position] = ALPHABET[int(rng.choice(20, p=weights))]
    return "".join(mutated)


def _shift_conserved_positions(
    sequence: str,
    site_rates: np.ndarray,
    rng: np.random.Generator,
) -> str:
    """Introduce family-specific residues at otherwise conserved sites."""
    conserved = np.flatnonzero(site_rates < 0.3)
    count = min(len(conserved), max(3, len(conserved) // 10))
    if count == 0:
        return sequence
    positions = rng.choice(conserved, size=count, replace=False)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    shifted = list(sequence)
    for position in positions:
        source = ALPHABET.index(shifted[int(position)])
        weights = np.exp2(matrix[source].astype(np.float64) / 2.0) * background
        weights[source] = 0.0
        weights /= weights.sum()
        shifted[int(position)] = ALPHABET[int(rng.choice(20, p=weights))]
    return "".join(shifted)


def generate_synthetic_dataset(seed: int, families: int = 36) -> SyntheticDataset:
    """Generate paired protein families and labeled held-out queries."""
    if families < 2 or families % 2:
        raise ValueError("families must be an even integer >= 2")
    rng = np.random.default_rng(seed)
    training: dict[int, tuple[str, ...]] = {}
    queries: list[SyntheticQuery] = []

    for pair_index in range(families // 2):
        length = int(rng.integers(80, 451))
        superfamily_root = _random_sequence(length, rng)
        site_rates = rng.choice(
            np.asarray([0.15, 0.8, 1.6]),
            size=length,
            p=np.asarray([0.25, 0.55, 0.20]),
        )

        for side in range(2):
            family_id = pair_index * 2 + side
            ancestor = _mutate_sequence(
                superfamily_root,
                site_rates,
                0.22 + 0.02 * side,
                rng,
            )
            ancestor = _shift_conserved_positions(ancestor, site_rates, rng)
            members = tuple(
                _mutate_sequence(ancestor, site_rates, distance, rng)
                for distance in TRAIN_DISTANCES
            )
            training[family_id] = members

            positive_sequences = {
                "close": _mutate_sequence(ancestor, site_rates, 0.16, rng),
                "standard": _mutate_sequence(ancestor, site_rates, 0.25, rng),
                "distant": _mutate_sequence(ancestor, site_rates, 0.38, rng),
            }
            fragment_source = _mutate_sequence(ancestor, site_rates, 0.23, rng)
            positive_sequences["fragment"] = fragment_source[
                int(0.10 * length):int(0.90 * length)
            ]
            for category, sequence in positive_sequences.items():
                queries.append(SyntheticQuery(
                    query_id=f"family_{family_id}_{category}",
                    sequence=sequence,
                    family_id=family_id,
                    category=category,
                ))

        close_paralog = _mutate_sequence(
            superfamily_root, site_rates, 0.34, rng
        )
        close_paralog = _shift_conserved_positions(
            close_paralog, site_rates, rng
        )
        queries.append(SyntheticQuery(
            query_id=f"pair_{pair_index}_close_paralog",
            sequence=close_paralog,
            family_id=None,
            category="close_paralog",
        ))
        queries.append(SyntheticQuery(
            query_id=f"pair_{pair_index}_unrelated",
            sequence=_random_sequence(length, rng),
            family_id=None,
            category="unrelated",
        ))

    return SyntheticDataset(seed, training, tuple(queries))


def _build_profiles(dataset: SyntheticDataset) -> dict[int, MSAProfile]:
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    profiles = {}
    for family_id, sequences in dataset.training_sequences.items():
        profile = build_msa_profile(
            list(sequences),
            [f"family_{family_id}_member_{index}" for index in range(len(sequences))],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build synthetic profile {family_id}")
        profiles[family_id] = profile
    return profiles


def _score_pairs(
    profiles: Mapping[int, MSAProfile],
    database: SpeciesSequences,
    pairs: np.ndarray,
    cpu: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    profile_ids, packed = _pack_profiles(profiles)
    if packed is None:
        raise ValueError("at least one profile is required")
    profile_lengths, profile_offsets, emissions, _consensus, inserts, transitions = packed
    scores = _score_profile_pairs(
        emissions,
        inserts,
        transitions,
        profile_offsets,
        profile_lengths,
        database,
        np.asarray(pairs, dtype=np.int32),
        band_width=64,
        cpu=cpu,
    )
    return profile_ids, profile_lengths, scores.astype(np.float64)


def _calibration_thresholds(
    dataset: SyntheticDataset,
    full_scores: np.ndarray,
    member_families: np.ndarray,
    cpu: int,
) -> tuple[dict[str, dict[int, float]], dict]:
    """Calculate in-sample, one-jackknife, and full LOO thresholds."""
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    in_sample: dict[int, float] = {}
    weakest_indices: dict[int, int] = {}
    for family_id in sorted(dataset.training_sequences):
        member_columns = np.flatnonzero(member_families == family_id)
        family_scores = full_scores[family_id, member_columns]
        weakest_local = int(np.argmin(family_scores))
        weakest_indices[family_id] = weakest_local
        in_sample[family_id] = float(family_scores[weakest_local])

    jackknife_profiles: dict[int, MSAProfile] = {}
    jackknife_sequences: list[str] = []
    for family_id, sequences in dataset.training_sequences.items():
        held_index = weakest_indices[family_id]
        retained = list(sequences[:held_index] + sequences[held_index + 1:])
        profile = build_msa_profile(
            retained,
            [f"retained_{index}" for index in range(len(retained))],
            matrix,
            background,
        )
        if profile is None:
            raise RuntimeError(f"could not build jackknife profile {family_id}")
        jackknife_profiles[family_id] = profile
        jackknife_sequences.append(sequences[held_index])
    jackknife_database = species_sequences(
        [str(family_id) for family_id in sorted(jackknife_profiles)],
        jackknife_sequences,
        "weakest_jackknife",
    )
    family_count = len(jackknife_profiles)
    _, _, jackknife_scores = _score_pairs(
        jackknife_profiles,
        jackknife_database,
        np.column_stack((np.arange(family_count), np.arange(family_count))),
        cpu,
    )
    weakest_jackknife = {
        family_id: min(in_sample[family_id], float(jackknife_scores[index]))
        for index, family_id in enumerate(sorted(jackknife_profiles))
    }

    loo_profiles: dict[int, MSAProfile] = {}
    loo_sequences: list[str] = []
    loo_families: list[int] = []
    profile_index = 0
    for family_id, sequences in dataset.training_sequences.items():
        for held_index, held_sequence in enumerate(sequences):
            retained = list(sequences[:held_index] + sequences[held_index + 1:])
            profile = build_msa_profile(
                retained,
                [f"retained_{index}" for index in range(len(retained))],
                matrix,
                background,
            )
            if profile is None:
                raise RuntimeError(
                    f"could not build leave-one-out profile {family_id}:{held_index}"
                )
            loo_profiles[profile_index] = profile
            loo_sequences.append(held_sequence)
            loo_families.append(family_id)
            profile_index += 1
    loo_database = species_sequences(
        [str(index) for index in range(len(loo_sequences))],
        loo_sequences,
        "full_leave_one_out",
    )
    _, _, loo_scores = _score_pairs(
        loo_profiles,
        loo_database,
        np.column_stack((np.arange(profile_index), np.arange(profile_index))),
        cpu,
    )
    full_leave_one_out = {
        family_id: min(
            in_sample[family_id],
            min(
                float(score)
                for score, score_family in zip(loo_scores, loo_families)
                if score_family == family_id
            ),
        )
        for family_id in sorted(dataset.training_sequences)
    }

    jackknife_ratios = np.asarray([
        weakest_jackknife[family_id] / in_sample[family_id]
        for family_id in sorted(in_sample)
    ])
    loo_ratios = np.asarray([
        full_leave_one_out[family_id] / in_sample[family_id]
        for family_id in sorted(in_sample)
    ])
    diagnostics = {
        "profiles": family_count,
        "jackknife_profiles_built": family_count,
        "leave_one_out_profiles_built": profile_index,
        "weakest_jackknife_to_in_sample_ratio": _distribution(jackknife_ratios),
        "full_leave_one_out_to_in_sample_ratio": _distribution(loo_ratios),
    }
    return {
        "in_sample": in_sample,
        "weakest_jackknife": weakest_jackknife,
        "full_leave_one_out": full_leave_one_out,
    }, diagnostics


def _distribution(values: np.ndarray) -> dict[str, float]:
    return {
        "minimum": float(np.min(values)),
        "median": float(np.median(values)),
        "maximum": float(np.max(values)),
    }


def select_profile(
    eligible: Sequence[tuple[int, float, float, float]],
    method: str,
) -> int | None:
    """Select one profile from (id, score, threshold, model_length) rows."""
    if not eligible:
        return None
    if method == "raw_score":
        value = lambda row: row[1]
    elif method == "margin_per_model_position":
        value = lambda row: (row[1] - row[2]) / max(row[3], 1.0)
    else:
        raise ValueError(f"unknown profile-selection method: {method}")
    return int(max(eligible, key=lambda row: (value(row), -row[0]))[0])


def classification_metrics(
    queries: Sequence[SyntheticQuery],
    predictions: Sequence[int | None],
) -> dict:
    """Return family-assignment and negative-rejection metrics."""
    if len(queries) != len(predictions):
        raise ValueError("queries and predictions must have equal lengths")
    true_positives = 0
    positive_count = 0
    negative_count = 0
    predictions_made = 0
    negative_rejections = 0
    correct_or_rejected = 0
    by_category: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    for query, prediction in zip(queries, predictions):
        assigned = prediction is not None
        correct = prediction == query.family_id if query.family_id is not None else not assigned
        by_category[query.category][2] += 1
        by_category[query.category][1] += int(assigned)
        by_category[query.category][0] += int(correct)
        predictions_made += int(assigned)
        correct_or_rejected += int(correct)
        if query.family_id is None:
            negative_count += 1
            negative_rejections += int(not assigned)
        else:
            positive_count += 1
            true_positives += int(prediction == query.family_id)

    precision = true_positives / predictions_made if predictions_made else 0.0
    recall = true_positives / positive_count if positive_count else 0.0
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "true_family_assignments": true_positives,
        "positive_queries": positive_count,
        "negative_queries": negative_count,
        "predictions_made": predictions_made,
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "negative_rejection_rate": (
            negative_rejections / negative_count if negative_count else 0.0
        ),
        "overall_accuracy": correct_or_rejected / len(queries) if queries else 0.0,
        "by_category": {
            category: {
                "correct": values[0],
                "assigned": values[1],
                "queries": values[2],
                "accuracy": values[0] / values[2],
            }
            for category, values in sorted(by_category.items())
        },
    }


def _prefilter_candidate_sets(
    profiles: Mapping[int, MSAProfile],
    queries: Sequence[SyntheticQuery],
    cpu: int,
) -> tuple[dict[int, list[set[int]]], dict]:
    profile_ids, packed = _pack_profiles(profiles)
    if packed is None:
        raise ValueError("at least one profile is required")
    lengths, offsets, _emissions, consensus, _inserts, _transitions = packed
    profile_database = SpeciesSequences(
        "synthetic_profiles",
        [str(profile_id) for profile_id in profile_ids],
        consensus,
        offsets,
        lengths,
    )
    query_database = species_sequences(
        [query.query_id for query in queries],
        [query.sequence for query in queries],
        "synthetic_queries",
    )
    matrix = get_matrix(MATRIX_NAME)
    candidate_sets = {}
    diagnostics = {}
    for kmer_k in (4, 3):
        candidates, candidate_offsets = prefilter_candidates(
            query_database,
            profile_database,
            k=kmer_k,
            use_reduced_alphabet=False,
            min_total_hits=4,
            min_diag_hits=1,
            max_candidates_per_query=20,
            n_threads=cpu,
            sub_matrix=matrix,
        )
        per_query = []
        true_found = 0
        positive_queries = 0
        for query_index, query in enumerate(queries):
            local = candidates[
                candidate_offsets[query_index]:candidate_offsets[query_index + 1]
            ]
            selected = set(int(item) for item in profile_ids[local])
            per_query.append(selected)
            if query.family_id is not None:
                positive_queries += 1
                true_found += int(query.family_id in selected)
        candidate_sets[kmer_k] = per_query
        diagnostics[str(kmer_k)] = {
            "candidate_pairs": int(len(candidates)),
            "mean_candidates_per_query": len(candidates) / len(queries),
            "true_profile_recall": true_found / positive_queries,
            "max_candidates_per_query": 20,
        }
    return candidate_sets, diagnostics


def evaluate_split(
    seed: int,
    families: int = 36,
    cpu: int = 4,
) -> dict:
    """Build and evaluate one synthetic family split."""
    dataset = generate_synthetic_dataset(seed, families=families)
    profiles = _build_profiles(dataset)

    member_ids: list[str] = []
    member_sequences: list[str] = []
    member_families: list[int] = []
    for family_id, sequences in dataset.training_sequences.items():
        for member_index, sequence in enumerate(sequences):
            member_ids.append(f"family_{family_id}_member_{member_index}")
            member_sequences.append(sequence)
            member_families.append(family_id)
    all_ids = member_ids + [query.query_id for query in dataset.queries]
    all_sequences = member_sequences + [query.sequence for query in dataset.queries]
    database = species_sequences(all_ids, all_sequences, "synthetic_all")
    profile_count = len(profiles)
    pairs = np.asarray([
        (profile_index, sequence_index)
        for profile_index in range(profile_count)
        for sequence_index in range(len(all_sequences))
    ], dtype=np.int32)
    profile_ids, profile_lengths, flat_scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    score_matrix = flat_scores.reshape(profile_count, len(all_sequences))
    thresholds, calibration_diagnostics = _calibration_thresholds(
        dataset,
        score_matrix,
        np.asarray(member_families, dtype=np.int32),
        cpu,
    )

    query_scores = score_matrix[:, len(member_sequences):]
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        query_scores.reshape(-1),
        np.repeat(profile_lengths, len(dataset.queries)),
        int(database.lengths.sum()),
        lam,
        k_parameter,
    ).reshape(query_scores.shape)
    significant = evalues < 1e-4
    candidate_sets, prefilter_diagnostics = _prefilter_candidate_sets(
        profiles, dataset.queries, cpu
    )

    evaluations = {}
    for threshold_name, family_thresholds in thresholds.items():
        for selection_method in ("raw_score", "margin_per_model_position"):
            for candidate_scope in ("all_profiles", "k4", "k3_fallback"):
                predictions = []
                ambiguous = 0
                fallback_queries = 0
                for query_index, _query in enumerate(dataset.queries):
                    if candidate_scope == "all_profiles":
                        allowed = set(int(item) for item in profile_ids)
                    else:
                        allowed = candidate_sets[4][query_index]
                    eligible = []
                    for local_profile, family_id in enumerate(profile_ids):
                        family_id = int(family_id)
                        if family_id not in allowed or not significant[local_profile, query_index]:
                            continue
                        score = float(query_scores[local_profile, query_index])
                        threshold = family_thresholds[family_id]
                        if score >= threshold:
                            eligible.append((
                                family_id,
                                score,
                                threshold,
                                float(profile_lengths[local_profile]),
                            ))
                    if candidate_scope == "k3_fallback" and not eligible:
                        fallback_queries += 1
                        allowed = candidate_sets[3][query_index]
                        for local_profile, family_id in enumerate(profile_ids):
                            family_id = int(family_id)
                            if family_id not in allowed or not significant[local_profile, query_index]:
                                continue
                            score = float(query_scores[local_profile, query_index])
                            threshold = family_thresholds[family_id]
                            if score >= threshold:
                                eligible.append((
                                    family_id,
                                    score,
                                    threshold,
                                    float(profile_lengths[local_profile]),
                                ))
                    ambiguous += int(len(eligible) > 1)
                    predictions.append(select_profile(eligible, selection_method))
                key = f"{threshold_name}__{selection_method}__{candidate_scope}"
                evaluations[key] = {
                    **classification_metrics(dataset.queries, predictions),
                    "ambiguous_queries": ambiguous,
                    "fallback_queries": fallback_queries,
                }

    return {
        "seed": seed,
        "families": families,
        "training_members_per_family": len(TRAIN_DISTANCES),
        "queries": len(dataset.queries),
        "positive_queries": sum(
            query.family_id is not None for query in dataset.queries
        ),
        "negative_queries": sum(
            query.family_id is None for query in dataset.queries
        ),
        "profile_length": _distribution(profile_lengths.astype(np.float64)),
        "calibration": calibration_diagnostics,
        "prefilter": prefilter_diagnostics,
        "evaluations": evaluations,
    }


def _git_output(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], check=False, capture_output=True, text=True
    )
    return completed.stdout.strip()


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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = {
        "schema_version": 1,
        "description": (
            "Benchmark-independent synthetic protein-family development and "
            "holdout evaluation for OrthoHMM profile calibration."
        ),
        "label_policy": (
            "Uses generated family labels only; no OrthoBench, QfO, "
            "Three Kingdoms, or OrthoFinder outputs are read."
        ),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "configuration": {
            "matrix": MATRIX_NAME,
            "development_seed": args.development_seed,
            "holdout_seed": args.holdout_seed,
            "families": args.families,
            "cpu": args.cpu,
        },
        "git": {
            "commit": _git_output("rev-parse", "HEAD"),
            "status_short": _git_output("status", "--short"),
        },
        "development": evaluate_split(
            args.development_seed, args.families, args.cpu
        ),
        "holdout": evaluate_split(
            args.holdout_seed, args.families, args.cpu
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "development_evaluations": len(payload["development"]["evaluations"]),
        "holdout_evaluations": len(payload["holdout"]["evaluations"]),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
