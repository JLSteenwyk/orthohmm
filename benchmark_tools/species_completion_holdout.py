#!/usr/bin/env python3
"""Validate missing-species HMM completion without benchmark labels."""

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
    build_strict_profile_edges,
    compute_profile_self_thresholds,
)


MATRIX_NAME = "BLOSUM62"
EVALUE_THRESHOLD = 1e-4


def _assignments(edges, query_start: int, training_family) -> dict[int, int]:
    assignments = {}
    for source, target in zip(edges.sources, edges.targets):
        source = int(source)
        target = int(target)
        if source >= query_start:
            query_gene, member_gene = source, target
        else:
            query_gene, member_gene = target, source
        assignments[query_gene] = training_family[member_gene]
    return assignments


def evaluate_split(seed: int, families: int, cpu: int) -> dict:
    dataset = generate_synthetic_dataset(seed, families=families)
    matrix = get_matrix(MATRIX_NAME)
    background = get_background_freqs(MATRIX_NAME)
    gene_names = []
    sequences = []
    gene_to_species = []
    profiles = {}
    clusters = []
    training_family = {}

    for family_id, members in dataset.training_sequences.items():
        member_ids = [
            f"family_{family_id}_member_{member_index}"
            for member_index in range(3)
        ]
        member_genes = []
        for species_id, (member_id, sequence) in enumerate(
            zip(member_ids, members[:3])
        ):
            member_genes.append(len(gene_names))
            training_family[len(gene_names)] = family_id
            gene_names.append(member_id)
            sequences.append(sequence)
            gene_to_species.append(species_id)
        clusters.append(member_genes)
        profile = build_msa_profile(
            list(members[:3]), member_ids, matrix, background
        )
        if profile is None:
            raise RuntimeError(f"could not build family profile {family_id}")
        profiles[family_id] = profile

    candidates = []
    for family_id, members in dataset.training_sequences.items():
        for species_id in (3, 4):
            candidates.append((
                family_id,
                "missing_species",
                species_id,
                members[species_id],
            ))
    for query in dataset.queries:
        if query.category == "close":
            candidates.append((
                query.family_id,
                "represented_species",
                0,
                query.sequence,
            ))
        elif query.category in ("close_paralog", "unrelated"):
            candidates.append((
                None,
                query.category,
                4,
                query.sequence,
            ))

    query_start = len(gene_names)
    candidate_genes = []
    for candidate_index, (_family, _category, species_id, sequence) in enumerate(
        candidates
    ):
        candidate_genes.append(len(gene_names))
        gene_names.append(f"candidate_{candidate_index}")
        sequences.append(sequence)
        gene_to_species.append(species_id)
        clusters.append([len(gene_names) - 1])
    database = species_sequences(gene_names, sequences, "species_completion")
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
        for gene in candidate_genes
    ], dtype=np.int32)
    _profile_ids, profile_lengths, scores = _score_pairs(
        profiles, database, pairs, cpu
    )
    lam, k_parameter = get_ka_params(MATRIX_NAME)
    evalues = batch_estimate_evalues(
        scores,
        np.repeat(profile_lengths, len(candidate_genes)),
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

    hit_queries = significant_pairs[:, 1]
    hit_targets = np.fromiter(
        (clusters[int(profile_id)][0] for profile_id in significant_pairs[:, 0]),
        dtype=np.int32,
        count=len(significant_pairs),
    )
    hit_scores = np.ones(len(significant_pairs), dtype=np.float64)
    unrestricted_edges = build_strict_profile_edges(
        gene_names,
        clusters,
        profile_hits,
        thresholds,
        hit_queries,
        hit_targets,
        hit_scores,
    )
    completion_edges = build_strict_profile_edges(
        gene_names,
        clusters,
        profile_hits,
        thresholds,
        hit_queries,
        hit_targets,
        hit_scores,
        candidate_missing_species_only=True,
        gene_to_species=np.asarray(gene_to_species, dtype=np.int32),
    )
    unrestricted = _assignments(
        unrestricted_edges, query_start, training_family
    )
    completion = _assignments(
        completion_edges, query_start, training_family
    )

    by_category = {}
    for category in sorted(set(candidate[1] for candidate in candidates)):
        rows = [
            index
            for index, candidate in enumerate(candidates)
            if candidate[1] == category
        ]
        unrestricted_assigned = sum(
            candidate_genes[index] in unrestricted for index in rows
        )
        completion_assigned = sum(
            candidate_genes[index] in completion for index in rows
        )
        completion_correct = sum(
            (
                completion.get(candidate_genes[index])
                == candidates[index][0]
            )
            if candidates[index][0] is not None
            else candidate_genes[index] not in completion
            for index in rows
        )
        by_category[category] = {
            "candidates": len(rows),
            "unrestricted_assigned": unrestricted_assigned,
            "completion_assigned": completion_assigned,
            "completion_correct_or_rejected": completion_correct,
        }

    missing = by_category["missing_species"]
    represented = by_category["represented_species"]
    passed = (
        missing["completion_assigned"] == missing["unrestricted_assigned"]
        and missing["completion_correct_or_rejected"] == missing[
            "completion_assigned"
        ]
        and represented["unrestricted_assigned"] > 0
        and represented["completion_assigned"] == 0
        and by_category["close_paralog"]["completion_assigned"] == 0
        and by_category["unrelated"]["completion_assigned"] == 0
    )
    return {
        "seed": seed,
        "families": families,
        "source_species_per_profile": 3,
        "candidate_categories": by_category,
        "unrestricted_edges": len(unrestricted_edges),
        "species_completion_edges": len(completion_edges),
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
            "Synthetic development/holdout test of a selective second HMM "
            "pass for missing species in incomplete single-copy groups."
        ),
        "label_policy": (
            "Uses generated family and species labels only; no OrthoBench, "
            "QfO, Three Kingdoms, or competing-tool output is read."
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
        "development_edges": development["species_completion_edges"],
        "holdout_edges": holdout["species_completion_edges"],
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
