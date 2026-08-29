#!/usr/bin/env python3
"""Validate profile-HMM species support on synthetic development/holdout data."""

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
    _mutate_sequence,
    generate_synthetic_dataset,
    species_sequences,
)
from orthohmm.search.profile_expansion import build_cluster_profiles


def evaluate_species_gate(seed: int, families: int, cpu: int) -> dict:
    """Build true and lineage-specific profiles with and without the gate."""
    dataset = generate_synthetic_dataset(seed, families=families)
    rng = np.random.default_rng(seed + 7919)
    ids: list[str] = []
    sequences: list[str] = []
    gene_to_species: list[int] = []
    clusters: list[list[int]] = []

    for family_id, members in dataset.training_sequences.items():
        cluster = []
        for member_index, sequence in enumerate(members):
            cluster.append(len(ids))
            ids.append(f"family_{family_id}_species_{member_index}")
            sequences.append(sequence)
            gene_to_species.append(member_index)
        clusters.append(cluster)

    for pair_index in range(families // 2):
        paralog = next(
            query.sequence
            for query in dataset.queries
            if query.query_id == f"pair_{pair_index}_close_paralog"
        )
        rates = np.ones(len(paralog), dtype=np.float64)
        paralogs = [
            paralog,
            _mutate_sequence(paralog, rates, 0.08, rng),
            _mutate_sequence(paralog, rates, 0.10, rng),
        ]
        cluster = []
        for copy_index, sequence in enumerate(paralogs):
            cluster.append(len(ids))
            ids.append(f"pair_{pair_index}_lineage_paralog_{copy_index}")
            sequences.append(sequence)
            gene_to_species.append(0)
        clusters.append(cluster)

    database = species_sequences(ids, sequences, "species_support")
    ungated = build_cluster_profiles(
        clusters,
        ids,
        database,
        "BLOSUM62",
        cpu=cpu,
        gene_to_species=gene_to_species,
        min_species_count=1,
    )
    gated = build_cluster_profiles(
        clusters,
        ids,
        database,
        "BLOSUM62",
        cpu=cpu,
        gene_to_species=gene_to_species,
        min_species_count=3,
    )
    true_ids = set(range(families))
    lineage_specific_ids = set(range(families, len(clusters)))
    gated_ids = set(gated)
    passed = (
        true_ids <= gated_ids
        and not (lineage_specific_ids & gated_ids)
        and len(ungated) == len(clusters)
    )
    return {
        "seed": seed,
        "families": families,
        "minimum_species": 3,
        "ungated_profiles": len(ungated),
        "gated_profiles": len(gated),
        "true_profiles": len(true_ids),
        "true_profiles_retained": len(true_ids & gated_ids),
        "lineage_specific_profiles": len(lineage_specific_ids),
        "lineage_specific_profiles_retained": len(
            lineage_specific_ids & gated_ids
        ),
        "passed": passed,
    }


def _git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        check=False,
        capture_output=True,
        text=True,
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

    development = evaluate_species_gate(
        args.development_seed, args.families, args.cpu
    )
    holdout = evaluate_species_gate(
        args.holdout_seed, args.families, args.cpu
    )
    payload = {
        "schema_version": 1,
        "description": (
            "Synthetic development/holdout construction test for requiring "
            "cross-species support before building a profile HMM."
        ),
        "label_policy": (
            "Uses generated family and species labels only; no final benchmark "
            "labels or tool outputs are read."
        ),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git": {
            "commit": _git_output("rev-parse", "HEAD"),
            "dirty": bool(_git_output("status", "--short")),
        },
        "configuration": {
            "development_seed": args.development_seed,
            "holdout_seed": args.holdout_seed,
            "families": args.families,
            "cpu": args.cpu,
            "matrix": "BLOSUM62",
            "minimum_species": 3,
        },
        "development": development,
        "holdout": holdout,
        "passed": development["passed"] and holdout["passed"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "passed": payload["passed"],
    }, indent=2))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
