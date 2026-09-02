#!/usr/bin/env python3
"""Audit candidate-superfamily recall ceilings and seed-family growth.

This post-inference diagnostic may read benchmark labels, but candidate
construction never does.  It verifies that candidate families are a strict
coarsening of immutable seed groups and quantifies the recall/contamination
tradeoff before gene-tree reconciliation is run.
"""

from __future__ import annotations

import argparse
import itertools
import json
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
    read_clusters,
    run_official_benchmark,
)


def _partition_index(clusters: Sequence[set[str]], label: str) -> dict[str, int]:
    index: dict[str, int] = {}
    for cluster_index, cluster in enumerate(clusters):
        for gene in cluster:
            if gene in index:
                raise ValueError(
                    f"{label} contains gene {gene!r} in more than one family"
                )
            index[gene] = cluster_index
    return index


def _quantiles(values: Sequence[int]) -> dict[str, float | int]:
    if not values:
        return {"minimum": 0, "median": 0.0, "p90": 0.0, "maximum": 0}
    ordered = sorted(values)

    def percentile(fraction: float) -> float:
        position = fraction * (len(ordered) - 1)
        lower = int(position)
        upper = min(lower + 1, len(ordered) - 1)
        weight = position - lower
        return ordered[lower] * (1.0 - weight) + ordered[upper] * weight

    return {
        "minimum": ordered[0],
        "median": statistics.median(ordered),
        "p90": percentile(0.9),
        "maximum": ordered[-1],
    }


def audit_candidate_partition(
    seeds: Sequence[set[str]],
    candidates: Sequence[set[str]],
    refogs: Sequence[set[str]],
) -> dict:
    """Return label-aware diagnostics for a frozen, label-blind partition."""

    seed_index = _partition_index(seeds, "seed partition")
    candidate_index = _partition_index(candidates, "candidate partition")
    if seed_index.keys() != candidate_index.keys():
        missing = sorted(seed_index.keys() - candidate_index.keys())
        extra = sorted(candidate_index.keys() - seed_index.keys())
        raise ValueError(
            "candidate and seed gene universes differ "
            f"(missing={missing[:5]}, extra={extra[:5]})"
        )

    candidate_to_seeds: list[set[int]] = [set() for _ in candidates]
    for seed_id, seed in enumerate(seeds):
        destinations = {candidate_index[gene] for gene in seed}
        if len(destinations) != 1:
            raise ValueError(
                f"candidate partition splits immutable seed family {seed_id}"
            )
        candidate_to_seeds[destinations.pop()].add(seed_id)

    gene_refogs: dict[str, set[int]] = {}
    for refog_id, group in enumerate(refogs):
        for gene in group:
            gene_refogs.setdefault(gene, set()).add(refog_id)
    reference_genes = set(gene_refogs)

    true_pairs: set[tuple[int, str, str]] = set()
    cross_refog_pairs: set[tuple[str, str]] = set()
    reference_nonreference_relations = 0
    mixed_reference_candidates = 0
    for candidate in candidates:
        retained = sorted(candidate & reference_genes)
        reference_nonreference_relations += len(retained) * (
            len(candidate) - len(retained)
        )
        represented_refogs = set().union(
            *(gene_refogs[gene] for gene in retained)
        ) if retained else set()
        if len(represented_refogs) > 1:
            mixed_reference_candidates += 1
        for left, right in itertools.combinations(retained, 2):
            common = gene_refogs[left] & gene_refogs[right]
            if common:
                true_pairs.update((refog, left, right) for refog in common)
            else:
                cross_refog_pairs.add((left, right))

    possible_pairs = sum(
        len(group) * (len(group) - 1) // 2 for group in refogs
    )
    pair_denominator = len(true_pairs) + len(cross_refog_pairs)
    refog_records = []
    macro_recalls = []
    fully_colocated = 0
    for refog_id, group in enumerate(refogs):
        possible = len(group) * (len(group) - 1) // 2
        recovered = sum(
            candidate_index.get(left) == candidate_index.get(right)
            for left, right in itertools.combinations(sorted(group), 2)
            if left in candidate_index and right in candidate_index
        )
        recall = recovered / possible if possible else 1.0
        macro_recalls.append(recall)
        observed = group & candidate_index.keys()
        candidate_components = {
            candidate_index[gene] for gene in observed
        }
        seed_components = {seed_index[gene] for gene in observed}
        if observed == group and len(candidate_components) == 1:
            fully_colocated += 1
        refog_records.append({
            "refog": f"RefOG{refog_id + 1:03d}",
            "genes": len(group),
            "genes_observed": len(observed),
            "seed_components": len(seed_components),
            "candidate_components": len(candidate_components),
            "possible_pairs": possible,
            "recovered_pairs": recovered,
            "pair_recall": recall,
        })

    seed_counts = [len(group) for group in candidate_to_seeds]
    family_sizes = [len(group) for group in candidates]
    merged = [count for count in seed_counts if count > 1]
    return {
        "genes": len(seed_index),
        "seed_families": len(seeds),
        "candidate_families": len(candidates),
        "seed_merges": len(seeds) - len(candidates),
        "candidates_combining_seeds": len(merged),
        "seed_families_per_candidate": _quantiles(seed_counts),
        "candidate_family_sizes": _quantiles(family_sizes),
        "largest_candidate_fraction": (
            max(family_sizes, default=0) / len(seed_index) if seed_index else 0.0
        ),
        "reference": {
            "refogs": len(refogs),
            "reference_genes": len(reference_genes),
            "possible_pairs": possible_pairs,
            "recovered_pairs": len(true_pairs),
            "pair_recall_micro": (
                len(true_pairs) / possible_pairs if possible_pairs else 0.0
            ),
            "pair_recall_macro": statistics.fmean(macro_recalls)
            if macro_recalls else 0.0,
            "pair_precision_among_reference_genes": (
                len(true_pairs) / pair_denominator if pair_denominator else 0.0
            ),
            "cross_refog_pairs": len(cross_refog_pairs),
            "mixed_reference_candidates": mixed_reference_candidates,
            "reference_to_nonreference_relations": (
                reference_nonreference_relations
            ),
            "fully_colocated_refogs": fully_colocated,
        },
        "refog_records": refog_records,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-clusters", required=True, type=Path)
    parser.add_argument("--candidate-clusters", required=True, type=Path)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--official-benchmark", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    seeds = read_clusters(args.seed_clusters)
    candidates = read_clusters(args.candidate_clusters)
    refogs = load_refogs(args.refogs)
    audit = audit_candidate_partition(seeds, candidates, refogs)
    payload = {
        "schema_version": 1,
        "description": "Candidate-superfamily pre-reconciliation audit",
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "inputs": {
            "seed_clusters": file_provenance(args.seed_clusters),
            "candidate_clusters": file_provenance(args.candidate_clusters),
            "refog_directory": str(args.refogs.resolve()),
        },
        "audit": audit,
    }
    if args.official_benchmark:
        payload["official_orthobench"] = run_official_benchmark(
            args.official_benchmark.resolve(), args.candidate_clusters.resolve()
        )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)
    summary = audit["reference"]
    print(
        f"candidate ceiling: recall={summary['pair_recall_micro']:.6f} "
        f"reference-pair precision="
        f"{summary['pair_precision_among_reference_genes']:.6f} "
        f"seed merges={audit['seed_merges']}"
    )
    print(args.json.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
