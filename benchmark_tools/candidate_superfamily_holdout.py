#!/usr/bin/env python3
"""Independent synthetic holdout for reciprocal candidate-family merging."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import random
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
)
from orthohmm.refinement import (
    merge_reciprocal_candidate_clusters,
    merge_supported_satellite_candidate_clusters,
)


def build_case(seed: int, families: int = 24):
    """Create fragmented families plus weaker asymmetric decoy evidence."""
    rng = random.Random(seed)
    clusters = []
    truth = []
    species = []
    family_cluster_ids = []
    for family_id in range(families):
        fragments = rng.randint(2, 5)
        cluster_ids = []
        family_genes = []
        for fragment in range(fragments):
            gene = len(species)
            species.append((family_id + fragment) % 6)
            cluster_ids.append(len(clusters))
            clusters.append([gene])
            family_genes.append(gene)
        family_cluster_ids.append(cluster_ids)
        truth.append(family_genes)

    queries = []
    targets = []
    scores = []
    for cluster_ids in family_cluster_ids:
        for offset, (left, right) in enumerate(zip(cluster_ids, cluster_ids[1:])):
            score = 8.0 - offset * 0.5
            left_gene = clusters[left][0]
            right_gene = clusters[right][0]
            queries.extend([left_gene, right_gene])
            targets.extend([right_gene, left_gene])
            scores.extend([score, score])

    # A directed cycle supplies realistic weak decoys without reciprocal pairs.
    for family_id, cluster_ids in enumerate(family_cluster_ids):
        next_family = family_cluster_ids[(family_id + 1) % families]
        queries.append(clusters[cluster_ids[-1]][0])
        targets.append(clusters[next_family[0]][0])
        scores.append(0.1 + rng.random() * 0.01)
    return (
        clusters,
        np.asarray(queries, dtype=np.int32),
        np.asarray(targets, dtype=np.int32),
        np.asarray(scores, dtype=np.float64),
        np.asarray(species, dtype=np.int32),
        truth,
    )


def pair_score(predicted, expected) -> dict:
    def pairs(groups):
        return {
            tuple(sorted(pair))
            for group in groups
            for pair in combinations(group, 2)
        }

    predicted_pairs = pairs(predicted)
    expected_pairs = pairs(expected)
    true_positive = len(predicted_pairs & expected_pairs)
    precision = true_positive / len(predicted_pairs) if predicted_pairs else 0.0
    recall = true_positive / len(expected_pairs) if expected_pairs else 0.0
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "predicted_pairs": len(predicted_pairs),
        "expected_pairs": len(expected_pairs),
        "true_positive_pairs": true_positive,
    }


def evaluate(seed: int, iterations: int) -> dict:
    clusters, queries, targets, scores, species, truth = build_case(seed)
    predicted, merges, relations, completed = (
        merge_reciprocal_candidate_clusters(
            clusters,
            queries,
            targets,
            scores,
            species,
            max_component_genes=5,
            min_avg_score=0.05,
            min_max_score=0.05,
            min_coverage=0.0,
            min_norm=0.0,
            max_iterations=iterations,
        )
    )
    return {
        "seed": seed,
        "requested_iterations": iterations,
        "completed_iterations": completed,
        "seed_fragments": len(clusters),
        "truth_families": len(truth),
        "candidate_families": len(predicted),
        "merges": merges,
        "directed_cluster_relations": relations,
        "pair_score": pair_score(predicted, truth),
    }


def evaluate_satellite(seed: int, iterations: int = 1) -> dict:
    clusters, queries, targets, scores, species, truth = build_case(seed)
    predicted, merges, relations, completed = (
        merge_supported_satellite_candidate_clusters(
            clusters,
            queries,
            targets,
            scores,
            species,
            max_component_genes=5,
            max_satellite_genes=1,
            max_satellite_to_anchor_ratio=1.0,
            min_margin=1.0,
            max_satellites_per_anchor=4,
            min_avg_score=0.05,
            min_max_score=0.05,
            min_coverage=0.0,
            min_norm=0.0,
            max_iterations=iterations,
        )
    )
    return {
        "seed": seed,
        "requested_iterations": iterations,
        "completed_iterations": completed,
        "seed_fragments": len(clusters),
        "truth_families": len(truth),
        "candidate_families": len(predicted),
        "merges": merges,
        "directed_cluster_relations": relations,
        "pair_score": pair_score(predicted, truth),
    }


def evaluate_two_hop_satellite(iterations: int) -> dict:
    """Exercise aggregate second-pass support with an unsupported decoy."""

    sizes = (6, 4, 4, 4, 6)
    clusters = []
    species = []
    for size in sizes:
        start = len(species)
        clusters.append(list(range(start, start + size)))
        species.extend(range(start, start + size))

    queries = []
    targets = []
    scores = []

    def add_relation(source, target, norm, reverse_norm=None):
        reverse_norm = norm if reverse_norm is None else reverse_norm
        source_cluster = clusters[source]
        target_cluster = clusters[target]
        scale = float(np.sqrt(len(source_cluster) * len(target_cluster)))
        for offset in range(2):
            queries.append(source_cluster[offset])
            targets.append(target_cluster[offset])
            scores.append(norm * scale / 2.0)
            queries.append(target_cluster[offset])
            targets.append(source_cluster[offset])
            scores.append(reverse_norm * scale / 2.0)

    # Clusters 1 and 2 first attach to anchor 0. Cluster 3 has individually
    # weaker links to all three, but their aggregate support clears the 1.5
    # margin only after the first-pass component is rebuilt.
    add_relation(1, 0, 3.0)
    add_relation(2, 0, 3.0)
    add_relation(3, 0, 0.9)
    add_relation(3, 1, 0.9)
    add_relation(3, 2, 0.9)
    add_relation(3, 4, 1.0, reverse_norm=0.001)

    predicted, merges, relations, completed = (
        merge_supported_satellite_candidate_clusters(
            clusters,
            np.asarray(queries, dtype=np.int32),
            np.asarray(targets, dtype=np.int32),
            np.asarray(scores, dtype=np.float64),
            np.asarray(species, dtype=np.int32),
            max_component_genes=500,
            max_satellite_genes=12,
            max_satellite_to_anchor_ratio=0.75,
            min_margin=1.5,
            max_satellites_per_anchor=2,
            max_species_overlap_fraction=1.0,
            min_avg_score=0.0,
            min_max_score=0.0,
            min_coverage=0.5,
            min_norm=0.02,
            max_iterations=iterations,
        )
    )
    truth = [sum(clusters[:4], []), clusters[4]]
    return {
        "requested_iterations": iterations,
        "completed_iterations": completed,
        "seed_fragments": len(clusters),
        "candidate_families": len(predicted),
        "merges": merges,
        "directed_cluster_relations": relations,
        "pair_score": pair_score(predicted, truth),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json", type=Path)
    parser.add_argument("--seeds", type=int, nargs="+", default=[101, 211, 307])
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    results = []
    for seed in args.seeds:
        results.append({
            "seed": seed,
            "one_pass": evaluate(seed, 1),
            "bounded_four_pass": evaluate(seed, 4),
            "satellite_one_pass": evaluate_satellite(seed),
        })
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": (
            "Independent fragmented-family graph holdout with asymmetric "
            "cross-family decoys"
        ),
        "selection_rule": {
            "method": "reciprocal_best_cluster_normalized_support",
            "max_component_genes": 5,
            "min_average_score": 0.05,
            "min_maximum_score": 0.05,
            "max_iterations": 4,
        },
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "results": results,
        "two_hop_satellite": {
            "one_pass": evaluate_two_hop_satellite(1),
            "two_pass": evaluate_two_hop_satellite(2),
        },
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)
    for result in results:
        one = result["one_pass"]["pair_score"]
        bounded = result["bounded_four_pass"]["pair_score"]
        print(
            f"seed {result['seed']}: one-pass F={one['f_score']:.3f}; "
            f"bounded F={bounded['f_score']:.3f}"
        )
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
