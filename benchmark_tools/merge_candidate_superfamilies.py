#!/usr/bin/env python3
"""Merge precise OrthoHMM groups into label-blind candidate superfamilies.

All significant HMM hits are aggregated between seed groups. Each group can
join only its mutually best-supported partner, subject to generic evidence
gates and a family-size cap. Optional OrthoBench scoring occurs only after the
candidate partition has been finalized.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    run_official_benchmark,
)
from benchmark_tools.replay_high_sensitivity import hit_arrays, load_hits
from orthohmm.accuracy import read_index_clusters
from orthohmm.metrics import PipelineMetrics
from orthohmm.refinement import (
    merge_reciprocal_candidate_clusters,
    merge_supported_satellite_candidate_clusters,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-pickle", required=True, type=Path)
    parser.add_argument("--seed-clusters", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--max-family-genes", type=int, default=500)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument(
        "--method", choices=("reciprocal", "satellite"), default="reciprocal"
    )
    parser.add_argument("--max-satellite-genes", type=int, default=5)
    parser.add_argument("--max-satellite-ratio", type=float, default=0.75)
    parser.add_argument("--min-margin", type=float, default=1.5)
    parser.add_argument("--iteration-margin-increment", type=float, default=0.0)
    parser.add_argument("--max-satellites-per-anchor", type=int, default=2)
    parser.add_argument("--max-species-overlap-fraction", type=float, default=1.0)
    parser.add_argument("--min-average-score", type=float, default=0.0)
    parser.add_argument("--min-maximum-score", type=float, default=0.0)
    parser.add_argument("--min-coverage", type=float, default=0.0)
    parser.add_argument("--min-normalized-support", type=float, default=0.0)
    parser.add_argument("--trace-json", type=Path)
    parser.add_argument("--official-benchmark", type=Path)
    return parser


def validate_partition(clusters, gene_count: int) -> None:
    members = np.fromiter(
        (int(gene) for cluster in clusters for gene in cluster),
        dtype=np.int64,
    )
    if len(members) != gene_count:
        raise ValueError(
            f"seed clusters contain {len(members)} memberships for "
            f"{gene_count} genes"
        )
    if gene_count and not np.array_equal(np.sort(members), np.arange(gene_count)):
        raise ValueError("seed clusters must contain every gene exactly once")


def write_clusters(path: Path, clusters, gene_names: list[str]) -> None:
    ordered = [sorted(int(gene) for gene in cluster) for cluster in clusters]
    ordered.sort(key=lambda cluster: (cluster[0], len(cluster), cluster))
    with path.open("w") as handle:
        for cluster in ordered:
            handle.write(" ".join(gene_names[gene] for gene in cluster) + "\n")


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    gates = (
        args.min_average_score,
        args.min_maximum_score,
        args.min_coverage,
        args.min_normalized_support,
        args.iteration_margin_increment,
    )
    if args.max_family_genes < 1:
        raise SystemExit("--max-family-genes must be positive")
    if args.iterations < 1:
        raise SystemExit("--iterations must be positive")
    if any(not np.isfinite(value) or value < 0.0 for value in gates):
        raise SystemExit("evidence gates must be finite and nonnegative")
    if args.trace_json is not None and args.method != "satellite":
        raise SystemExit("--trace-json requires --method satellite")

    output_directory = args.output_directory.resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise SystemExit(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)
    output_path = output_directory / "candidate_clusters.txt"
    result_json = args.json.resolve()

    with PipelineMetrics(str(result_json)) as metrics:
        metrics.add_metadata(
            harness="benchmark_tools.merge_candidate_superfamilies"
        )
        with metrics.stage("load_hits_and_seeds"):
            payload = load_hits(args.hits_pickle)
            gene_names = sorted({str(gene) for gene in payload["all_gene_ids"]})
            gene_to_id = {gene: idx for idx, gene in enumerate(gene_names)}
            species_names = sorted({
                str(payload["gene_to_species"][gene]) for gene in gene_names
            })
            species_to_id = {
                species: idx for idx, species in enumerate(species_names)
            }
            gene_to_species = np.fromiter(
                (
                    species_to_id[str(payload["gene_to_species"][gene])]
                    for gene in gene_names
                ),
                dtype=np.int32,
                count=len(gene_names),
            )
            queries, targets, scores = hit_arrays(
                payload["all_hits"], gene_to_id
            )
            del payload, gene_to_id
            seed_clusters = read_index_clusters(
                str(args.seed_clusters), gene_names
            )
            validate_partition(seed_clusters, len(gene_names))
        with metrics.stage("reciprocal_candidate_merges"):
            merge_function = (
                merge_supported_satellite_candidate_clusters
                if args.method == "satellite"
                else merge_reciprocal_candidate_clusters
            )
            method_parameters = {}
            merge_trace = [] if args.trace_json is not None else None
            if args.method == "satellite":
                method_parameters = {
                    "max_satellite_genes": args.max_satellite_genes,
                    "max_satellite_to_anchor_ratio": args.max_satellite_ratio,
                    "min_margin": args.min_margin,
                    "iteration_margin_increment": (
                        args.iteration_margin_increment
                    ),
                    "max_satellites_per_anchor": args.max_satellites_per_anchor,
                    "max_species_overlap_fraction": (
                        args.max_species_overlap_fraction
                    ),
                    "merge_trace": merge_trace,
                }
            clusters, merge_count, relation_count, completed_iterations = (
                merge_function(
                    seed_clusters,
                    queries,
                    targets,
                    scores,
                    gene_to_species,
                    max_component_genes=args.max_family_genes,
                    min_avg_score=args.min_average_score,
                    min_max_score=args.min_maximum_score,
                    min_coverage=args.min_coverage,
                    min_norm=args.min_normalized_support,
                    max_iterations=args.iterations,
                    **method_parameters,
                )
            )
        with metrics.stage("write_candidates"):
            write_clusters(output_path, clusters, gene_names)
            if args.trace_json is not None:
                trace_payload = {
                    "schema_version": 1,
                    "description": "Accepted label-blind satellite merge evidence",
                    "command": [sys.executable, *sys.argv],
                    "git": git_state(),
                    "input": {
                        "hits": file_provenance(args.hits_pickle),
                        "seed_clusters": file_provenance(args.seed_clusters),
                    },
                    "merges": [],
                }
                for record in merge_trace or []:
                    named_record = dict(record)
                    for field in ("source_genes", "target_genes"):
                        named_record[field] = [
                            gene_names[int(gene)] for gene in record[field]
                        ]
                    trace_payload["merges"].append(named_record)
                _atomic_json(args.trace_json.resolve(), trace_payload)
        metrics.add_counts(
            genes=len(gene_names),
            species=len(species_names),
            significant_hits=len(scores),
            seed_families=len(seed_clusters),
            cross_family_directed_relations=relation_count,
            reciprocal_family_merges=merge_count,
            completed_merge_iterations=completed_iterations,
            candidate_families=len(clusters),
        )

    payload = json.loads(result_json.read_text())
    payload.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "input": {
            "hits": file_provenance(args.hits_pickle),
            "seed_clusters": file_provenance(args.seed_clusters),
        },
        "parameters": {
            "selection": args.method,
            "max_family_genes": args.max_family_genes,
            "max_iterations": args.iterations,
            "min_average_score": args.min_average_score,
            "min_maximum_score": args.min_maximum_score,
            "min_coverage": args.min_coverage,
            "min_normalized_support": args.min_normalized_support,
            "max_satellite_genes": args.max_satellite_genes,
            "max_satellite_ratio": args.max_satellite_ratio,
            "min_margin": args.min_margin,
            "iteration_margin_increment": args.iteration_margin_increment,
            "max_satellites_per_anchor": args.max_satellites_per_anchor,
            "max_species_overlap_fraction": args.max_species_overlap_fraction,
        },
        "outputs": {"candidate_clusters": file_provenance(output_path)},
    })
    if args.trace_json is not None:
        payload["outputs"]["merge_trace"] = file_provenance(
            args.trace_json.resolve()
        )
    if args.official_benchmark is not None:
        payload["official_orthobench"] = run_official_benchmark(
            args.official_benchmark.resolve(), output_path
        )
    _atomic_json(result_json, payload)
    score = payload.get("official_orthobench", {})
    print(
        f"candidates: F={score.get('f_score', 'NA')} "
        f"P={score.get('precision', 'NA')} R={score.get('recall', 'NA')}"
    )
    print(result_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
