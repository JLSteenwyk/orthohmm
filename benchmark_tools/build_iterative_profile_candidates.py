#!/usr/bin/env python3
"""Build phylogeny candidates with a second, self-calibrated profile-HMM pass.

This inference harness reads no benchmark labels. Accuracy scoring must be run
after its candidate and provenance files have been finalized.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
)
from orthohmm.accuracy import read_index_clusters
from orthohmm.files import fetch_fasta_files
from orthohmm.helpers import get_sequence_lengths
from orthohmm.metrics import PipelineMetrics
from orthohmm.search.profile_expansion import build_iterative_profile_candidates


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-directory", required=True, type=Path)
    parser.add_argument("--seed-clusters", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--matrix", default="BLOSUM62")
    parser.add_argument("--evalue", type=float, default=1e-4)
    parser.add_argument("--min-profile-cluster-size", type=int, default=3)
    parser.add_argument("--max-profile-cluster-size", type=int, default=200)
    parser.add_argument("--profile-kmer-k", type=int, default=4)
    parser.add_argument("--profile-reduced-alphabet", action="store_true")
    parser.add_argument("--profile-min-total-hits", type=int, default=4)
    parser.add_argument("--profile-min-diag-hits", type=int, default=1)
    parser.add_argument("--max-candidates-per-gene", type=int, default=20)
    parser.add_argument("--max-component-genes", type=int, default=500)
    parser.add_argument("--max-satellite-genes", type=int, default=12)
    parser.add_argument("--max-satellite-ratio", type=float, default=0.75)
    parser.add_argument("--min-margin", type=float, default=1.5)
    parser.add_argument("--max-satellites-per-anchor", type=int, default=4)
    parser.add_argument("--max-species-overlap", type=float, default=1.0)
    parser.add_argument("--min-coverage", type=float, default=0.5)
    parser.add_argument("--min-score-ratio", type=float, default=1.0)
    return parser


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _write_clusters(path: Path, clusters, gene_names: list[str]) -> None:
    with path.open("w") as handle:
        for cluster in clusters:
            handle.write(
                " ".join(sorted(gene_names[int(gene)] for gene in cluster)) + "\n"
            )


def _validate_partition(clusters, gene_count: int) -> None:
    members = [int(gene) for cluster in clusters for gene in cluster]
    if len(members) != gene_count or set(members) != set(range(gene_count)):
        raise ValueError("seed clusters are not an exact partition of FASTA genes")


def _gene_table(fasta_directory: Path, files: list[str]):
    gene_lengths = get_sequence_lengths(str(fasta_directory), files)
    order = np.argsort(
        np.asarray([str(row["name"]) for row in gene_lengths]), kind="stable"
    )
    gene_lengths = gene_lengths[order]
    gene_names = [str(row["name"]) for row in gene_lengths]
    if len(gene_names) != len(set(gene_names)):
        raise ValueError("FASTA gene IDs must be globally unique")
    species_names = sorted({str(row["spp"]) for row in gene_lengths})
    species_to_id = {species: index for index, species in enumerate(species_names)}
    gene_to_species = np.fromiter(
        (species_to_id[str(row["spp"])] for row in gene_lengths),
        dtype=np.int32,
        count=len(gene_lengths),
    )
    return gene_names, species_names, gene_to_species


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.cpu < 1:
        raise SystemExit("--cpu must be positive")
    fasta_directory = args.fasta_directory.resolve()
    seed_path = args.seed_clusters.resolve()
    output_directory = args.output_directory.resolve()
    result_json = args.json.resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise SystemExit(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)
    candidate_path = output_directory / "candidate_clusters.txt"
    trace_path = output_directory / "profile_candidate_merges.json"
    files = sorted(fetch_fasta_files(str(fasta_directory)))
    if not files:
        raise SystemExit(f"no FASTA files found: {fasta_directory}")

    parameters = {
        "matrix": args.matrix,
        "evalue_threshold": args.evalue,
        "min_profile_cluster_size": args.min_profile_cluster_size,
        "max_profile_cluster_size": args.max_profile_cluster_size,
        "profile_kmer_k": args.profile_kmer_k,
        "profile_use_reduced_alphabet": args.profile_reduced_alphabet,
        "profile_min_total_hits": args.profile_min_total_hits,
        "profile_min_diag_hits": args.profile_min_diag_hits,
        "max_candidates_per_gene": args.max_candidates_per_gene,
        "max_component_genes": args.max_component_genes,
        "max_satellite_genes": args.max_satellite_genes,
        "max_satellite_to_anchor_ratio": args.max_satellite_ratio,
        "min_margin": args.min_margin,
        "max_satellites_per_anchor": args.max_satellites_per_anchor,
        "max_species_overlap_fraction": args.max_species_overlap,
        "min_coverage": args.min_coverage,
        "min_score_ratio": args.min_score_ratio,
    }
    with PipelineMetrics(str(result_json)) as metrics:
        metrics.add_metadata(
            harness="benchmark_tools.build_iterative_profile_candidates",
            method="rebuilt_self_calibrated_profile_hmm_satellite_assignment",
            cpu_budget=args.cpu,
            parameters=parameters,
        )
        with metrics.stage("load_inputs"):
            gene_names, species_names, gene_to_species = _gene_table(
                fasta_directory, files
            )
            seed_clusters = read_index_clusters(str(seed_path), gene_names)
            _validate_partition(seed_clusters, len(gene_names))
        with metrics.stage("iterative_profile_candidates"):
            result = build_iterative_profile_candidates(
                seed_clusters,
                gene_names,
                str(fasta_directory),
                files,
                args.matrix,
                args.cpu,
                evalue_threshold=args.evalue,
                gene_to_species=gene_to_species,
                min_profile_cluster_size=args.min_profile_cluster_size,
                max_profile_cluster_size=args.max_profile_cluster_size,
                profile_kmer_k=args.profile_kmer_k,
                profile_use_reduced_alphabet=args.profile_reduced_alphabet,
                profile_min_total_hits=args.profile_min_total_hits,
                profile_min_diag_hits=args.profile_min_diag_hits,
                max_candidates_per_gene=args.max_candidates_per_gene,
                max_component_genes=args.max_component_genes,
                max_satellite_genes=args.max_satellite_genes,
                max_satellite_to_anchor_ratio=args.max_satellite_ratio,
                min_margin=args.min_margin,
                max_satellites_per_anchor=args.max_satellites_per_anchor,
                max_species_overlap_fraction=args.max_species_overlap,
                min_coverage=args.min_coverage,
                min_score_ratio=args.min_score_ratio,
            )
        with metrics.stage("write_outputs"):
            _write_clusters(candidate_path, result.clusters, gene_names)
            named_trace = []
            for record in result.merge_trace:
                named = dict(record)
                for field in ("source_genes", "target_genes"):
                    named[field] = [gene_names[int(gene)] for gene in record[field]]
                named_trace.append(named)
            _atomic_json(trace_path, named_trace)
        metrics.add_counts(
            genes=len(gene_names),
            species=len(species_names),
            seed_families=len(seed_clusters),
            profiles_built=result.profiles_built,
            profile_candidates=result.profile_candidates,
            significant_profile_hits=result.significant_profile_hits,
            strict_profile_hits=result.strict_profile_hits,
            directed_profile_relations=result.directed_relations,
            candidate_merges=result.merges,
            candidate_families=len(result.clusters),
        )

    payload = json.loads(result_json.read_text())
    payload.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git": git_state(),
        "inputs": {
            "seed_clusters": file_provenance(seed_path),
            "fasta_files": [
                file_provenance(fasta_directory / filename) for filename in files
            ],
        },
        "outputs": {
            "candidate_clusters": file_provenance(candidate_path),
            "merge_trace": file_provenance(trace_path),
        },
        "selection_isolation": (
            "No benchmark labels or scores were read while building candidates."
        ),
    })
    _atomic_json(result_json, payload)
    print(
        f"profile candidates: {len(seed_clusters)} -> {len(result.clusters)} "
        f"families ({result.merges} merges)"
    )
    print(result_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
