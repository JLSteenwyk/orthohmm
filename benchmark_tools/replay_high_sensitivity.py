#!/usr/bin/env python3
"""Replay production high-sensitivity graph inference from a trusted hit cache.

This offline validation tool never reads benchmark labels during inference. It
exists to test production graph and refinement code without repeating the
expensive all-to-all sequence search.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import resource
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    run_official_benchmark,
)
from orthohmm.accuracy import (
    build_rbnh_edges,
    build_singleton_assignment_edges,
    combine_edges,
    read_index_clusters,
)
from orthohmm.externals import execute_leiden
from orthohmm.files import fetch_fasta_files
from orthohmm.refinement import refine_cluster_indices
from orthohmm.search.profile_expansion import expand_profiles


def load_hits(path: Path):
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    required = {"all_gene_ids", "all_hits", "gene_to_species"}
    if not isinstance(payload, dict) or not required.issubset(payload):
        raise ValueError(f"{path} is not an OrthoHMM normalized-hit cache")
    return payload


def hit_arrays(all_hits, gene_to_id):
    count = len(all_hits)
    queries = np.fromiter(
        (gene_to_id[str(pair[0])] for pair in all_hits),
        dtype=np.int32,
        count=count,
    )
    targets = np.fromiter(
        (gene_to_id[str(pair[1])] for pair in all_hits),
        dtype=np.int32,
        count=count,
    )
    scores = np.fromiter(all_hits.values(), dtype=np.float64, count=count)
    return queries, targets, scores


def write_clusters(path: Path, clusters, gene_names) -> None:
    with path.open("w") as handle:
        for cluster in clusters:
            if cluster:
                handle.write(
                    " ".join(sorted(gene_names[int(gene)] for gene in cluster))
                    + "\n"
                )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-pickle", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--official-benchmark", type=Path)
    parser.add_argument("--cpm-resolution", type=float, default=0.1)
    parser.add_argument("--fasta-directory", type=Path)
    parser.add_argument("--matrix", default="BLOSUM62")
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--leiden-seed", type=int, default=4)
    parser.add_argument(
        "--profile-iterations",
        type=int,
        choices=(1, 2),
        default=1,
        help="Number of strict profile rebuild/search passes (default: 1)",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.profile_iterations > 1 and not args.fasta_directory:
        raise SystemExit(
            "--profile-iterations 2 requires --fasta-directory"
        )
    started = time.perf_counter()
    timings = {}
    args.output_directory.mkdir(parents=True, exist_ok=True)
    working = args.output_directory / "orthohmm_working_res"
    working.mkdir(parents=True, exist_ok=True)

    stage_started = time.perf_counter()
    payload = load_hits(args.hits_pickle)
    all_hits = payload["all_hits"]
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
    queries, targets, scores = hit_arrays(all_hits, gene_to_id)
    del payload, all_hits, gene_to_id
    timings["load_and_index_hits_s"] = round(
        time.perf_counter() - stage_started, 6
    )

    stage_started = time.perf_counter()
    rbnh_edges = build_rbnh_edges(
        gene_names, gene_to_species, queries, targets, scores
    )
    timings["rbnh_edges_s"] = round(time.perf_counter() - stage_started, 6)

    stage_started = time.perf_counter()
    execute_leiden(
        args.cpm_resolution,
        str(args.output_directory),
        edges=rbnh_edges,
        include_isolates=True,
        seed=args.leiden_seed,
    )
    clustered_path = working / "orthohmm_edges_clustered.txt"
    initial_clusters = read_index_clusters(str(clustered_path), gene_names)
    singleton_edges = build_singleton_assignment_edges(
        gene_names, initial_clusters, queries, targets, scores
    )
    multipass_edges = combine_edges(rbnh_edges, singleton_edges)
    execute_leiden(
        args.cpm_resolution,
        str(args.output_directory),
        edges=multipass_edges,
        include_isolates=True,
        seed=args.leiden_seed,
    )
    multipass_path = args.output_directory / "orthogroups_multipass.txt"
    shutil.copyfile(clustered_path, multipass_path)
    multipass_clusters = read_index_clusters(str(clustered_path), gene_names)
    timings["singleton_reassignment_s"] = round(
        time.perf_counter() - stage_started, 6
    )

    stage_started = time.perf_counter()
    refined_clusters = refine_cluster_indices(
        multipass_clusters,
        queries,
        targets,
        scores,
        multipass_edges.sources,
        multipass_edges.targets,
        gene_to_species,
        rbnh_scores=multipass_edges.weights,
    )
    refined_path = args.output_directory / "orthogroups_multipass_refined.txt"
    write_clusters(refined_path, refined_clusters, gene_names)
    timings["production_refinement_s"] = round(
        time.perf_counter() - stage_started, 6
    )

    profile_results = []
    profile_stage_paths = []
    profile_edges = None
    profile_base_edges = rbnh_edges
    final_singleton_edges = None
    if args.fasta_directory:
        stage_started = time.perf_counter()
        files = fetch_fasta_files(str(args.fasta_directory))
        profile_clusters = multipass_clusters
        for iteration in range(1, args.profile_iterations + 1):
            profile_result = expand_profiles(
                profile_clusters,
                gene_names,
                str(args.fasta_directory),
                files,
                args.matrix,
                args.cpu,
                1e-4,
                queries,
                targets,
                scores,
            )
            profile_results.append(profile_result)
            profile_base_edges = combine_edges(
                profile_base_edges, profile_result.edges
            )
            execute_leiden(
                args.cpm_resolution,
                str(args.output_directory),
                edges=profile_base_edges,
                include_isolates=True,
                seed=args.leiden_seed,
            )
            profile_base_clusters = read_index_clusters(
                str(clustered_path), gene_names
            )
            final_singleton_edges = build_singleton_assignment_edges(
                gene_names,
                profile_base_clusters,
                queries,
                targets,
                scores,
            )
            profile_edges = combine_edges(
                profile_base_edges, final_singleton_edges
            )
            execute_leiden(
                args.cpm_resolution,
                str(args.output_directory),
                edges=profile_edges,
                include_isolates=True,
                seed=args.leiden_seed,
            )
            profile_clusters = read_index_clusters(
                str(clustered_path), gene_names
            )
            suffix = "" if iteration == 1 else f"_iteration_{iteration}"
            profile_path = (
                args.output_directory / f"orthogroups_profiles{suffix}.txt"
            )
            shutil.copyfile(clustered_path, profile_path)
            refined_profile_clusters = refine_cluster_indices(
                profile_clusters,
                queries,
                targets,
                scores,
                profile_edges.sources,
                profile_edges.targets,
                gene_to_species,
                rbnh_scores=profile_edges.weights,
            )
            profile_refined_path = (
                args.output_directory
                / f"orthogroups_profiles{suffix}_refined.txt"
            )
            write_clusters(
                profile_refined_path, refined_profile_clusters, gene_names
            )
            label_suffix = "" if iteration == 1 else f"_iteration_{iteration}"
            profile_stage_paths.extend([
                (f"strict_profiles{label_suffix}", profile_path),
                (
                    f"strict_profiles{label_suffix}_refined",
                    profile_refined_path,
                ),
            ])
        timings["strict_profile_expansion_and_refinement_s"] = round(
            time.perf_counter() - stage_started, 6
        )

    stages = []
    stage_paths = [
        ("multipass", multipass_path),
        ("multipass_refined", refined_path),
    ]
    stage_paths.extend(profile_stage_paths)
    for label, path in stage_paths:
        record = {
            "label": label,
            "clusters": sum(1 for line in path.open() if line.strip()),
            "output": file_provenance(path),
        }
        if args.official_benchmark:
            record["official_orthobench"] = run_official_benchmark(
                args.official_benchmark, path
            )
        stages.append(record)

    max_rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "input": file_provenance(args.hits_pickle),
        "parameters": {
            "accuracy_profile": "high_sensitivity",
            "cpm_resolution": args.cpm_resolution,
            "profile_expansion": bool(profile_results),
            "profile_iterations": args.profile_iterations,
            "matrix": args.matrix,
            "leiden_seed": args.leiden_seed,
        },
        "counts": {
            "genes": len(gene_names),
            "species": len(species_names),
            "significant_hits": len(scores),
            "rbnh_edges": len(rbnh_edges),
            "singleton_assignment_edges": len(singleton_edges),
            "multipass_edges": len(multipass_edges),
        },
        "timings": timings,
        "wall_s": round(time.perf_counter() - started, 6),
        "peak_process_rss_gib": round(max_rss_kib / (1024 ** 2), 6),
        "stages": stages,
    }
    if profile_results:
        result["counts"].update({
            "profiles_built": sum(
                item.profiles_built for item in profile_results
            ),
            "profile_candidates": sum(
                item.profile_candidates for item in profile_results
            ),
            "significant_profile_hits": (
                sum(item.significant_profile_hits for item in profile_results)
            ),
            "strict_profile_edges": sum(
                len(item.edges) for item in profile_results
            ),
            "final_singleton_assignment_edges": len(final_singleton_edges),
            "profile_multipass_edges": len(profile_edges),
        })
        result["profile_iterations"] = [
            {
                "iteration": iteration,
                "profiles_built": item.profiles_built,
                "profile_candidates": item.profile_candidates,
                "significant_profile_hits": item.significant_profile_hits,
                "strict_profile_edges": len(item.edges),
            }
            for iteration, item in enumerate(profile_results, start=1)
        ]
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)

    for stage in stages:
        score = stage.get("official_orthobench", {})
        print(
            f"{stage['label']}: F={score.get('f_score', 'NA')} "
            f"P={score.get('precision', 'NA')} R={score.get('recall', 'NA')}"
        )
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
