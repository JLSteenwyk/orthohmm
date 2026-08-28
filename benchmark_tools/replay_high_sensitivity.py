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
from orthohmm.refinement import refine_cluster_indices


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
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
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

    stages = []
    for label, path in (
        ("multipass", multipass_path),
        ("multipass_refined", refined_path),
    ):
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
            "profile_expansion": False,
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
