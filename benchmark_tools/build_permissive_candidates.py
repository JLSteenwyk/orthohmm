#!/usr/bin/env python3
"""Build label-blind candidate superfamilies from cached OrthoHMM hits.

This development harness varies only the reciprocal-best normalized-hit
support factor and MCL inflation. Reference labels are never read during graph
construction or clustering; an optional OrthoBench score is computed only
after the candidate file has been finalized.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
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
from orthohmm.accuracy import build_rbnh_edges
from orthohmm.metrics import PipelineMetrics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-pickle", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--mcl", default="mcl")
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--threshold-factor", type=float, default=1.0)
    parser.add_argument("--inflation", type=float, default=1.2)
    parser.add_argument("--official-benchmark", type=Path)
    return parser


def write_binary_abc(path: Path, sources, targets) -> None:
    """Write integer-ID MCL input with constant edge weights."""
    with path.open("w") as handle:
        for source, target in zip(sources, targets):
            handle.write(f"{int(source)}\t{int(target)}\t1\n")


def materialize_clusters(
    mcl_path: Path,
    output_path: Path,
    gene_names: list[str],
) -> int:
    """Map MCL integer IDs to genes and restore omitted isolates."""
    clusters: list[list[int]] = []
    seen: set[int] = set()
    with mcl_path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                members = sorted({int(value) for value in line.split()})
            except ValueError as exc:
                raise ValueError(
                    f"{mcl_path}:{line_number}: non-integer cluster member"
                ) from exc
            if not members:
                continue
            if members[0] < 0 or members[-1] >= len(gene_names):
                raise ValueError(
                    f"{mcl_path}:{line_number}: gene ID outside input table"
                )
            duplicate = seen.intersection(members)
            if duplicate:
                raise ValueError(
                    f"{mcl_path}:{line_number}: repeated gene IDs across clusters"
                )
            seen.update(members)
            clusters.append(members)
    clusters.extend(
        [[gene_id] for gene_id in range(len(gene_names)) if gene_id not in seen]
    )
    clusters.sort(key=lambda cluster: (cluster[0], len(cluster), cluster))
    with output_path.open("w") as handle:
        for cluster in clusters:
            handle.write(" ".join(gene_names[gene] for gene in cluster) + "\n")
    return len(clusters)


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.cpu < 1:
        raise SystemExit("--cpu must be positive")
    if not np.isfinite(args.inflation) or args.inflation <= 1.0:
        raise SystemExit("--inflation must be finite and greater than 1")
    if not np.isfinite(args.threshold_factor) or args.threshold_factor <= 0.0:
        raise SystemExit("--threshold-factor must be finite and positive")
    mcl = shutil.which(args.mcl) if os.sep not in args.mcl else args.mcl
    if not mcl or not Path(mcl).is_file():
        raise SystemExit(f"MCL executable not found: {args.mcl}")

    output_directory = args.output_directory.resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise SystemExit(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)
    edges_path = output_directory / "edges.binary.abc"
    raw_clusters = output_directory / "mcl_clusters.txt"
    candidate_clusters = output_directory / "candidate_clusters.txt"
    log_path = output_directory / "mcl.log"
    result_json = args.json.resolve()

    with PipelineMetrics(str(result_json)) as metrics:
        metrics.add_metadata(
            harness="benchmark_tools.build_permissive_candidates",
            cpu_budget=args.cpu,
        )
        with metrics.stage("load_hits"):
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
        with metrics.stage("relative_support_graph"):
            edges = build_rbnh_edges(
                gene_names,
                gene_to_species,
                queries,
                targets,
                scores,
                threshold_factor=args.threshold_factor,
            )
        with metrics.stage("write_binary_graph"):
            write_binary_abc(edges_path, edges.sources, edges.targets)
        command = [
            str(mcl),
            str(edges_path),
            "--abc",
            "-I",
            str(args.inflation),
            "-te",
            str(args.cpu),
            "-o",
            str(raw_clusters),
        ]
        with metrics.stage("mcl"):
            with log_path.open("w") as log:
                subprocess.run(
                    command,
                    check=True,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
        with metrics.stage("materialize_candidates"):
            cluster_count = materialize_clusters(
                raw_clusters, candidate_clusters, gene_names
            )
        metrics.add_counts(
            genes=len(gene_names),
            species=len(species_names),
            significant_hits=len(scores),
            candidate_edges=len(edges),
            candidate_families=cluster_count,
        )

    payload = json.loads(result_json.read_text())
    payload.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "input": file_provenance(args.hits_pickle),
        "parameters": {
            "threshold_factor": args.threshold_factor,
            "edge_weights": "binary",
            "inflation": args.inflation,
            "mcl": str(Path(mcl).resolve()),
            "cpu": args.cpu,
        },
        "mcl_command": command,
        "outputs": {
            "candidate_clusters": file_provenance(candidate_clusters),
            "binary_edges": file_provenance(edges_path),
            "mcl_log": file_provenance(log_path),
        },
    })
    if args.official_benchmark is not None:
        payload["official_orthobench"] = run_official_benchmark(
            args.official_benchmark.resolve(), candidate_clusters
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
