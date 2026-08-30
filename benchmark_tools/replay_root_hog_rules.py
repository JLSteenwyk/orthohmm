#!/usr/bin/env python3
"""Replay root-HOG evidence rules from fixed rooted gene trees."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    run_official_benchmark,
)
from orthohmm.files import fetch_fasta_files
from orthohmm.phylogeny import parse_gene_tree, parse_species_tree, reconcile_gene_tree
from orthohmm.phylogeny_pipeline import (
    _load_sequence_data,
    family_requires_reconciliation,
)


RULES = ("supported_children", "species_overlap", "mapped_event")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-directory", required=True, type=Path)
    parser.add_argument("--candidate-clusters", required=True, type=Path)
    parser.add_argument("--phylogeny-output", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--official-benchmark", type=Path)
    parser.add_argument("--rules", nargs="+", choices=RULES, default=list(RULES))
    return parser


def _read_clusters(path: Path) -> list[tuple[str, ...]]:
    return [
        tuple(sorted(line.split()))
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def write_groups(path: Path, groups) -> None:
    with path.open("w") as handle:
        for group in groups:
            handle.write(" ".join(group) + "\n")


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    output_directory = args.output_directory.resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise SystemExit(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)
    phylogeny_output = args.phylogeny_output.resolve()
    species_tree_path = phylogeny_output / "species_tree.rooted.nwk"
    fasta_directory = args.fasta_directory.resolve()
    files = fetch_fasta_files(str(fasta_directory))
    sequences, gene_species = _load_sequence_data(
        str(fasta_directory), files
    )
    del sequences
    species_tree = parse_species_tree(species_tree_path, files).tree
    clusters = _read_clusters(args.candidate_clusters)

    outputs = {}
    for rule in args.rules:
        started = time.perf_counter()
        groups = []
        reconciled = 0
        for family_index, genes in enumerate(clusters):
            family_id = f"Family{family_index:07d}"
            if not family_requires_reconciliation(genes, gene_species):
                groups.append(genes)
                continue
            tree_path = (
                phylogeny_output / "gene_trees" / f"{family_id}.rooted.nwk"
            )
            if not tree_path.is_file():
                raise SystemExit(f"rooted gene tree not found: {tree_path}")
            tree = parse_gene_tree(tree_path.read_text())
            result = reconcile_gene_tree(
                tree,
                species_tree,
                {gene: gene_species[gene] for gene in genes},
                family_id=family_id,
                root_duplication_rule=rule,
            )
            groups.extend(result.root_groups)
            reconciled += 1
        output_path = output_directory / f"root_hogs_{rule}.txt"
        write_groups(output_path, groups)
        record = {
            "rule": rule,
            "candidate_families": len(clusters),
            "reconciled_families": reconciled,
            "root_hogs": len(groups),
            "wall_s": round(time.perf_counter() - started, 6),
            "output": file_provenance(output_path),
        }
        if args.official_benchmark is not None:
            record["official_orthobench"] = run_official_benchmark(
                args.official_benchmark.resolve(), output_path
            )
        outputs[rule] = record

    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": "Root-HOG rule replay from immutable rooted gene trees",
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "inputs": {
            "candidate_clusters": file_provenance(args.candidate_clusters),
            "species_tree": file_provenance(species_tree_path),
            "phylogeny_output": str(phylogeny_output),
        },
        "rules": outputs,
    }
    _atomic_json(args.json.resolve(), payload)
    for rule, record in outputs.items():
        score = record.get("official_orthobench", {})
        print(
            f"{rule}: F={score.get('f_score', 'NA')} "
            f"P={score.get('precision', 'NA')} "
            f"R={score.get('recall', 'NA')}"
        )
    print(args.json.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
