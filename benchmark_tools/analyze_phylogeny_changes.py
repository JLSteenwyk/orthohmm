#!/usr/bin/env python3
"""Stratify OrthoBench changes made by phylogenetic reconciliation.

This is a post-inference diagnostic: reference labels and sequences are read
only after both prediction files are frozen. Divergence strata use the median
mean pairwise identity of independently aligned reference orthogroups.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import itertools
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
from typing import Iterable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    load_refogs,
    read_clusters,
)


def load_sequences(
    fasta_directory: Path,
) -> tuple[dict[str, str], dict[str, str]]:
    """Return sequence and source-species maps keyed by FASTA identifier."""
    sequences: dict[str, str] = {}
    gene_species: dict[str, str] = {}
    for path in sorted(item for item in fasta_directory.iterdir() if item.is_file()):
        species = path.name
        header: str | None = None
        parts: list[str] = []
        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        if header in sequences:
                            raise ValueError(f"Duplicate FASTA identifier: {header}")
                        sequences[header] = "".join(parts)
                        gene_species[header] = species
                    header = line[1:].split()[0]
                    parts = []
                elif header is not None:
                    parts.append(line)
        if header is not None:
            if header in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {header}")
            sequences[header] = "".join(parts)
            gene_species[header] = species
    return sequences, gene_species


def read_root_hogs(path: Path) -> tuple[list[set[str]], list[str]]:
    clusters: list[set[str]] = []
    sources: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source_family", "genes"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                f"{path} must contain tab-separated source_family and genes columns"
            )
        for row in reader:
            genes = {gene for gene in row["genes"].split(",") if gene}
            if not genes:
                raise ValueError(f"Empty root HOG in {path}")
            clusters.append(genes)
            sources.append(row["source_family"])
    return clusters, sources


def _gene_cluster_index(clusters: Sequence[set[str]]) -> dict[str, int]:
    index: dict[str, int] = {}
    for cluster_id, cluster in enumerate(clusters):
        for gene in cluster:
            if gene in index:
                raise ValueError(f"Gene occurs in more than one cluster: {gene}")
            index[gene] = cluster_id
    return index


def relation_sets(
    clusters: Sequence[set[str]],
    gene_refogs: dict[str, set[int]],
) -> tuple[set[tuple[int, str, str]], set[tuple[str, str]]]:
    """Return true within-RefOG and false cross-RefOG predicted pairs."""
    within: set[tuple[int, str, str]] = set()
    cross: set[tuple[str, str]] = set()
    reference_genes = set(gene_refogs)
    for cluster in clusters:
        retained = sorted(cluster & reference_genes)
        for left, right in itertools.combinations(retained, 2):
            common_refogs = gene_refogs[left] & gene_refogs[right]
            if common_refogs:
                within.update((refog, left, right) for refog in common_refogs)
            else:
                cross.add((left, right))
    return within, cross


def pair_metrics(
    refogs: Sequence[set[str]],
    within: set[tuple[int, str, str]],
    cross: set[tuple[str, str]],
) -> dict[str, float | int]:
    possible = sum(len(group) * (len(group) - 1) // 2 for group in refogs)
    true_positive = len(within)
    precision = true_positive / (true_positive + len(cross)) if within or cross else 0.0
    recall = true_positive / possible if possible else 0.0
    f_score = (
        2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    )
    return {
        "true_pairs": true_positive,
        "cross_refog_pairs": len(cross),
        "possible_true_pairs": possible,
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
    }


def mean_pairwise_identity(aligned_sequences: Iterable[str]) -> float:
    records = list(aligned_sequences)
    if len(records) < 2:
        return 1.0
    identities = []
    for left, right in itertools.combinations(records, 2):
        if len(left) != len(right):
            raise ValueError("Aligned sequences have different lengths")
        comparable = [
            (left_aa, right_aa)
            for left_aa, right_aa in zip(left, right)
            if left_aa != "-" and right_aa != "-"
        ]
        if comparable:
            identities.append(
                sum(left_aa == right_aa for left_aa, right_aa in comparable)
                / len(comparable)
            )
    return statistics.fmean(identities) if identities else 0.0


def _read_fasta_records(path: Path) -> list[str]:
    records: list[str] = []
    parts: list[str] = []
    seen_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    records.append("".join(parts))
                seen_header = True
                parts = []
            elif seen_header:
                parts.append(line)
    if seen_header:
        records.append("".join(parts))
    return records


def align_refog(
    index: int,
    genes: set[str],
    sequences: dict[str, str],
    aligner: str,
    alignment_directory: Path,
) -> float:
    alignment = alignment_directory / f"RefOG{index + 1:03d}.faa"
    if not alignment.is_file() or alignment.stat().st_size == 0:
        missing = sorted(genes - sequences.keys())
        if missing:
            raise ValueError(
                f"RefOG{index + 1:03d} genes missing from FASTA input: "
                + ", ".join(missing[:5])
            )
        input_path = alignment.with_suffix(".input.faa")
        input_path.write_text(
            "".join(f">{gene}\n{sequences[gene]}\n" for gene in sorted(genes)),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [aligner, "--auto", "--thread", "1", str(input_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        temporary = alignment.with_suffix(".faa.tmp")
        temporary.write_text(completed.stdout, encoding="utf-8")
        temporary.replace(alignment)
        input_path.unlink()
    return mean_pairwise_identity(_read_fasta_records(alignment))


def split_summary(
    candidates: Sequence[set[str]],
    root_hogs: Sequence[set[str]],
    sources: Sequence[str],
    reference_genes: set[str],
) -> dict[str, int | bool]:
    if len(root_hogs) != len(sources):
        raise ValueError("Root HOG and source-family counts differ")
    source_counts: dict[int, int] = {}
    root_genes: set[str] = set()
    reference_split_sources = 0
    for root_hog, source_name in zip(root_hogs, sources):
        if not source_name.startswith("Family"):
            raise ValueError(f"Invalid source family identifier: {source_name}")
        source_index = int(source_name.removeprefix("Family"))
        if source_index >= len(candidates):
            raise ValueError(f"Unknown source family: {source_name}")
        if not root_hog <= candidates[source_index]:
            raise ValueError(
                f"{source_name} root HOG contains genes from another family"
            )
        if root_genes & root_hog:
            raise ValueError("A gene occurs in more than one root HOG")
        root_genes.update(root_hog)
        source_counts[source_index] = source_counts.get(source_index, 0) + 1
    split_sources = {index for index, count in source_counts.items() if count > 1}
    for index in split_sources:
        if candidates[index] & reference_genes:
            reference_split_sources += 1
    candidate_genes = set().union(*candidates) if candidates else set()
    return {
        "candidate_families": len(candidates),
        "root_hogs": len(root_hogs),
        "split_source_families": len(split_sources),
        "reference_bearing_split_source_families": reference_split_sources,
        "cross_source_merges": 0,
        "all_candidate_genes_preserved": root_genes == candidate_genes,
        "candidate_genes": len(candidate_genes),
        "root_hog_genes": len(root_genes),
    }


def _category_recall(
    refogs: Sequence[set[str]],
    within: set[tuple[int, str, str]],
    selected: Sequence[int],
) -> dict[str, float | int]:
    possible = 0
    recovered = 0
    for index in selected:
        pairs = {
            (index, left, right)
            for left, right in itertools.combinations(sorted(refogs[index]), 2)
        }
        possible += len(pairs)
        recovered += len(pairs & within)
    return {
        "refogs": len(selected),
        "possible_pairs": possible,
        "recovered_pairs": recovered,
        "recall": recovered / possible if possible else 0.0,
    }


def external_relation_count(
    group: set[str],
    clusters: Sequence[set[str]],
    gene_clusters: dict[str, int],
) -> int:
    """Count predicted relations from one RefOG to genes outside that RefOG."""
    cluster_ids = {gene_clusters[gene] for gene in group if gene in gene_clusters}
    relations = 0
    for cluster_id in cluster_ids:
        retained = len(clusters[cluster_id] & group)
        relations += retained * (len(clusters[cluster_id]) - retained)
    return relations


def reference_nonreference_relations(
    clusters: Sequence[set[str]], reference_genes: set[str]
) -> int:
    relations = 0
    for cluster in clusters:
        retained = len(cluster & reference_genes)
        relations += retained * (len(cluster) - retained)
    return relations


def _stratified_recall(
    refogs: Sequence[set[str]],
    baseline_within: set[tuple[int, str, str]],
    phylogeny_within: set[tuple[int, str, str]],
    categories: dict[str, list[int]],
) -> dict:
    result = {}
    for label, indices in categories.items():
        baseline = _category_recall(refogs, baseline_within, indices)
        phylogeny = _category_recall(refogs, phylogeny_within, indices)
        result[label] = {
            "baseline": baseline,
            "phylogeny": phylogeny,
            "recall_delta": phylogeny["recall"] - baseline["recall"],
        }
    return result


def _stratified_external_relations(
    refogs: Sequence[set[str]],
    candidates: Sequence[set[str]],
    root_hogs: Sequence[set[str]],
    categories: dict[str, list[int]],
) -> dict:
    candidate_index = _gene_cluster_index(candidates)
    root_hog_index = _gene_cluster_index(root_hogs)
    result = {}
    for label, indices in categories.items():
        baseline = sum(
            external_relation_count(refogs[index], candidates, candidate_index)
            for index in indices
        )
        phylogeny = sum(
            external_relation_count(refogs[index], root_hogs, root_hog_index)
            for index in indices
        )
        result[label] = {
            "baseline_external_relations": baseline,
            "phylogeny_external_relations": phylogeny,
            "removed_external_relations": baseline - phylogeny,
            "fraction_removed": (
                (baseline - phylogeny) / baseline if baseline else 0.0
            ),
        }
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--fasta-directory", required=True, type=Path)
    parser.add_argument("--baseline-clusters", required=True, type=Path)
    parser.add_argument("--phylogeny-root-hogs", required=True, type=Path)
    parser.add_argument("--aligner", default="mafft")
    parser.add_argument("--alignment-directory", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--json", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cpu < 1:
        raise SystemExit("--cpu must be at least 1")
    started = time.monotonic()
    refogs = load_refogs(args.refogs)
    reference_genes = set().union(*refogs)
    gene_refogs: dict[str, set[int]] = {}
    for index, group in enumerate(refogs):
        for gene in group:
            gene_refogs.setdefault(gene, set()).add(index)
    sequences, gene_species = load_sequences(args.fasta_directory)
    candidates = read_clusters(args.baseline_clusters)
    root_hogs, sources = read_root_hogs(args.phylogeny_root_hogs)
    candidate_index = _gene_cluster_index(candidates)
    root_hog_index = _gene_cluster_index(root_hogs)
    baseline_within, baseline_cross = relation_sets(candidates, gene_refogs)
    phylogeny_within, phylogeny_cross = relation_sets(root_hogs, gene_refogs)

    args.alignment_directory.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.cpu) as executor:
        identities = list(
            executor.map(
                lambda item: align_refog(
                    item[0],
                    item[1],
                    sequences,
                    args.aligner,
                    args.alignment_directory,
                ),
                enumerate(refogs),
            )
        )
    identity_cutoff = statistics.median(identities)
    copy_categories = {"single_copy": [], "multi_copy": []}
    identity_categories = {"lower_identity": [], "higher_identity": []}
    refog_records = []
    for index, (group, identity) in enumerate(zip(refogs, identities)):
        species_counts: dict[str, int] = {}
        for gene in group:
            species = gene_species[gene]
            species_counts[species] = species_counts.get(species, 0) + 1
        copy_label = (
            "multi_copy"
            if max(species_counts.values(), default=0) > 1
            else "single_copy"
        )
        identity_label = (
            "lower_identity" if identity <= identity_cutoff else "higher_identity"
        )
        copy_categories[copy_label].append(index)
        identity_categories[identity_label].append(index)
        possible = {
            (index, left, right)
            for left, right in itertools.combinations(sorted(group), 2)
        }
        refog_records.append(
            {
                "refog": f"RefOG{index + 1:03d}",
                "genes": len(group),
                "species": len(species_counts),
                "max_species_copies": max(species_counts.values(), default=0),
                "copy_category": copy_label,
                "mean_pairwise_identity": identity,
                "identity_category": identity_label,
                "baseline_pair_recall": len(possible & baseline_within) / len(possible),
                "phylogeny_pair_recall": len(possible & phylogeny_within)
                / len(possible),
                "baseline_external_relations": external_relation_count(
                    group, candidates, candidate_index
                ),
                "phylogeny_external_relations": external_relation_count(
                    group, root_hogs, root_hog_index
                ),
            }
        )

    removed_cross = baseline_cross - phylogeny_cross
    added_cross = phylogeny_cross - baseline_cross
    payload = {
        "schema_version": 1,
        "description": "Post-inference OrthoBench phylogeny change analysis",
        "command": [sys.executable, *sys.argv],
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "inputs": {
            "baseline_clusters": file_provenance(args.baseline_clusters),
            "phylogeny_root_hogs": file_provenance(args.phylogeny_root_hogs),
            "refog_directory": str(args.refogs.resolve()),
            "fasta_directory": str(args.fasta_directory.resolve()),
            "aligner": args.aligner,
        },
        "partition_changes": split_summary(
            candidates, root_hogs, sources, reference_genes
        ),
        "reference_pair_metrics": {
            "baseline": pair_metrics(refogs, baseline_within, baseline_cross),
            "phylogeny": pair_metrics(refogs, phylogeny_within, phylogeny_cross),
            "removed_cross_refog_pairs": len(removed_cross),
            "added_cross_refog_pairs": len(added_cross),
            "removed_true_pairs": len(baseline_within - phylogeny_within),
            "added_true_pairs": len(phylogeny_within - baseline_within),
            "baseline_reference_to_nonreference_relations": (
                reference_nonreference_relations(candidates, reference_genes)
            ),
            "phylogeny_reference_to_nonreference_relations": (
                reference_nonreference_relations(root_hogs, reference_genes)
            ),
        },
        "copy_number_strata": {
            "pair_recall": _stratified_recall(
                refogs, baseline_within, phylogeny_within, copy_categories
            ),
            "external_relations": _stratified_external_relations(
                refogs, candidates, root_hogs, copy_categories
            ),
        },
        "sequence_identity_strata": {
            "median_mean_pairwise_identity": identity_cutoff,
            "method": "MAFFT alignment; pairwise identity excludes gap-bearing positions",
            "strata": _stratified_recall(
                refogs, baseline_within, phylogeny_within, identity_categories
            ),
            "external_relations": _stratified_external_relations(
                refogs, candidates, root_hogs, identity_categories
            ),
        },
        "refogs": refog_records,
        "wall_s": time.monotonic() - started,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)
    print(
        json.dumps(
            {
                "partition_changes": payload["partition_changes"],
                "reference_pair_metrics": payload["reference_pair_metrics"],
                "copy_number_strata": payload["copy_number_strata"],
                "sequence_identity_strata": payload["sequence_identity_strata"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
