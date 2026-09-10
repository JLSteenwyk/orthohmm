#!/usr/bin/env python3
"""Convert an OrthoMCL similarity matrix to canonical cross-species pairs."""

from __future__ import annotations

import argparse
import hashlib
import math
import re
import sys
from pathlib import Path
from typing import TextIO


DIMENSIONS_RE = re.compile(r"^dimensions\s+(\d+)x(\d+)$")


def _accession(identifier: str) -> str:
    fields = identifier.split("|")
    return fields[1] if len(fields) >= 2 and fields[1] else identifier


def load_species(gg_path: Path) -> dict[str, str]:
    """Return the OrthoMCL gene-to-species mapping from an all.gg file."""
    species_by_gene: dict[str, str] = {}
    with gg_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            taxon, separator, member_text = line.rstrip("\n").partition(":")
            members = member_text.split()
            if not separator or not taxon or not members:
                raise ValueError(f"Malformed OrthoMCL GG row at {gg_path}:{line_number}")
            for gene in members:
                if gene in species_by_gene:
                    raise ValueError(f"Duplicate OrthoMCL GG gene: {gene}")
                species_by_gene[gene] = taxon
    if not species_by_gene:
        raise ValueError(f"No genes found in {gg_path}")
    return species_by_gene


def load_index(
    index_path: Path, species_by_gene: dict[str, str]
) -> tuple[list[str], list[str]]:
    """Load the dense matrix index and its corresponding species labels."""
    genes: list[str] = []
    species: list[str] = []
    accessions: set[str] = set()
    with index_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"Malformed OrthoMCL index row at {index_path}:{line_number}"
                )
            try:
                index = int(fields[0])
            except ValueError as error:
                raise ValueError(
                    f"Invalid OrthoMCL index at {index_path}:{line_number}"
                ) from error
            if index != len(genes):
                raise ValueError(
                    f"Non-contiguous OrthoMCL index at {index_path}:{line_number}: "
                    f"expected {len(genes)}, found {index}"
                )
            gene = fields[1]
            if gene not in species_by_gene:
                raise ValueError(f"OrthoMCL index gene absent from GG file: {gene}")
            accession = _accession(gene)
            if accession in accessions:
                raise ValueError(f"Duplicate normalized OrthoMCL accession: {accession}")
            accessions.add(accession)
            genes.append(accession)
            species.append(species_by_gene[gene])
    if not genes:
        raise ValueError(f"No genes found in {index_path}")
    return genes, species


def _edge_fingerprint(source: int, target: int, score: str) -> int:
    low, high = sorted((source, target))
    digest = hashlib.blake2b(
        f"{low}\t{high}\t{score}".encode(), digest_size=16
    ).digest()
    return int.from_bytes(digest, "big")


def write_pairs(
    matrix_path: Path,
    index_path: Path,
    gg_path: Path,
    output: TextIO,
) -> dict[str, int]:
    """Validate a native OrthoMCL graph and write each cross-species edge once."""
    species_by_gene = load_species(gg_path)
    genes, species = load_index(index_path, species_by_gene)
    dimension: int | None = None
    in_matrix = False
    rows_seen = bytearray(len(genes))
    row_count = 0
    directed_edges = 0
    same_species_edges = 0
    pair_count = 0
    reciprocal_fingerprint = 0

    with matrix_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            match = DIMENSIONS_RE.match(stripped)
            if match:
                rows, columns = map(int, match.groups())
                if rows != columns:
                    raise ValueError(f"Non-square OrthoMCL matrix in {matrix_path}")
                if dimension is not None:
                    raise ValueError(f"Duplicate dimensions declaration in {matrix_path}")
                dimension = rows
                continue
            if stripped == "begin":
                if in_matrix:
                    raise ValueError(f"Duplicate matrix begin marker in {matrix_path}")
                in_matrix = True
                continue
            if not in_matrix or not stripped or stripped.startswith("("):
                continue
            if stripped == ")":
                in_matrix = False
                continue
            if not stripped.endswith("$"):
                raise ValueError(
                    f"Malformed OrthoMCL matrix row at {matrix_path}:{line_number}"
                )
            fields = stripped[:-1].split()
            if not fields:
                raise ValueError(f"Empty OrthoMCL matrix row at {matrix_path}:{line_number}")
            try:
                source = int(fields[0])
            except ValueError as error:
                raise ValueError(
                    f"Invalid OrthoMCL matrix row at {matrix_path}:{line_number}"
                ) from error
            if not 0 <= source < len(genes):
                raise ValueError(
                    f"OrthoMCL matrix source out of range at {matrix_path}:{line_number}"
                )
            if rows_seen[source]:
                raise ValueError(f"Duplicate OrthoMCL matrix row: {source}")
            rows_seen[source] = 1
            row_count += 1
            targets_seen: set[int] = set()

            for edge in fields[1:]:
                target_text, separator, score = edge.partition(":")
                if not separator or not target_text or not score:
                    raise ValueError(
                        f"Malformed OrthoMCL edge at {matrix_path}:{line_number}: {edge}"
                    )
                try:
                    target = int(target_text)
                    numeric_score = float(score)
                except ValueError as error:
                    raise ValueError(
                        f"Invalid OrthoMCL edge at {matrix_path}:{line_number}: {edge}"
                    ) from error
                if not 0 <= target < len(genes):
                    raise ValueError(
                        f"OrthoMCL matrix target out of range at "
                        f"{matrix_path}:{line_number}"
                    )
                if source == target:
                    raise ValueError(f"Self-edge in OrthoMCL matrix: {source}")
                if target in targets_seen:
                    raise ValueError(
                        f"Duplicate OrthoMCL edge in row {source}: target {target}"
                    )
                if not math.isfinite(numeric_score) or numeric_score < 0:
                    raise ValueError(
                        f"Invalid OrthoMCL edge score at {matrix_path}:{line_number}: "
                        f"{score}"
                    )
                targets_seen.add(target)
                directed_edges += 1
                reciprocal_fingerprint ^= _edge_fingerprint(source, target, score)
                if species[source] == species[target]:
                    same_species_edges += 1
                    continue
                if source < target:
                    first, second = sorted((genes[source], genes[target]))
                    if first == second:
                        raise ValueError(f"Self-pair after accession normalization: {first}")
                    output.write(f"{first}\t{second}\n")
                    pair_count += 1

    if dimension is None:
        raise ValueError(f"Missing dimensions declaration in {matrix_path}")
    if dimension != len(genes):
        raise ValueError(
            f"OrthoMCL matrix/index dimension mismatch: {dimension} != {len(genes)}"
        )
    if row_count != dimension or any(value == 0 for value in rows_seen):
        raise ValueError(
            f"Incomplete OrthoMCL matrix: found {row_count} of {dimension} rows"
        )
    if reciprocal_fingerprint:
        raise ValueError("OrthoMCL matrix edges or scores are not reciprocal")
    if directed_edges % 2 or same_species_edges % 2:
        raise ValueError("OrthoMCL matrix has an odd directed-edge count")
    expected_pairs = (directed_edges - same_species_edges) // 2
    if pair_count != expected_pairs:
        raise ValueError(
            f"OrthoMCL cross-species edge mismatch: {pair_count} != {expected_pairs}"
        )
    return {
        "matrix_genes": dimension,
        "directed_edges": directed_edges,
        "same_species_directed_edges": same_species_edges,
        "cross_species_pairs": pair_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("index", type=Path)
    parser.add_argument("gg", type=Path)
    args = parser.parse_args()
    stats = write_pairs(args.matrix, args.index, args.gg, sys.stdout)
    print(
        f"Converted {stats['cross_species_pairs']:,} native cross-species pairs "
        f"from {stats['matrix_genes']:,} matrix genes and "
        f"{stats['directed_edges']:,} validated directed edges; skipped "
        f"{stats['same_species_directed_edges'] // 2:,} within-species pairs",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
