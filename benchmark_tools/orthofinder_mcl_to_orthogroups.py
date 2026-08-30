#!/usr/bin/env python3
"""Convert OrthoFinder's pre-phylogenetic MCL clusters to orthogroups."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TextIO


def load_sequence_ids(path: Path) -> dict[str, str]:
    """Return OrthoFinder internal sequence ID to original FASTA token."""
    sequence_ids = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            internal_id, separator, title = line.rstrip("\n").partition(": ")
            if not separator or not internal_id or not title:
                raise ValueError(
                    f"Malformed sequence ID mapping at {path}:{line_number}"
                )
            if internal_id in sequence_ids:
                raise ValueError(f"Duplicate internal sequence ID: {internal_id}")
            sequence_ids[internal_id] = title.split(maxsplit=1)[0]
    return sequence_ids


def iter_mcl_clusters(path: Path) -> Iterator[tuple[str, ...]]:
    """Yield member IDs from an MCL matrix file with wrapped cluster rows."""
    in_matrix = False
    members: list[str] | None = None

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not in_matrix:
                if stripped == "begin":
                    in_matrix = True
                continue
            if stripped == ")":
                break
            if not stripped:
                continue

            tokens = stripped.split()
            if members is None:
                if not tokens[0].isdigit():
                    raise ValueError(
                        f"Expected numeric MCL row ID at {path}:{line_number}"
                    )
                members = []
                tokens = tokens[1:]

            row_finished = tokens[-1] == "$"
            if row_finished:
                tokens = tokens[:-1]
            elif tokens[-1].endswith("$"):
                tokens[-1] = tokens[-1][:-1]
                row_finished = True
            members.extend(token for token in tokens if token)

            if row_finished:
                if not members:
                    raise ValueError(f"Empty MCL cluster at {path}:{line_number}")
                yield tuple(members)
                members = None

    if not in_matrix:
        raise ValueError(f"MCL matrix has no begin marker: {path}")
    if members is not None:
        raise ValueError(f"Unterminated MCL cluster row: {path}")


def write_orthogroups(
    clusters: Iterable[tuple[str, ...]],
    sequence_ids: dict[str, str],
    output: TextIO,
) -> tuple[int, int]:
    """Write one space-delimited orthogroup per MCL cluster."""
    cluster_count = 0
    gene_count = 0
    observed = set()
    for cluster in clusters:
        genes = []
        for internal_id in cluster:
            try:
                gene = sequence_ids[internal_id]
            except KeyError as error:
                raise ValueError(
                    f"MCL cluster references unknown sequence ID: {internal_id}"
                ) from error
            if internal_id in observed:
                raise ValueError(f"Sequence occurs in multiple clusters: {internal_id}")
            observed.add(internal_id)
            genes.append(gene)
        output.write(" ".join(genes) + "\n")
        cluster_count += 1
        gene_count += len(genes)
    return cluster_count, gene_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "clusters",
        type=Path,
        help="OrthoFinder clusters_OrthoFinder_I*.txt_id_pairs.txt",
    )
    parser.add_argument("sequence_ids", type=Path, help="OrthoFinder SequenceIDs.txt")
    parser.add_argument("output", type=Path, help="Output orthogroups file")
    args = parser.parse_args()

    sequence_ids = load_sequence_ids(args.sequence_ids)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        cluster_count, gene_count = write_orthogroups(
            iter_mcl_clusters(args.clusters), sequence_ids, output
        )
    if gene_count != len(sequence_ids):
        raise ValueError(
            f"MCL clusters contain {gene_count} genes, expected {len(sequence_ids)}"
        )
    print(f"Wrote {cluster_count} clusters containing {gene_count} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
