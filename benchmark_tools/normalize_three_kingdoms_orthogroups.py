#!/usr/bin/env python3
"""Normalize orthogroup formats used by the Three Kingdoms parity runs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from pathlib import Path


def iter_fastoma(path: Path) -> Iterator[tuple[str, tuple[str, ...]]]:
    """Yield groups from FastOMA's two-column OrthologousGroups.tsv."""
    groups: OrderedDict[str, list[str]] = OrderedDict()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Group", "Protein"]:
            raise ValueError(f"Unexpected FastOMA columns in {path}: {reader.fieldnames}")
        for line_number, row in enumerate(reader, start=2):
            group = row["Group"].strip()
            protein = row["Protein"].strip()
            if not group or not protein:
                raise ValueError(f"Empty FastOMA field at {path}:{line_number}")
            groups.setdefault(group, []).append(protein)
    for group, proteins in groups.items():
        yield group, tuple(proteins)


def iter_root_hogs(path: Path) -> Iterator[tuple[str, tuple[str, ...]]]:
    """Yield groups from OrthoHMM's root-HOG table."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"root_hog", "genes"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected root-HOG columns in {path}: {reader.fieldnames}")
        for line_number, row in enumerate(reader, start=2):
            group = row["root_hog"].strip()
            genes = tuple(gene for gene in row["genes"].split(",") if gene)
            if not group or not genes:
                raise ValueError(f"Empty root-HOG field at {path}:{line_number}")
            yield group, genes


_ORTHOMCL_MEMBER = re.compile(r"^(?P<gene>.+)\([^()]+\)$")


def iter_orthomcl(path: Path) -> Iterator[tuple[str, tuple[str, ...]]]:
    """Yield groups from OrthoMCL 1.4's annotated output."""
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            header, separator, payload = line.rstrip("\n").partition(":")
            if not separator or not header or not payload.strip():
                raise ValueError(f"Malformed OrthoMCL row at {path}:{line_number}")
            genes = []
            for member in payload.split():
                match = _ORTHOMCL_MEMBER.fullmatch(member)
                if match is None:
                    raise ValueError(
                        f"Malformed OrthoMCL member at {path}:{line_number}: {member}"
                    )
                genes.append(match.group("gene"))
            yield header.split("(", 1)[0], tuple(genes)


def write_groups(
    groups: Iterable[tuple[str, tuple[str, ...]]], output: Path
) -> tuple[int, int]:
    """Write one unlabeled, space-delimited group per line."""
    group_count = 0
    gene_count = 0
    observed: set[str] = set()
    with output.open("w") as handle:
        for group, genes in groups:
            duplicate = observed.intersection(genes)
            if duplicate:
                raise ValueError(
                    f"Genes occur in multiple groups; first duplicate in {group}: "
                    f"{min(duplicate)}"
                )
            observed.update(genes)
            handle.write(" ".join(genes) + "\n")
            group_count += 1
            gene_count += len(genes)
    return group_count, gene_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=("fastoma", "root-hogs", "orthomcl"))
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    readers = {
        "fastoma": iter_fastoma,
        "root-hogs": iter_root_hogs,
        "orthomcl": iter_orthomcl,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    group_count, gene_count = write_groups(readers[args.format](args.input), args.output)
    print(f"Wrote {group_count} groups containing {gene_count} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
