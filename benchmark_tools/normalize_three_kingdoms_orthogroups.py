#!/usr/bin/env python3
"""Normalize orthogroup formats used by the Three Kingdoms parity runs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from pathlib import Path


_SONICPARANOID_METADATA = (
    "group_id",
    "group_size",
    "sp_in_grp",
    "seed_ortholog_cnt",
)


def iter_sonicparanoid(path: Path) -> Iterator[tuple[str, tuple[str, ...]]]:
    """Yield groups from SonicParanoid's species-column output."""
    observed_groups: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or tuple(reader.fieldnames[:4]) != _SONICPARANOID_METADATA:
            raise ValueError(
                f"Unexpected SonicParanoid columns in {path}: {reader.fieldnames}"
            )
        species_columns = reader.fieldnames[4:]
        if not species_columns:
            raise ValueError(f"No species columns in SonicParanoid table: {path}")

        for line_number, row in enumerate(reader, start=2):
            group = row["group_id"].strip()
            if not group or group in observed_groups:
                raise ValueError(
                    f"Empty or duplicate SonicParanoid group at {path}:{line_number}: "
                    f"{group!r}"
                )
            observed_groups.add(group)

            genes: list[str] = []
            occupied_species = 0
            for column in species_columns:
                cell = row[column]
                if cell is None:
                    raise ValueError(f"Truncated SonicParanoid row at {path}:{line_number}")
                cell = cell.strip()
                if not cell or cell == "*":
                    continue
                occupied_species += 1
                genes.extend(gene.strip() for gene in cell.split(",") if gene.strip())

            try:
                expected_size = int(row["group_size"])
                expected_species = int(row["sp_in_grp"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid SonicParanoid counts at {path}:{line_number}"
                ) from error
            if len(genes) != expected_size or occupied_species != expected_species:
                raise ValueError(
                    f"SonicParanoid count mismatch at {path}:{line_number}: "
                    f"found {len(genes)} genes in {occupied_species} species, expected "
                    f"{expected_size} genes in {expected_species} species"
                )
            yield group, tuple(genes)


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
    parser.add_argument(
        "format", choices=("fastoma", "root-hogs", "orthomcl", "sonicparanoid")
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    readers = {
        "fastoma": iter_fastoma,
        "root-hogs": iter_root_hogs,
        "orthomcl": iter_orthomcl,
        "sonicparanoid": iter_sonicparanoid,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    group_count, gene_count = write_groups(readers[args.format](args.input), args.output)
    print(f"Wrote {group_count} groups containing {gene_count} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
