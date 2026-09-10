#!/usr/bin/env python3
"""Convert SonicParanoid's native species-pair tables to QfO pairs."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterator
from itertools import product
from pathlib import Path
from typing import TextIO


EXPECTED_HEADER = ["Size", "Relations", "OrthoA", "OrthoB"]


def _accession(identifier: str) -> str:
    fields = identifier.split("|")
    return fields[1] if len(fields) >= 2 and fields[1] else identifier


def _members(field: str, path: Path, line_number: int) -> tuple[str, ...]:
    tokens = field.split()
    if not tokens or len(tokens) % 2:
        raise ValueError(f"Malformed SonicParanoid members at {path}:{line_number}")
    members = tuple(tokens[::2])
    if len(members) != len(set(members)):
        raise ValueError(f"Duplicate SonicParanoid member at {path}:{line_number}")
    try:
        for score in tokens[1::2]:
            float(score)
    except ValueError as error:
        raise ValueError(
            f"Invalid SonicParanoid confidence at {path}:{line_number}"
        ) from error
    return members


def iter_pairs(
    path: Path, stats: dict[str, int] | None = None
) -> Iterator[tuple[str, str]]:
    """Yield validated ortholog pairs from one species-pair table."""
    observed: set[tuple[str, str]] = set()
    with path.open() as handle:
        header = next(handle, "").rstrip("\n").split("\t")
        if header != EXPECTED_HEADER:
            raise ValueError(f"Unexpected SonicParanoid pair columns in {path}: {header}")
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Malformed SonicParanoid row at {path}:{line_number}")
            try:
                expected_size = int(fields[0])
                expected_relations = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"Invalid SonicParanoid counts at {path}:{line_number}"
                ) from error
            left = _members(fields[2], path, line_number)
            right = _members(fields[3], path, line_number)
            if expected_size != len(left) + len(right):
                raise ValueError(f"SonicParanoid size mismatch at {path}:{line_number}")
            if expected_relations != len(left) * len(right):
                raise ValueError(
                    f"SonicParanoid relation mismatch at {path}:{line_number}"
                )
            for first, second in product(left, right):
                pair = tuple(sorted((_accession(first), _accession(second))))
                if pair[0] == pair[1]:
                    raise ValueError(f"Self-pair in SonicParanoid output: {pair[0]}")
                if pair in observed:
                    if stats is not None:
                        stats["duplicates"] = stats.get("duplicates", 0) + 1
                    continue
                observed.add(pair)
                yield pair


def find_pair_directory(path: Path) -> Path:
    if path.name == "species_to_species_orthologs" and path.is_dir():
        return path
    candidates = sorted(path.glob("runs/*/species_to_species_orthologs"))
    if len(candidates) != 1:
        raise ValueError(
            f"Expected one species_to_species_orthologs directory under {path}, "
            f"found {len(candidates)}"
        )
    return candidates[0]


def write_pairs(pair_directory: Path, output: TextIO) -> tuple[int, int, int]:
    files = sorted(path for path in pair_directory.rglob("*") if path.is_file())
    if not files:
        raise ValueError(f"No SonicParanoid pair tables found in {pair_directory}")
    pair_count = 0
    stats = {"duplicates": 0}
    observed_species_pairs: set[tuple[str, str]] = set()
    observed_species: set[str] = set()
    for path in files:
        left_species = path.parent.name
        prefix = f"{left_species}-"
        if not path.name.startswith(prefix):
            raise ValueError(
                f"Unexpected SonicParanoid pair-table path under {pair_directory}: {path}"
            )
        species_pair = tuple(sorted((left_species, path.name[len(prefix) :])))
        if not all(species_pair) or species_pair[0] == species_pair[1]:
            raise ValueError(f"Invalid SonicParanoid species pair: {path}")
        if species_pair in observed_species_pairs:
            raise ValueError(f"Duplicate SonicParanoid species-pair table: {species_pair}")
        observed_species_pairs.add(species_pair)
        observed_species.update(species_pair)
        for first, second in iter_pairs(path, stats):
            output.write(f"{first}\t{second}\n")
            pair_count += 1
    expected_files = len(observed_species) * (len(observed_species) - 1) // 2
    if len(files) != expected_files:
        raise ValueError(
            f"Incomplete SonicParanoid species-pair matrix: found {len(files)} "
            f"tables for {len(observed_species)} species, expected {expected_files}"
        )
    return len(files), pair_count, stats["duplicates"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    files, pairs, duplicates = write_pairs(
        find_pair_directory(args.output_directory), sys.stdout
    )
    print(
        f"Converted {pairs:,} distinct pairs from {files:,} species-pair tables; "
        f"removed {duplicates:,} duplicate relations",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
