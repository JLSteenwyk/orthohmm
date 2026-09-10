#!/usr/bin/env python3
"""Convert a ProteinOrtho graph to canonical QfO protein pairs."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def _accession(identifier: str) -> str:
    fields = identifier.split("|")
    return fields[1] if len(fields) >= 2 and fields[1] else identifier


def iter_pairs(path: Path) -> Iterator[tuple[str, str]]:
    """Yield validated, distinct pairs from a ProteinOrtho graph."""
    section: tuple[str, str] | None = None
    observed: set[tuple[str, str]] = set()
    observed_sections: set[tuple[str, str]] = set()
    observed_species: set[str] = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                fields = line[1:].strip().split("\t")
                if len(fields) == 2 and all(field.endswith(".fasta") for field in fields):
                    section = (fields[0], fields[1])
                    canonical_section = tuple(sorted(section))
                    if section[0] == section[1]:
                        raise ValueError(
                            f"ProteinOrtho self-species section at {path}:{line_number}"
                        )
                    if canonical_section in observed_sections:
                        raise ValueError(
                            f"Duplicate ProteinOrtho species section at "
                            f"{path}:{line_number}: {canonical_section}"
                        )
                    observed_sections.add(canonical_section)
                    observed_species.update(section)
                    observed.clear()
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(
                    f"Expected six ProteinOrtho graph columns at {path}:{line_number}"
                )
            if section is None:
                raise ValueError(
                    f"ProteinOrtho relation precedes species header at {path}:{line_number}"
                )
            try:
                for value in fields[2:]:
                    float(value)
            except ValueError as error:
                raise ValueError(
                    f"Invalid ProteinOrtho score at {path}:{line_number}"
                ) from error

            pair = tuple(sorted((_accession(fields[0]), _accession(fields[1]))))
            if pair[0] == pair[1]:
                raise ValueError(f"Self-pair in ProteinOrtho graph: {pair[0]}")
            if pair in observed:
                raise ValueError(
                    f"Duplicate ProteinOrtho pair in section {section}: {pair}"
                )
            observed.add(pair)
            yield pair

    expected_sections = len(observed_species) * (len(observed_species) - 1) // 2
    if len(observed_sections) != expected_sections:
        raise ValueError(
            f"Incomplete ProteinOrtho species-pair matrix: found "
            f"{len(observed_sections)} sections for {len(observed_species)} species, "
            f"expected {expected_sections}"
        )


def write_pairs(path: Path, output: TextIO) -> int:
    count = 0
    for first, second in iter_pairs(path):
        output.write(f"{first}\t{second}\n")
        count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path)
    args = parser.parse_args()
    count = write_pairs(args.graph, sys.stdout)
    print(f"Converted {count:,} native ProteinOrtho pairs", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
