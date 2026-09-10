#!/usr/bin/env python3
"""Restore original FASTA identifiers in OrthoFinder orthogroup output."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import TextIO


FASTA_SUFFIXES = (".fa", ".faa", ".fasta", ".fas")


def orthofinder_identifier(identifier: str) -> str:
    """Apply the accession sanitization used by OrthoFinder 3.1.5."""
    return (
        identifier.replace(":", "_")
        .replace(",", "_")
        .replace("(", "_")
        .replace(")", "_")
    )


def iter_fasta_identifiers(paths: Iterable[Path]) -> Iterable[str]:
    for path in paths:
        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.startswith(">"):
                    continue
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA header: {path}:{line_number}")
                yield fields[0]


def build_restoration_map(input_dir: Path) -> dict[str, str]:
    fasta_paths = sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )
    if not fasta_paths:
        raise ValueError(f"No FASTA files found in {input_dir}")

    restoration: dict[str, str] = {}
    for original in iter_fasta_identifiers(fasta_paths):
        sanitized = orthofinder_identifier(original)
        previous = restoration.setdefault(sanitized, original)
        if previous != original:
            raise ValueError(
                "OrthoFinder identifier sanitization is ambiguous: "
                f"{previous!r} and {original!r} both become {sanitized!r}"
            )
    return restoration


def restore_orthogroups(
    source: TextIO, destination: TextIO, restoration: Mapping[str, str]
) -> tuple[int, int, int]:
    """Write orthogroups with restored IDs and return group/gene/change counts."""
    groups = genes = changed = 0
    seen: set[str] = set()

    for line_number, line in enumerate(source, start=1):
        fields = line.split()
        if not fields:
            continue
        label = fields[0] if fields[0].endswith(":") else None
        members = fields[1:] if label else fields
        if not members:
            raise ValueError(f"Orthogroup has no members on line {line_number}")

        restored_members = []
        for member in members:
            if member not in restoration:
                raise ValueError(
                    f"Unknown OrthoFinder identifier on line {line_number}: {member}"
                )
            original = restoration[member]
            if original in seen:
                raise ValueError(
                    f"Sequence identifier appears in multiple orthogroups: {original}"
                )
            seen.add(original)
            restored_members.append(original)
            genes += 1
            changed += original != member

        output_fields = ([label] if label else []) + restored_members
        destination.write(" ".join(output_fields) + "\n")
        groups += 1

    return groups, genes, changed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("orthogroups", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    restoration = build_restoration_map(args.input_dir)
    with args.orthogroups.open() as source, args.output.open("w") as destination:
        groups, genes, changed = restore_orthogroups(source, destination, restoration)
    print(
        f"Restored {changed:,} of {genes:,} identifiers across {groups:,} "
        "OrthoFinder orthogroups"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
