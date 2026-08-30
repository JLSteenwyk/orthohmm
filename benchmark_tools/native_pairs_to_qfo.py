#!/usr/bin/env python3
"""Convert OrthoHMM native reconciled pairs to two-column QfO input."""

from __future__ import annotations

import argparse
from pathlib import Path


def strip_to_uniprot(identifier: str) -> str:
    if "|" not in identifier:
        return identifier
    fields = identifier.split("|")
    return fields[1] if len(fields) >= 2 and fields[1] else identifier


def convert(input_path: Path, output_path: Path) -> int:
    emitted = 0
    with input_path.open() as input_handle, output_path.open("w") as output_handle:
        header = input_handle.readline().rstrip("\n").split("\t")
        if header != ["gene_a", "species_a", "gene_b", "species_b"]:
            raise ValueError(f"unexpected native pair header in {input_path}")
        for line_number, line in enumerate(input_handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(
                    f"{input_path}:{line_number}: expected four tab-separated fields"
                )
            left = strip_to_uniprot(fields[0])
            right = strip_to_uniprot(fields[2])
            if left == right:
                continue
            if left > right:
                left, right = right, left
            output_handle.write(f"{left}\t{right}\n")
            emitted += 1
    return emitted


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    emitted = convert(args.input, args.output)
    print(f"Emitted {emitted} native reconciled QfO pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
