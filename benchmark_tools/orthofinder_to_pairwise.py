#!/usr/bin/env python3
"""Convert OrthoFinder's native orthologue tables to QfO gene pairs."""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import product
from pathlib import Path
from typing import Iterable, TextIO


def _genes(cell: str) -> list[str]:
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def _accession(gene: str) -> str:
    parts = gene.split("|")
    return parts[1] if len(parts) >= 2 and parts[1] else gene


def iter_pairs(results_dir: Path) -> Iterable[tuple[str, str]]:
    """Yield each native OrthoFinder pair once across unordered species pairs."""
    orthologues_dir = results_dir / "Orthologues"
    if not orthologues_dir.is_dir():
        raise FileNotFoundError(f"OrthoFinder Orthologues directory not found: {orthologues_dir}")

    pair_files = sorted(orthologues_dir.glob("Orthologues_*/*__v__*.tsv"))
    if not pair_files:
        raise FileNotFoundError(f"No pairwise OrthoFinder tables found in {orthologues_dir}")

    for path in pair_files:
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None or len(header) != 3 or header[0] != "Orthogroup":
                raise ValueError(f"Unexpected OrthoFinder table header in {path}: {header}")

            species_a, species_b = header[1], header[2]
            if species_a >= species_b:
                continue

            for row in reader:
                if len(row) != 3:
                    raise ValueError(f"Expected three columns in {path}: {row}")
                for gene_a, gene_b in product(_genes(row[1]), _genes(row[2])):
                    acc_a, acc_b = _accession(gene_a), _accession(gene_b)
                    yield (acc_a, acc_b) if acc_a < acc_b else (acc_b, acc_a)


def write_pairs(results_dir: Path, output: TextIO) -> int:
    count = 0
    for gene_a, gene_b in iter_pairs(results_dir):
        output.write(f"{gene_a}\t{gene_b}\n")
        count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path, help="OrthoFinder Results_* directory")
    args = parser.parse_args()
    count = write_pairs(args.results_dir, sys.stdout)
    print(f"Emitted {count:,} native OrthoFinder pairs", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
