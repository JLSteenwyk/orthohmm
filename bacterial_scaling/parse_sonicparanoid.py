#!/usr/bin/env python3
"""Convert SonicParanoid's ortholog_groups.tsv to one-OG-per-line format."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark_tools.normalize_three_kingdoms_orthogroups import (  # noqa: E402
    iter_sonicparanoid,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()

    observed: set[str] = set()
    for group, genes in iter_sonicparanoid(args.input):
        duplicate = observed.intersection(genes)
        if duplicate:
            raise ValueError(
                f"Genes occur in multiple groups; first duplicate in {group}: "
                f"{min(duplicate)}"
            )
        observed.update(genes)
        if len(genes) >= 2:
            print(" ".join(genes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
