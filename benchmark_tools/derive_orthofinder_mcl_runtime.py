#!/usr/bin/env python3
"""Derive elapsed time to OrthoFinder's pre-phylogenetic MCL checkpoint."""

from __future__ import annotations

import argparse
import datetime as dt
import re
from collections.abc import Iterable
from pathlib import Path


TIMESTAMP = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) : (.+)$")


def derive_mcl_runtime(log_paths: Iterable[Path]) -> int:
    """Return seconds from OrthoFinder startup through completed MCL."""
    starts: list[dt.datetime] = []
    mcl_finishes: list[dt.datetime] = []
    for path in log_paths:
        with path.open() as handle:
            for line in handle:
                match = TIMESTAMP.match(line)
                if not match:
                    continue
                when = dt.datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S")
                message = match.group(2)
                if "Started OrthoFinder version 3.1.5" in message:
                    starts.append(when)
                if "Ran MCL" in message:
                    mcl_finishes.append(when)

    if not starts or not mcl_finishes:
        raise ValueError("Could not derive OrthoFinder sequence-only runtime")
    start = min(starts)
    finishes_after_start = [finish for finish in mcl_finishes if finish > start]
    if not finishes_after_start:
        raise ValueError("Could not derive OrthoFinder sequence-only runtime")
    return int((min(finishes_after_start) - start).total_seconds())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()

    runtime = derive_mcl_runtime(args.logs)
    args.output.write_text(f"{runtime}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
