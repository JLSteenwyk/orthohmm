#!/usr/bin/env python3
"""Freeze a label-independent OrthoBench development/validation partition."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
)


PARTITION_SALT = "orthohmm-phylogeny-candidates-v1"


def partition_names(names: Sequence[str]) -> tuple[list[str], list[str]]:
    """Split names deterministically without reading benchmark labels."""

    unique = sorted(set(names))
    ranked = sorted(
        unique,
        key=lambda name: (
            hashlib.sha256(f"{PARTITION_SALT}:{name}".encode()).hexdigest(),
            name,
        ),
    )
    midpoint = (len(ranked) + 1) // 2
    return sorted(ranked[:midpoint]), sorted(ranked[midpoint:])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refogs", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    refog_directory = args.refogs.resolve()
    paths = sorted(refog_directory.glob("RefOG*.txt"))
    if len(paths) < 2:
        raise SystemExit(f"fewer than two RefOG files found in {refog_directory}")
    development, validation = partition_names([path.name for path in paths])
    records = {path.name: file_provenance(path) for path in paths}
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": (
            "Predefined OrthoBench development/validation partition; assignment "
            "uses filenames only and does not inspect RefOG membership labels."
        ),
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "partition_salt": PARTITION_SALT,
        "development": development,
        "validation": validation,
        "refogs": records,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.json)
    print(
        f"development={len(development)} validation={len(validation)} "
        f"output={args.json.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
