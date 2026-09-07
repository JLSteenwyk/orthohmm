#!/usr/bin/env python3
"""Configure an isolated OrthoMCL 1.4 module for a benchmark run."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def replace_assignment(source: str, name: str, value: str) -> str:
    """Replace exactly one top-level Perl variable assignment."""
    pattern = re.compile(rf"^our \${re.escape(name)}\s*=.*$", re.MULTILINE)
    updated, count = pattern.subn(lambda _: f"our ${name} = {value};", source)
    if count != 1:
        raise ValueError(f"Expected one ${name} assignment, found {count}")
    return updated


def configure_module(path: Path, tool_dir: Path, data_dir: Path, threads: int) -> None:
    """Write benchmark-local paths and CPU count to an OrthoMCL module."""
    if threads < 1:
        raise ValueError("threads must be positive")
    source = path.read_text()
    source = replace_assignment(source, "PATH_TO_ORTHOMCL", f'"{tool_dir}/"')
    source = replace_assignment(source, "ORTHOMCL_DATA_DIR", f'"{data_dir}/"')
    source = replace_assignment(source, "BLAST_NOCPU", str(threads))
    path.write_text(source)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("module", type=Path)
    parser.add_argument("tool_dir", type=Path)
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("threads", type=int)
    args = parser.parse_args()
    configure_module(args.module, args.tool_dir, args.data_dir, args.threads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
