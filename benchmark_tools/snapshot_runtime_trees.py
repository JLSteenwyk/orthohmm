"""Inventory explicit runtime trees for prospective before/after identity checks."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def inventory(roots):
    roots = [Path(os.path.abspath(path)) for path in roots]
    if not roots or len(set(roots)) != len(roots):
        raise ValueError("Require unique nonempty runtime roots")
    records = []
    external = []

    def visit(path):
        before = path.lstat()
        mode = stat.S_IMODE(before.st_mode)
        row = {"path": str(path), "mode": mode}
        if path.is_symlink():
            target = path.resolve()
            row.update(kind="symlink", link=os.readlink(path), resolved=str(target),
                       target_exists=target.exists())
            covered = any(target == root or target.is_relative_to(root) for root in roots)
            if not covered:
                external.append(str(path))
            if target.is_file():
                row.update(target_bytes=target.stat().st_size, target_sha256=digest(target))
        elif stat.S_ISREG(before.st_mode):
            row.update(kind="file", bytes=before.st_size, sha256=digest(path))
        elif stat.S_ISDIR(before.st_mode):
            row.update(kind="directory")
        else:
            raise ValueError("Nonregular runtime entry: " + str(path))
        after = path.lstat()
        if (before.st_ino, before.st_size, before.st_mtime_ns, before.st_mode) != (
                after.st_ino, after.st_size, after.st_mtime_ns, after.st_mode):
            raise ValueError("Runtime entry changed while reading: " + str(path))
        records.append(row)
        if row["kind"] == "directory":
            children = sorted(path.iterdir())
            for child in children:
                if child.name not in {"__pycache__", ".git"} and child.suffix not in {".pyc", ".pyo"}:
                    visit(child)
            if children != sorted(path.iterdir()):
                raise ValueError("Runtime directory changed while reading: " + str(path))

    for root in sorted(roots):
        visit(root)
    if len({row["path"] for row in records}) != len(records):
        raise ValueError("Overlapping runtime roots")
    return {"schema": 1, "roots": sorted(map(str, roots)), "records": records,
            "external_symlinks": sorted(external),
            "exclusions": ["__pycache__", ".git", "*.pyc", "*.pyo"],
            "scientific_execution_authorized": False,
            "limitations": [
                "Explicit tree identity, not a hermetic operating-system or loader snapshot.",
                "Directory symlinks are recorded but not traversed; external targets require separate review.",
                "Python bytecode caches are excluded; execution must isolate cache lookup to a fresh prefix and disable writes.",
                "Before/after equality cannot rule out a temporary change during execution."]}


def verify(expected):
    observed = inventory(expected["roots"])
    if observed != expected:
        raise ValueError("Runtime inventory changed")
    return {"status": "runtime_tree_identity_matches", "records": len(observed["records"]),
            "scientific_execution_authorized": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--root", type=Path, action="append")
    group.add_argument("--verify", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = verify(json.loads(args.verify.read_text())) if args.verify else inventory(args.root)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
