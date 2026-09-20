"""Verify the source-only scientific baseline export against immutable Git blobs."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

REVISION = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
PREFIX = "orthohmm-0.5.0-scientific-7f3a9e4/"
PATHS = ["LICENSE.md", "README.md", "requirements.txt", "setup.py", "orthohmm"]


def git_files(repo):
    raw = subprocess.check_output(["git", "-C", str(repo), "ls-tree", "-rz", REVISION, "--", *PATHS])
    files = {}
    for entry in raw.split(b"\0"):
        if not entry:
            continue
        header, name = entry.split(b"\t", 1)
        mode, kind, oid = header.decode().split()
        if kind != "blob" or mode not in ("100644", "100755"):
            raise ValueError("Expected only regular source blobs")
        name = name.decode()
        if name in files:
            raise ValueError("Duplicate source path")
        files[name] = dict(git_blob=oid, mode=int(mode, 8) & 0o777,
            content=subprocess.check_output(["git", "-C", str(repo), "cat-file", "blob", oid]))
    if not files:
        raise ValueError("Empty source inventory")
    return files


def verify(repo, archive_path):
    identity = record(archive_path)
    expected = git_files(repo)
    directories = {str(parent) for name in expected for parent in PurePosixPath(PREFIX + name).parents
                   if str(parent) != "."}
    seen, rows, syntax = set(), [], 0
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            if member.name in seen:
                raise ValueError("Duplicate archive member")
            seen.add(member.name)
            if member.isdir():
                if member.name.rstrip("/") not in directories:
                    raise ValueError("Unexpected archive directory")
                continue
            if not member.isfile() or not member.name.startswith(PREFIX):
                raise ValueError("Require regular source files under the fixed archive prefix")
            name = member.name[len(PREFIX):]
            if name not in expected or member.mode != expected[name]["mode"]:
                raise ValueError("Unexpected archive file or mode")
            content = archive.extractfile(member).read()
            if content != expected[name]["content"]:
                raise ValueError("Archive source differs from frozen Git blob")
            if name.endswith(".py"):
                compile(content, name, "exec")
                syntax += 1
            rows.append(dict(path=name, mode=member.mode, bytes=len(content),
                sha256=hashlib.sha256(content).hexdigest(), git_blob=expected[name]["git_blob"]))
    if {row["path"] for row in rows} != set(expected):
        raise ValueError("Source inventory incomplete")
    check(identity)
    return dict(status="frozen_source_archive_matches_git", scientific_revision=REVISION,
        prefix=PREFIX, selected_paths=PATHS, archive=identity, files=sorted(rows, key=lambda r:r["path"]),
        python_files_syntax_checked=syntax, python=sys.version, source=record(__file__),
        executable_benchmark_reproduced=False, redistribution_clearance=False, publication_ready=False,
        limitations=["Source-byte and mode verification only; syntax compilation does not execute or import code.",
            "Not a full repository, dependency/runtime archive, installed build or numerical-equivalence test.",
            "Historical setup/build behavior and dependency declarations remain unchanged.",
            "No public release, external deposition, complete rights review or publication readiness asserted."])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = verify(args.repo, args.archive)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
