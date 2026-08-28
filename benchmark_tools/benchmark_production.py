#!/usr/bin/env python3
"""Run and record a reproducible benchmark of the production OrthoHMM CLI."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


FASTA_EXTENSIONS = (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, base: Path) -> dict:
    return {
        "path": str(path.relative_to(base)),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def count_nonempty_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def git_value(repo_root: Path, *args: str) -> str | None:
    completed = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_directory", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("result_json", type=Path)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--threads-per-worker", type=int, default=8)
    parser.add_argument("--matrix", default="BLOSUM62")
    parser.add_argument("--evalue", type=float, default=1e-4)
    parser.add_argument("--clustering", choices=("leiden", "mcl"), default="leiden")
    parser.add_argument("--cpm-resolution", default="0.1")
    parser.add_argument("--refinement-profile", choices=("default", "qfo"), default="default")
    parser.add_argument(
        "--accuracy-profile",
        choices=("standard", "high_sensitivity"),
        default="standard",
    )
    parser.add_argument("--full-output", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parent.parent
    input_directory = args.input_directory.resolve()
    output_directory = args.output_directory.resolve()
    result_json = args.result_json.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    result_json.parent.mkdir(parents=True, exist_ok=True)

    existing_outputs = list(output_directory.glob("orthohmm_*"))
    if existing_outputs:
        raise SystemExit(
            f"output directory already contains OrthoHMM outputs: {output_directory}"
        )

    command = [
        sys.executable,
        "-m",
        "orthohmm",
        str(input_directory),
        "-o",
        str(output_directory),
        "-c",
        str(args.cpu),
        "--threads_per_worker",
        str(args.threads_per_worker),
        "-x",
        args.matrix,
        "-e",
        str(args.evalue),
        "--clustering",
        args.clustering,
        "--cpm_resolution",
        args.cpm_resolution,
        "--refinement_profile",
        args.refinement_profile,
        "--accuracy_profile",
        args.accuracy_profile,
        "--metrics_json",
        str(result_json),
    ]
    if not args.full_output:
        command.extend(["--stop", "infer"])

    log_path = output_directory / "benchmark.log"
    environment = os.environ.copy()
    environment["PYTHONUNBUFFERED"] = "1"
    with log_path.open("w") as log:
        completed = subprocess.run(
            command,
            cwd=repo_root,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if result_json.exists():
        result = json.loads(result_json.read_text())
    else:
        result = {"schema_version": 1, "status": "failed"}

    fasta_files = sorted(
        path for path in input_directory.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )
    output_candidates = [
        output_directory / "orthohmm_orthogroups.txt",
        output_directory / "orthohmm_working_res" / "orthohmm_edges.txt",
        output_directory / "orthohmm_working_res" / "orthohmm_edges_clustered.txt",
        log_path,
    ]
    source_files = sorted((repo_root / "orthohmm").rglob("*.py"))
    source_files.append(Path(__file__).resolve())
    git_status = git_value(repo_root, "status", "--porcelain=v1")
    result["harness"] = {
        "name": "benchmark_tools.benchmark_production",
        "schema_version": 1,
        "exit_code": completed.returncode,
        "command": command,
        "git_commit": git_value(repo_root, "rev-parse", "HEAD"),
        "git_dirty": bool(git_status),
        "source_manifest": [
            file_record(path, repo_root) for path in source_files
        ],
        "input_manifest": [
            file_record(path, input_directory) for path in fasta_files
        ],
        "output_manifest": [
            file_record(path, output_directory)
            for path in output_candidates
            if path.exists()
        ],
        "output_line_counts": {
            "orthogroups": count_nonempty_lines(output_candidates[0]),
            "edges": count_nonempty_lines(output_candidates[1]),
            "clusters_before_singletons": count_nonempty_lines(output_candidates[2]),
        },
    }
    result_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result_json)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
