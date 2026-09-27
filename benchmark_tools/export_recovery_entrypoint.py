"""Export the committed recovery entrypoint without copying the full checkout."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


REVISION = "ebac5b26a06fa389192bbd54a01737ceb3f7a425"
FILES = (
    "benchmark_tools/run_publication_pipeline.py",
    "benchmark_tools/candidate_hit_order_policy.py",
    "benchmark_tools/probe_ob_canonical_candidates.py",
    "benchmark_tools/results/publication_recovery_requirements_20260926.txt",
)


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def verify(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest["revision"] != REVISION or tuple(manifest["files"]) != FILES:
        raise ValueError("Unexpected export revision or inventory")
    for name, expected in manifest["files"].items():
        path = directory / name
        if path.is_symlink() or identity(path.read_bytes()) != expected:
            raise ValueError("Changed export: " + name)
    return manifest


def export(repo, directory):
    directory = Path(directory).absolute()
    if directory.exists() or directory.is_symlink():
        raise FileExistsError(directory)
    payload = {name: subprocess.check_output(
        ["git", "show", REVISION + ":" + name], cwd=repo) for name in FILES}
    directory.mkdir(parents=True)
    manifest = dict(revision=REVISION, files={}, publication_ready=False,
        scope="Relocatable runtime harness only; no native execution implied",
        prerequisites=["Validated installed recovery Python environment",
            "MAFFT with helper executables and FastTree",
            "Input FASTA files and a fresh output directory"],
        limitations=["Wheels, native tools and datasets are not included",
            "The diagnostic helper's other entrypoints are not exported workflows",
            "Hashes detect changes relative to this manifest, not manifest authenticity"])
    for name, data in payload.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        manifest["files"][name] = identity(data)
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    verify(directory)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.repo, args.output), indent=2))
