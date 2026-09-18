"""Export the exact historical method-figure helper without executing it."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
SOURCE = "benchmark_tools/replay_high_sensitivity.py"
HISTORICAL_PATH = "benchmarks/work/publication_method_native_v2/" + SOURCE
SHA256 = "852c3e4fc1a53de7e6046aa78324da587376c0df6a0db1cd8265af65c75bea0f"
SIZE = 14723
BLOB = "55dc223310a7ba764bff0955c6ef7f23185c1951"


def export(repo, output):
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    spec = f"{COMMIT}:{SOURCE}"
    blob = subprocess.check_output(
        ["git", "rev-parse", "--verify", spec], cwd=repo, text=True).strip()
    payload = subprocess.check_output(["git", "show", spec], cwd=repo)
    if blob != BLOB or len(payload) != SIZE or hashlib.sha256(payload).hexdigest() != SHA256:
        raise ValueError("Frozen helper identity mismatch")
    manifest = {
        "schema_version": 1,
        "status": "frozen_figure_helper_exported",
        "publication_ready": False,
        "source_commit": COMMIT,
        "source_git_path": SOURCE,
        "source_git_blob": blob,
        "historical_repository_relative_path": HISTORICAL_PATH,
        "export": {"path": SOURCE, "bytes": SIZE, "sha256": SHA256},
        "limitations": [
            "Source recovery only; the exported helper is not executed.",
            "Historical figure manifests are unchanged; relocation is explicit in this manifest.",
            "Imported modules, runtimes and raw data are not bundled or validated.",
            "This is not the complete publication archive or a runnable method distribution.",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    target = output / SOURCE
    target.parent.mkdir(parents=True)
    with target.open("xb") as stream:
        stream.write(payload)
    recovered = target.read_bytes()
    if recovered != payload:
        raise ValueError("Export read-back mismatch")
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.repo, args.output), indent=2, sort_keys=True))
