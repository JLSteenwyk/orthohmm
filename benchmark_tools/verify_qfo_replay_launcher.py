"""Verify the isolated numeric QfO replay against the frozen publication core."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.build_publication_runtime import COMMIT, verify_runtime
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.validate_profile_runtime import require_profile_runtime

LAUNCHER_COMMIT = "49ab110358c0b4c73806a640de9068494a311f63"
SOURCE_SUFFIXES = {".py", ".c", ".cu", ".h"}


def compare_trees(frozen, launcher):
    def inventory(root):
        return {p.relative_to(root): p for p in (root / "orthohmm").rglob("*")
                if p.is_file() and p.suffix in SOURCE_SUFFIXES | {".so"}}
    expected, observed = inventory(frozen), inventory(launcher)
    if expected.keys() != observed.keys():
        raise ValueError("Frozen and launcher source/native file sets differ")
    evidence = []
    for relative in sorted(expected):
        a, b = file_provenance(expected[relative]), file_provenance(observed[relative])
        if a["sha256"] != b["sha256"] or a["bytes"] != b["bytes"]:
            raise ValueError("Frozen and launcher file differs: " + str(relative))
        evidence.append({"relative_path": str(relative), "frozen": a, "launcher": b})
    return evidence


def verify(frozen, launcher, runtime_manifest):
    frozen, launcher = Path(frozen).resolve(), Path(launcher).resolve()
    for root, revision in ((frozen, COMMIT), (launcher, LAUNCHER_COMMIT)):
        actual = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        if actual != revision:
            raise ValueError("Wrong pinned checkout revision")
        subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--",
                        "orthohmm", "benchmark_tools", "setup.py"], check=True, capture_output=True)
    verify_runtime(runtime_manifest, frozen)
    evidence = compare_trees(frozen, launcher)
    probe = require_profile_runtime(launcher)
    return {"schema_version": 1, "status": "verified", "core_commit": COMMIT,
            "launcher_commit": LAUNCHER_COMMIT, "files": evidence,
            "runtime_manifest": file_provenance(runtime_manifest), "profile_probe": probe,
            "replay_script": file_provenance(launcher / "benchmark_tools/replay_high_sensitivity.py"),
            "accuracy_evaluated": False, "replay_equivalence": "not yet evaluated",
            "limitations": ["The frozen core retains historical exception-to-None behavior; the development fix is not reverted.",
                            "A successful synthetic profile probe does not establish successful construction for every biological cluster."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-root", required=True, type=Path)
    parser.add_argument("--launcher-root", required=True, type=Path)
    parser.add_argument("--runtime-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite launcher provenance")
    report = verify(args.frozen_root, args.launcher_root, args.runtime_manifest)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
