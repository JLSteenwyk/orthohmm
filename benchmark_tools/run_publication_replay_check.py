"""Run a pinned, label-blind OrthoBench replay and verify stage equivalence."""

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import STAGES, read_partition, verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.validate_profile_runtime import require_profile_runtime
from benchmark_tools.build_publication_runtime import verify_runtime


OUTPUTS = dict(zip(STAGES, ("orthogroups_multipass.txt", "orthogroups_multipass_refined.txt",
                           "orthogroups_profiles.txt", "orthogroups_profiles_refined.txt")))


def compare_stages(expected, output, universe):
    result = {}
    for name, filename in OUTPUTS.items():
        source = Path(expected[name]["path"])
        verify_file(source, expected[name])
        target = output / filename
        reference = {frozenset(g) for g in read_partition(source, universe)}
        observed = {frozenset(g) for g in read_partition(target, universe)}
        actual = file_provenance(target)
        result[name] = {"partition_equal": reference == observed,
                        "byte_equal": actual["sha256"] == expected[name]["sha256"],
                        "expected_groups": len(reference), "observed_groups": len(observed),
                        "expected_only_groups": len(reference - observed),
                        "observed_only_groups": len(observed - reference), "output": actual}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--historical-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=32)
    parser.add_argument("--runtime-manifest", type=Path, required=True)
    args = parser.parse_args()
    root = args.frozen_root.resolve()
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if not commit.startswith("7f3a9e4") or args.cpu < 1:
        raise ValueError("Unexpected source revision or CPU budget")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--", "orthohmm",
                    "benchmark_tools/replay_high_sensitivity.py"], check=True, stdout=subprocess.DEVNULL)
    audit = json.loads(args.historical_audit.read_text())
    if not audit["fresh_partition_byte_identical"] or audit["genes"] != 251378:
        raise ValueError("Historical replay gate not satisfied")
    if args.output.exists():
        raise ValueError("Refusing to overwrite replay check")
    verify_runtime(args.runtime_manifest, root)
    runtime_manifest = file_provenance(args.runtime_manifest)
    profile_runtime = require_profile_runtime(root)
    cache = audit["inputs"]["cache"]
    verify_file(Path(cache["path"]), cache)
    universe = set()
    fasta_parents = set()
    for item in audit["inputs"]["fastas"]:
        path = Path(item["path"])
        verify_file(path, item)
        fasta_parents.add(path.parent)
        for record in SeqIO.parse(path, "fasta"):
            if record.id in universe:
                raise ValueError("Duplicate input ID")
            universe.add(record.id)
    if len(fasta_parents) != 1 or len(universe) != audit["genes"]:
        raise ValueError("Unexpected input universe")
    fasta_directory = next(iter(fasta_parents))
    observed = {p.resolve() for p in fasta_directory.iterdir() if p.suffix.lower() in {".fa", ".faa", ".fasta", ".fsa"}}
    if observed != {Path(r["path"]).resolve() for r in audit["inputs"]["fastas"]}:
        raise ValueError("Unmanifested FASTA input")
    for item in audit["inputs"]["stage_partitions"].values():
        verify_file(Path(item["path"]), item)
    output = args.output.resolve()
    output.mkdir(parents=True)
    command = [sys.executable, str(root / "benchmark_tools/replay_high_sensitivity.py"),
               "--hits-pickle", cache["path"], "--fasta-directory", str(fasta_directory),
               "--output-directory", str(output / "replay"), "--json", str(output / "replay.json"),
               "--cpu", str(args.cpu), "--matrix", "BLOSUM62", "--cpm-resolution", "0.1",
               "--leiden-seed", "4", "--profile-iterations", "1", "--profile-min-species", "1"]
    environment = os.environ.copy()
    environment.update(PYTHONPATH=str(root), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    metadata = {"source_commit": commit, "launcher": file_provenance(Path(__file__)),
                "profile_runtime": profile_runtime,
                "runtime_manifest": runtime_manifest,
                "historical_audit": file_provenance(args.historical_audit), "command": command,
                "source_manifest": [file_provenance(p) for p in sorted((root / "orthohmm").rglob("*"))
                                    if p.is_file() and p.suffix in {".py", ".c", ".cu", ".h"}],
                "replay_source": file_provenance(root / "benchmark_tools/replay_high_sensitivity.py"),
                "environment": {k: environment[k] for k in ("PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
                "python": sys.version, "packages": sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions()),
                "job_id": os.environ.get("SLURM_JOB_ID"), "started_at": datetime.now(timezone.utc).isoformat(),
                "accuracy_scoring_requested": False, "runtime_kind": "incremental_cached_replay"}
    (output / "preflight.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    with (output / "replay.log").open("w") as log:
        completed = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.log"), *command],
                                   cwd=root, env=environment, stdout=log, stderr=subprocess.STDOUT)
    result = {"exit_code": completed.returncode, "finished_at": datetime.now(timezone.utc).isoformat(),
              "preflight": file_provenance(output / "preflight.json"), "status": "inference_failed"}
    if completed.returncode == 0:
        try:
            verify_file(args.runtime_manifest, runtime_manifest)
            verify_runtime(args.runtime_manifest, root)
            result["stages"] = compare_stages(audit["inputs"]["stage_partitions"], output / "replay", universe)
            result["status"] = "equivalent" if all(v["partition_equal"] for v in result["stages"].values()) else "not_equivalent"
        except (ValueError, OSError) as error:
            result.update(status="verification_failed", error=str(error))
    (output / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"])
    return 0 if result["status"] == "equivalent" else 1


if __name__ == "__main__":
    raise SystemExit(main())
