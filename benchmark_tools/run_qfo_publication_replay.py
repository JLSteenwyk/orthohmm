"""Run the frozen QfO cached replay without reading accuracy reference labels."""

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.verify_qfo_replay_launcher import verify

CHECKPOINT_SHA = "b90c787f050a087adeeb81b1cace18c9cbb30e1724926a40a1df6d7e49bd549c"


def command_for(root, launcher, output):
    checkpoint = root / "qfo_benchmark/results/orthohmm_high_sensitivity_isolated/output/orthohmm_working_res/high_sensitivity_checkpoint"
    return [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py"),
            "--accuracy-checkpoint", str(checkpoint), "--checkpoint-sha256", CHECKPOINT_SHA,
            "--fasta-directory", str(root / "qfo_benchmark/input"),
            "--output-directory", str(output / "replay"), "--json", str(output / "replay.json"),
            "--cpu", "32", "--matrix", "BLOSUM62", "--cpm-resolution", "0.1",
            "--leiden-seed", "4", "--profile-iterations", "1", "--profile-min-species", "1"]


def compare_partition(expected, observed, universe):
    a = {frozenset(g) for g in read_partition(expected, universe)}
    b = {frozenset(g) for g in read_partition(observed, universe)}
    return {"partition_equal": a == b, "expected_groups": len(a), "observed_groups": len(b),
            "expected_only_groups": len(a - b), "observed_only_groups": len(b - a),
            "expected": file_provenance(expected), "observed": file_provenance(observed),
            "byte_equal": file_provenance(expected)["sha256"] == file_provenance(observed)["sha256"]}


def write_json(path, data):
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    executor = Path(__file__).resolve().parent.parent
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "preflight_failed", "accuracy_evaluated": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "started_at": datetime.now(timezone.utc).isoformat()}
    try:
        before = verify(frozen, launcher, runtime)
        audit_script = executor / "benchmark_tools/audit_qfo_replay_inputs.py"
        def audit(name):
            path = output / name
            subprocess.run([sys.executable, str(audit_script), "--root", str(root), "--output", str(path)], check=True)
            return json.loads(path.read_text())
        inputs = audit("inputs_before.json")
        fastas = {Path(item["path"]) for item in inputs["input_fastas"]}
        actual = {p.resolve() for p in (root / "qfo_benchmark/input").iterdir()
                  if p.suffix.lower() in {".fa", ".faa", ".fasta", ".fsa"}}
        if fastas != actual:
            raise ValueError("Unexpected FASTA files visible to inference")
        command = command_for(root, launcher, output)
        overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                         OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        environment = os.environ.copy()
        environment.update(overrides)
        packages = sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
        preflight = {"launcher": before, "command": command, "environment_overrides": overrides,
                     "executor": file_provenance(Path(__file__)), "python": file_provenance(Path(sys.executable)),
                     "executor_commit": subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip(),
                     "packages": packages, "inputs": file_provenance(output / "inputs_before.json"),
                     "runtime_kind": "incremental_cached_replay_shared_machine",
                     "resource_limitations": "GNU time max RSS is not sampled simultaneous process-tree RSS; shared search cost excluded.",
                     "accuracy_evaluated": False, "job_id": result["job_id"]}
        write_json(output / "preflight.json", preflight)
        result["status"] = "inference_failed"
        with (output / "replay.log").open("x") as log:
            run = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.log"), *command],
                                 cwd=launcher, env=environment, stdout=log, stderr=subprocess.STDOUT)
        result["exit_code"] = run.returncode
        result["postflight_launcher"] = verify(frozen, launcher, runtime)
        after = audit("inputs_after.json")
        if before != result["postflight_launcher"] or inputs != after:
            raise ValueError("Replay sources, runtime or inputs changed during execution")
        if packages != sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions()):
            raise ValueError("Installed package versions changed")
        if run.returncode == 0:
            result["status"] = "verification_failed"
            checkpoint = Path(command[command.index("--accuracy-checkpoint") + 1])
            universe = set((checkpoint / "gene_names.txt").read_text().splitlines())
            result["partition"] = compare_partition(Path(inputs["target_partition"]["path"]),
                                                     output / "replay/orthogroups_profiles_refined.txt", universe)
            result["status"] = "equivalent" if result["partition"]["partition_equal"] else "not_equivalent"
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
        if result.get("exit_code") == 0:
            result["status"] = "verification_failed"
        raise
    finally:
        result["finished_at"] = datetime.now(timezone.utc).isoformat()
        write_json(output / "verification.json", result)
    return 0 if result["status"] == "equivalent" else 1


if __name__ == "__main__":
    raise SystemExit(main())
