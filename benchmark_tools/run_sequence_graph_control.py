"""Replay frozen profile-off graph inference from validated sequence-search hits."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import read_partition, verify_file
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.verify_qfo_replay_launcher import verify as verify_launcher
from benchmark_tools.verify_ygob_validation import require_completed_job


def graph_command(launcher, checkpoint, sha256, output):
    # Omitting FASTAs is the frozen replay's existing profile-off path.
    return [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py"),
            "--accuracy-checkpoint", str(checkpoint), "--checkpoint-sha256", sha256,
            "--output-directory", str(output / "replay"), "--json", str(output / "replay.json"),
            "--cpu", "32", "--matrix", "BLOSUM62", "--cpm-resolution", "0.1", "--leiden-seed", "4"]


def check_profile_off(metrics):
    if metrics["parameters"]["profile_expansion"] is not False:
        raise ValueError("Sequence-only control unexpectedly used profile expansion")
    if [row["label"] for row in metrics["stages"]] != ["multipass", "multipass_refined"]:
        raise ValueError("Unexpected graph-only stages")
    if metrics["counts"]["genes"] != 251378 or metrics["counts"]["species"] != 12:
        raise ValueError("Graph control input universe differs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--variant", required=True, choices=("all_hits", "top100"))
    args = parser.parse_args()
    root = args.root.resolve()
    conversion = root / "benchmarks/results/ob_sequence_numeric_v1/manifest.json"
    accounting = subprocess.check_output(["sacct", "-j", "21292", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21292)
    converted = json.loads(conversion.read_text())
    if converted["status"] != "numeric_checkpoints_verified" or converted["accuracy_evaluated"] is not False:
        raise ValueError("Sequence conversion is not admitted")
    if converted["scheduler"]["JobIDRaw"] != "21291" or set(converted["variants"]) != {"all_hits", "top100"}:
        raise ValueError("Unexpected conversion provenance or variants")
    for item in (converted["source"], converted["execution"]):
        verify_file(Path(item["path"]), item)
    variant = converted["variants"][args.variant]
    if variant["cap"] != (None if args.variant == "all_hits" else 100):
        raise ValueError("Wrong hit-cap configuration")
    checkpoint = Path(variant["checkpoint"])
    verify_file(checkpoint / "manifest.json", variant["manifest"])
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    before = verify_launcher(frozen, launcher, runtime)
    output = root / "benchmarks/results/ob_sequence_graph_v1" / args.variant
    output.mkdir(parents=True, exist_ok=False)
    command = graph_command(launcher, checkpoint, variant["manifest"]["sha256"], output)
    overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                     OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env = os.environ.copy()
    env.update(overrides)
    provenance = {"variant": args.variant, "command": command, "conversion": file_provenance(conversion),
                  "conversion_scheduler": scheduler, "runtime": before, "environment_overrides": overrides,
                  "job_id": os.environ.get("SLURM_JOB_ID"), "source": file_provenance(Path(__file__)),
                  "accuracy_evaluated": False, "runtime_kind": "incremental graph-only shared-machine replay"}
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    result = {"status": "inference_failed", "accuracy_evaluated": False}
    try:
        with (output / "replay.log").open("x") as handle:
            process = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.log"), *command],
                                     cwd=launcher, env=env, stdout=handle, stderr=subprocess.STDOUT)
        result["exit_code"] = process.returncode
        process.check_returncode()
        result["status"] = "verification_failed"
        metrics = json.loads((output / "replay.json").read_text())
        check_profile_off(metrics)
        universe = set((checkpoint / "gene_names.txt").read_text().splitlines())
        prediction = output / "replay/orthogroups_multipass_refined.txt"
        groups = read_partition(prediction, universe)
        verify_file(conversion, provenance["conversion"])
        verify_file(checkpoint / "manifest.json", variant["manifest"])
        audit(checkpoint, variant["manifest"]["sha256"])
        if verify_launcher(frozen, launcher, runtime) != before:
            raise ValueError("Frozen runtime changed")
        result.update(status="graph_complete_pending_scoring", groups=len(groups), genes=len(universe),
                      prediction=file_provenance(prediction), metrics=file_provenance(output / "replay.json"))
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
