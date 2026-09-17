"""Rebuild HMM grouping for a fixed CPM neighborhood, with an unchanged control."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import PREPARED_SHA, check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.build_publication_runtime import verify_runtime
from benchmark_tools.validate_profile_runtime import require_profile_runtime
from benchmark_tools.run_publication_replay_check import OUTPUTS, compare_stages
from benchmark_tools.audit_historical_profile_ablation import read_partition

ARMS = (("control", "0.1"), ("cpm_low", "0.08"), ("cpm_high", "0.12"))
RUNTIME_SHA = "aebea83807356b02307473506fa30c2dbd2c511d7ba75a0655eb12180a474d74"


def replay_command(original, label, output):
    resolutions = dict(ARMS)
    if label not in resolutions:
        raise ValueError("Unknown CPM arm")
    argv = list(original)
    required = {"--cpm-resolution": "0.1", "--leiden-seed": "4", "--cpu": "32",
                "--matrix": "BLOSUM62", "--profile-iterations": "1", "--profile-min-species": "1"}
    for flag, value in required.items():
        if argv.count(flag) != 1 or argv[argv.index(flag) + 1] != value:
            raise ValueError("Changed frozen replay parameters")
    if "--official-benchmark" in argv:
        raise ValueError("Accuracy scoring is not allowed")
    for flag, value in (("--cpm-resolution", resolutions[label]),
                        ("--output-directory", str(output / "replay")),
                        ("--json", str(output / "replay.json"))):
        if argv.count(flag) != 1:
            raise ValueError("Missing or duplicated replay destination")
        argv[argv.index(flag) + 1] = value
    return argv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    prepared_path = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    check(prepared["replay_verification"])
    verification = json.loads(Path(prepared["replay_verification"]["path"]).read_text())
    if verification["status"] != "equivalent":
        raise ValueError("Baseline replay not verified")
    check(verification["preflight"])
    preflight = json.loads(Path(verification["preflight"]["path"]).read_text())
    frozen = root / "benchmarks/work/publication_method_native_v2"
    replay = frozen / "benchmark_tools/replay_high_sensitivity.py"
    if preflight["command"][:2] != [sys.executable, str(replay)]:
        raise ValueError("Unexpected replay executable")
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    read_frozen(runtime_path, RUNTIME_SHA)
    verify_runtime(runtime_path, frozen)
    native = require_profile_runtime(frozen)
    inputs = [record(prepared_path), record(runtime_path), prepared["replay_verification"],
        verification["preflight"], preflight["replay_source"], *prepared["core_sources"],
        *prepared["fasta_inputs"], prepared["cache"]]
    expected = {name: verification["stages"][name]["output"] for name in OUTPUTS}
    inputs.extend(expected.values())
    for item in inputs:
        check(item)
    universe = set()
    for item in prepared["fasta_inputs"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in universe:
                raise ValueError("Duplicate FASTA identifier")
            universe.add(sequence.id)
    env = os.environ.copy()
    env.update(prepared["environment_overrides"], PYTHONPATH=str(frozen))
    output.mkdir(parents=True)
    report = {"status": "running_unscored", "accuracy_evaluated": False,
        "job_id": os.environ.get("SLURM_JOB_ID"), "source": record(__file__),
        "inputs": inputs, "profile_runtime": native, "arms": [],
        "environment_overrides": {**prepared["environment_overrides"], "PYTHONPATH": str(frozen)},
        "limitations": ["CPM resolution alone varies; graph clustering, singleton assignment, profiles and refinement rebuilt from fixed cached hits.",
            "Candidate expansion and inferred phylogeny remain separate downstream requirements.",
            "Unchanged control must reproduce all four baseline stage partitions before variants execute.",
            "Incremental shared-node measurements exclude original all-to-all search; no accuracy-based selection."]}
    try:
        for label, resolution in ARMS:
            directory = output / label
            directory.mkdir()
            command = replay_command(preflight["command"], label, directory)
            row = {"label": label, "resolution": float(resolution), "command": command, "status": "running"}
            report["arms"].append(row)
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            started = time.monotonic()
            with (directory / "replay.log").open("x") as log:
                run = subprocess.run(["/usr/bin/time", "-v", "-o", str(directory / "time.log"), *command],
                    cwd=frozen, env=env, stdout=log, stderr=subprocess.STDOUT)
            row.update(exit_code=run.returncode, wall_s=time.monotonic() - started)
            if run.returncode != 0:
                raise RuntimeError("CPM replay failed; output preserved: " + label)
            row["stages"] = {}
            for name, filename in OUTPUTS.items():
                path = directory / "replay" / filename
                groups = read_partition(path, universe)
                row["stages"][name] = {"output": record(path), "groups": len(groups)}
            row["metrics"] = record(directory / "replay.json")
            if label == "control":
                row["control_equivalence"] = compare_stages(expected, directory / "replay", universe)
                if not all(item["byte_equal"] and item["partition_equal"] for item in row["control_equivalence"].values()):
                    raise ValueError("Unchanged CPM control failed stage equivalence; variants not started")
            row["status"] = "replayed_unscored"
        for item in inputs:
            check(item)
        verify_runtime(runtime_path, frozen)
        report["status"] = "three_cpm_arms_replayed_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
