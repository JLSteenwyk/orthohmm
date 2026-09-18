"""Execute a planned sequence-control replay through checked clustering workers."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def validate_completion(worker, metrics, variant):
    if (worker["status"] != "sequence_checked_worker_returned"
            or [call["stage"] for call in worker["calls"]] != variant["expected_clustering_calls"]
            or any(call["status"] != "checked" or call["exit_code"] != 0 for call in worker["calls"])
            or metrics["parameters"]["profile_expansion"] is not False
            or [stage["label"] for stage in metrics["stages"]] != variant["expected_stages"]
            or metrics["counts"]["genes"] != variant["expected_genes"]
            or metrics["counts"]["species"] != variant["expected_species"]):
        raise ValueError("Incomplete or unexpected sequence-control replay")


def replay_worker(root, path, sha, label):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    executor = Path(__file__).resolve().parent.parent
    sys.path[0] = str(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from orthohmm import externals
    if Path(replay.__file__).resolve() != launcher / "benchmark_tools/replay_high_sensitivity.py":
        raise ValueError("Wrong frozen replay import")
    sys.path.insert(0, str(executor / "benchmark_tools"))
    from sequence_graph_evidence import sequence_evidence
    from checked_replay_interceptor import CheckedReplaySubprocess
    from validate_checked_replay_payload import validate
    from prepare_ob_candidate_neighborhood import record
    plan, plan_record, _, _ = sequence_evidence(path, sha, label)
    variant = plan["variants"][label]
    output = Path(variant["output_root"])
    original = externals.subprocess
    proxy = CheckedReplaySubprocess(original, root, output / "clustering",
        executor / "benchmark_tools/checked_replay_payload_worker.py",
        lambda payload, manifest: validate(payload, manifest, root, executor,
                                          sequence_plan=plan_record, sequence_variant=label),
        worker_args=["--sequence-plan", str(path), "--sequence-plan-sha256", sha, "--sequence-variant", label])
    report = dict(status="running", source=record(__file__), plan=plan_record, variant=label,
                  replay_source=record(replay.__file__), calls=proxy.calls, accuracy_evaluated=False)
    externals.subprocess = proxy
    try:
        sys.argv = variant["native_command"][1:]
        replay.main(variant["native_command"][2:])
        report["status"] = "sequence_checked_worker_returned"
        validate_completion(report, json.loads((output / "replay.json").read_text()), variant)
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        externals.subprocess = original
        (output / "checked_worker.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def run(root, path, sha, label):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.sequence_graph_evidence import sequence_evidence
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    plan, plan_record, _, names_record = sequence_evidence(path, sha, label)
    variant = plan["variants"][label]
    if (os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID")
            or os.environ.get("SLURM_MEM_PER_NODE") != str(variant["requested_memory_gib"] * 1024)):
        raise ValueError("Require scheduled CPU/memory allocation matching the reviewed plan")
    output = Path(variant["output_root"])
    if output.exists():
        raise FileExistsError(output)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    if plan["cwd"] != str(launcher):
        raise ValueError("Wrong frozen working directory")
    environment = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if plan["environment_overrides"] != environment:
        raise ValueError("Wrong graph environment")
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Frozen runtime changed")
    executor = Path(__file__).resolve().parent.parent
    helpers = [record(executor / "benchmark_tools" / name) for name in (
        "sequence_graph_evidence.py", "checked_replay_interceptor.py", "checked_replay_payload_worker.py",
        "validate_checked_replay_payload.py", "repeat_qfo_saved_graph.py", "checked_python_pair_worker.py",
        "probe_leiden_boundary.py", "run_sequence_graph_control.py")]
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="running", source=record(__file__), helpers=helpers, plan=plan_record,
        variant=label, job_id=os.environ["SLURM_JOB_ID"], started_epoch=time.time(),
        runtime_before=runtime, accuracy_evaluated=False,
        executor_commit=subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip())
    try:
        command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--plan", str(path),
                   "--plan-sha256", sha, "--variant", label, "--replay-worker"]
        report["worker_command"] = command
        with (output / "replay.log").open("x") as log:
            process = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *command],
                cwd=launcher, env={**os.environ, **environment}, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        process.check_returncode()
        metrics = json.loads((output / "replay.json").read_text())
        worker = json.loads((output / "checked_worker.json").read_text())
        validate_completion(worker, metrics, variant)
        if metrics["command"] != variant["native_command"] or metrics["cwd"] != str(launcher):
            raise ValueError("Native replay command differs")
        universe = set(Path(names_record["path"]).read_text().splitlines())
        if len(universe) != variant["expected_genes"]:
            raise ValueError("Sequence universe count differs")
        coverage = []
        for stage in metrics["stages"]:
            expected = output / "replay" / ("orthogroups_" + stage["label"] + ".txt")
            if stage["output"]["path"] != str(expected):
                raise ValueError("Stage output outside planned location")
            check(stage["output"])
            groups = read_partition(expected, universe)
            if len(groups) != stage["clusters"]:
                raise ValueError("Stage group count differs")
            coverage.append(dict(label=stage["label"], groups=len(groups), output=stage["output"]))
        sequence_evidence(path, sha, label)
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during replay")
        for item in [report["source"], *helpers]:
            check(item)
        report.update(status="sequence_checked_graph_complete_pending_admission", coverage=coverage,
            worker=record(output / "checked_worker.json"), replay=record(output / "replay.json"),
            limitations=["Independent native-output admission and scoring remain required.",
                         "Shared-host incremental graph replay is not dedicated end-to-end timing.",
                         "All-hit and top100 outcomes must be retained separately, including failures."])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "plan"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--variant", required=True, choices=("all_hits", "top100"))
    parser.add_argument("--replay-worker", action="store_true")
    args = parser.parse_args()
    (replay_worker if args.replay_worker else run)(args.root.resolve(), args.plan.resolve(), args.plan_sha256, args.variant)
