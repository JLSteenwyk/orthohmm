"""Run the frozen corrected-input replay through checked clustering workers."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def replay_worker(root, plan_path, plan_sha):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    executor = Path(__file__).resolve().parent.parent
    sys.path[0] = str(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from orthohmm import externals
    if Path(replay.__file__).resolve() != launcher / "benchmark_tools/replay_high_sensitivity.py":
        raise ValueError("Wrong frozen replay import")
    sys.path.insert(0, str(executor / "benchmark_tools"))
    from checked_replay_payload_worker import corrected_evidence
    from checked_replay_interceptor import CheckedReplaySubprocess
    from validate_checked_replay_payload import validate
    from prepare_ob_candidate_neighborhood import record
    plan, plan_record, _, _ = corrected_evidence(plan_path, plan_sha)
    output = Path(plan["output_root"])
    command = plan["native_command"]
    if command[:2] != [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py")]:
        raise ValueError("Unexpected replay executable/module")
    original = externals.subprocess
    proxy = CheckedReplaySubprocess(original, root, output / "clustering",
        executor / "benchmark_tools/checked_replay_payload_worker.py",
        lambda payload, manifest: validate(payload, manifest, root, executor, corrected_plan=plan_record),
        worker_args=["--corrected-plan", str(plan_path), "--corrected-plan-sha256", plan_sha])
    report = {"status": "running", "source": record(__file__), "plan": plan_record,
              "replay_source": record(replay.__file__), "calls": proxy.calls, "accuracy_evaluated": False}
    externals.subprocess = proxy
    try:
        sys.argv = command[1:]
        replay.main(command[2:])
        if len(proxy.calls) != 4 or any(row["status"] != "checked" for row in proxy.calls):
            raise ValueError("Require all four checked clustering stages")
        report["status"] = "corrected_checked_replay_worker_returned"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        externals.subprocess = original
        (output / "checked_worker.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def validate_completion(worker, replay, expected_stages):
    if (worker.get("status") != "corrected_checked_replay_worker_returned"
            or len(worker["calls"]) != 4
            or [c["stage"] for c in worker["calls"]] != ["initial", "multipass", "profile_base", "profile_expanded"]
            or any(c["status"] != "checked" or c["exit_code"] != 0 for c in worker["calls"])):
        raise ValueError("Incomplete checked replay worker")
    counts = replay["counts"]
    if (counts["genes"] != 984137 or type(counts.get("profiles_built")) is not int
            or counts["profiles_built"] <= 0
            or [s["label"] for s in replay["stages"]] != expected_stages):
        raise ValueError("Incomplete corrected replay or missing profile construction")


def run(root, plan_path, plan_sha):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.prepare_qfo_corrected_replay import command_for
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.score_ygob_groups import read_predictions, membership
    if os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled 32-CPU replay")
    plan, plan_record, admission_record, names_record = corrected_evidence(plan_path, plan_sha)
    output = Path(plan["output_root"])
    if output.exists():
        raise FileExistsError(output)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    fasta_dirs = {Path(r["path"]).parent for r in plan["input_fastas"]}
    if len(fasta_dirs) != 1:
        raise ValueError("Ambiguous corrected input directory")
    expected = command_for(Path(sys.executable), launcher, output, next(iter(fasta_dirs)),
        Path(plan["checkpoint_manifest"]["path"]).parent, plan["checkpoint_manifest"]["sha256"])
    if plan["native_command"] != expected or plan["cwd"] != str(launcher):
        raise ValueError("Replay command differs from frozen scientific settings")
    expected_env = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if plan["environment_overrides"] != expected_env:
        raise ValueError("Replay environment differs from frozen settings")
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Frozen runtime differs from command plan")
    executor = Path(__file__).resolve().parent.parent
    helpers = [record(executor / "benchmark_tools" / name) for name in (
        "run_qfo_corrected_replay.py", "checked_replay_payload_worker.py", "checked_replay_interceptor.py",
        "validate_checked_replay_payload.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py",
        "probe_leiden_boundary.py", "prepare_qfo_corrected_replay.py")]
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "source": record(__file__), "plan": plan_record,
        "helpers": helpers, "runtime_before": runtime, "job_id": os.environ["SLURM_JOB_ID"],
        "executor_commit": subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip(),
        "started_epoch": time.time(), "accuracy_evaluated": False}
    try:
        env = {**os.environ, **plan["environment_overrides"]}
        command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root),
                   "--plan", str(plan_path), "--plan-sha256", plan_sha, "--replay-worker"]
        report["worker_command"] = command
        with (output / "replay.log").open("x") as log:
            result = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *command],
                cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = result.returncode
        if result.returncode:
            raise RuntimeError("Corrected replay failed; no implicit retry")
        worker = json.loads((output / "checked_worker.json").read_text())
        replay = json.loads((output / "replay.json").read_text())
        validate_completion(worker, replay, plan["expected_stages"])
        universe = set(Path(names_record["path"]).read_text().splitlines())
        coverage = []
        for stage in replay["stages"]:
            check(stage["output"])
            groups = read_partition(Path(stage["output"]["path"]), universe)
            coverage.append({"label": stage["label"], "groups": len(groups), "output": stage["output"]})
        admission = json.loads(Path(admission_record["path"]).read_text())
        native_record = admission["content"]["native_groups"]
        check(native_record)
        native_groups = read_predictions(Path(native_record["path"]), "named_groups")
        if set(membership(native_groups)) != universe:
            raise ValueError("Native comparison groups differ from admitted universe")
        native_set = {frozenset(g) for g in native_groups.values()}
        replay_set = {frozenset(g) for g in groups}
        report["native_partition_comparison"] = {"partition_equal": native_set == replay_set,
            "native_groups": len(native_set), "replay_groups": len(replay_set),
            "native_only_groups": len(native_set - replay_set), "replay_only_groups": len(replay_set - native_set)}
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during replay")
        corrected_evidence(plan_path, plan_sha)
        for item in helpers:
            check(item)
        report.update(status="corrected_checked_replay_complete_pending_admission", coverage=coverage,
            worker=record(output / "checked_worker.json"), replay=record(output / "replay.json"),
            limitations=["Incremental shared-host cached replay; not end-to-end matched timing.",
                         "Nonidentical native/replay partitions require separate scoring, never score transfer.",
                         "Independent admission, candidate preparation, reconciliation and assessment remain required."])
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
    parser.add_argument("--replay-worker", action="store_true")
    args = parser.parse_args()
    (replay_worker if args.replay_worker else run)(args.root.resolve(), args.plan.resolve(), args.plan_sha256)
