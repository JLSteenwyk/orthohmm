"""Repeat the unchanged corrected CPM replay before authorizing changed arms."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

BASELINE_SHA = "1eec5ffb675fff234e5ce0db7e65abfa9bec09bdea8f5bd13ca72ad72d7ec40e"


def worker(root, arm="control"):
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
    from cpm_replay_context import evidence, REPLAY_SHA
    from prepare_ob_candidate_neighborhood import record
    from validate_checked_replay_payload import validate
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, _, _ = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, arm)
    output = Path(context["output_root"])
    original = externals.subprocess
    proxy = CheckedReplaySubprocess(original, root, output / "clustering",
        executor / "benchmark_tools/checked_replay_payload_worker.py",
        lambda payload, manifest: validate(payload, manifest, root, executor,
            corrected_plan=plan_record, cpm_arm=arm),
        worker_args=["--corrected-plan", str(plan_path), "--corrected-plan-sha256", REPLAY_SHA, "--cpm-arm", arm])
    report = {"status": "running", "source": record(__file__), "context": context,
              "replay_source": record(replay.__file__), "calls": proxy.calls, "accuracy_evaluated": False}
    externals.subprocess = proxy
    try:
        sys.argv = context["native_command"][1:]
        replay.main(context["native_command"][2:])
        if len(proxy.calls) != 4 or any(row["status"] != "checked" for row in proxy.calls):
            raise ValueError("Require four checked clustering stages")
        report["status"] = "corrected_checked_replay_worker_returned"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        externals.subprocess = original
        (output / "checked_worker.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def compare_partitions(replay, baseline, universe):
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.prepare_ob_candidate_neighborhood import check
    if [s["label"] for s in replay["stages"]] != [s["label"] for s in baseline["coverage"]] or len(replay["stages"]) != 4:
        raise ValueError("Require all four baseline stages in order")
    comparisons = []
    for current, prior in zip(replay["stages"], baseline["coverage"]):
        check(current["output"])
        check(prior["output"])
        groups = read_partition(Path(current["output"]["path"]), universe)
        old_groups = read_partition(Path(prior["output"]["path"]), universe)
        equal = {frozenset(g) for g in groups} == {frozenset(g) for g in old_groups}
        comparisons.append({"label": current["label"], "partition_equal": equal,
            "bytes_equal": all(current["output"][k] == prior["output"][k] for k in ("bytes", "sha256")),
            "groups": len(groups), "baseline_groups": len(old_groups), "output": current["output"], "baseline": prior["output"]})
    return comparisons


def run(root):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.run_qfo_corrected_replay import validate_completion
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"):
        raise ValueError("Require scheduled 32-CPU bizon control replay")
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, _, names_record = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, "control")
    output = Path(context["output_root"])
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    baseline_path = root / "benchmark_tools/results/qfo_corrected_replay_admission_21757.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    if baseline["status"] != "corrected_checked_replay_admitted" or baseline["plan"] != plan_record:
        raise ValueError("Wrong corrected baseline admission")
    inputs = [record(baseline_path), names_record, *context["checked_records"], *baseline["checked_records"]]
    for item in inputs:
        check(item)
    launcher = Path(context["cwd"])
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Changed frozen scientific runtime")
    executor = Path(__file__).resolve().parent.parent
    helpers = [record(executor / "benchmark_tools" / name) for name in (
        "cpm_replay_context.py", "checked_replay_payload_worker.py", "checked_replay_interceptor.py",
        "validate_checked_replay_payload.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py",
        "probe_leiden_boundary.py", "run_qfo_corrected_replay.py", "audit_historical_profile_ablation.py")]
    source = record(__file__)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "source": source, "helpers": helpers, "checked_inputs": inputs,
        "context": context, "runtime_before": runtime, "job_id": os.environ["SLURM_JOB_ID"],
        "executor_commit": subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip(),
        "started_epoch": time.time(), "accuracy_evaluated": False, "changed_arms_authorized": False}
    try:
        command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--worker"]
        report["worker_command"] = command
        with (output / "replay.log").open("x") as log:
            process = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *command],
                cwd=launcher, env={**os.environ, **context["environment_overrides"]}, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        if process.returncode:
            raise RuntimeError("Control replay failed; no implicit retry")
        checked_worker = json.loads((output / "checked_worker.json").read_text())
        replay = json.loads((output / "replay.json").read_text())
        if checked_worker["context"] != context:
            raise ValueError("Worker CPM context changed")
        validate_completion(checked_worker, replay, context["expected_stages"])
        universe = set(Path(names_record["path"]).read_text().splitlines())
        comparisons = compare_partitions(replay, baseline, universe)
        report["stage_comparisons"] = comparisons
        if not all(row["partition_equal"] for row in comparisons):
            raise ValueError("Control does not reproduce all four baseline partitions")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during control replay")
        corrected_evidence(plan_path, REPLAY_SHA)
        evidence(root, plan, plan_record, "control")
        for item in [source, *helpers, *inputs, *[r["output"] for r in comparisons]]:
            check(item)
        report.update(status="cpm_control_reproduced_pending_independent_admission",
            worker=record(output / "checked_worker.json"), replay=record(output / "replay.json"),
            limitations=["Control repeat only; changed CPM arms remain unauthorized pending independent admission.",
                         "Incremental shared-host cached runtime is not controlled end-to-end timing."])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    (worker if args.worker else run)(args.root.resolve())
