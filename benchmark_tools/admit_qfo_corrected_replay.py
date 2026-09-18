"""Independently admit corrected cached replay evidence, not benchmark accuracy."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_corrected_replay_stages import audit as audit_stages
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_replay import command_for, validate_admission, PLAN_SHA
from benchmark_tools.run_qfo_corrected_replay import validate_completion
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "188fde21860a70da55fee1358177485c970a11a3"
STAGES = ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]
FILENAMES = ["orthogroups_multipass.txt", "orthogroups_multipass_refined.txt",
             "orthogroups_profiles.txt", "orthogroups_profiles_refined.txt"]
HELPERS = ("run_qfo_corrected_replay.py", "checked_replay_payload_worker.py", "checked_replay_interceptor.py",
           "validate_checked_replay_payload.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py",
           "probe_leiden_boundary.py", "prepare_qfo_corrected_replay.py")


def validate_parent(parent, worker, replay, plan, plan_record, scheduler, root, executor, helpers):
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"
            or parent["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Replay scheduler identity/allocation differs")
    if (parent["status"] != "corrected_checked_replay_complete_pending_admission"
            or type(parent["exit_code"]) is not int or parent["exit_code"] != 0
            or parent["accuracy_evaluated"] is not False or worker["accuracy_evaluated"] is not False
            or parent["executor_commit"] != EXECUTOR):
        raise ValueError("Replay parent has not completed under the frozen executor")
    timestamps = [parent["started_epoch"], parent["finished_epoch"]]
    if (any(type(t) not in (int, float) or not math.isfinite(t) for t in timestamps)
            or not 0 < timestamps[0] <= timestamps[1]):
        raise ValueError("Invalid parent timestamps")
    if (parent["plan"] != plan_record or worker["plan"] != plan_record
            or parent["source"] != helpers[0] or worker["source"] != helpers[0]
            or parent["helpers"] != helpers):
        raise ValueError("Parent/worker source or plan binding differs")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    expected_worker = [sys.executable, str(executor / "benchmark_tools/run_qfo_corrected_replay.py"),
                       "--root", str(root), "--plan", plan_record["path"],
                       "--plan-sha256", plan_record["sha256"], "--replay-worker"]
    if parent["worker_command"] != expected_worker or replay["command"] != plan["native_command"]:
        raise ValueError("Replay or parent worker command differs")
    if (replay["cwd"] != str(launcher) or plan["cwd"] != str(launcher)
            or replay["source"] != worker["replay_source"]):
        raise ValueError("Replay source/cwd differs")
    parameters = {"accuracy_profile": "high_sensitivity", "cpm_resolution": .1, "profile_expansion": True,
                  "profile_iterations": 1, "jackknife_profile_thresholds": False,
                  "jackknife_single_copy_profiles": False, "profile_min_species": 1,
                  "matrix": "BLOSUM62", "leiden_seed": 4}
    if replay["parameters"] != parameters or plan["expected_stages"] != STAGES:
        raise ValueError("Scientific replay settings differ")
    validate_completion(worker, replay, STAGES)
    if type(replay["counts"]["species"]) is not int or replay["counts"]["species"] != 78:
        raise ValueError("Wrong corrected species universe")
    if any("official_orthobench" in stage for stage in replay["stages"]):
        raise ValueError("Unexpected evaluation during replay")


def validate_partitions(directory, replay, universe, native_record):
    coverage = []
    if [stage["label"] for stage in replay["stages"]] != STAGES:
        raise ValueError("Unexpected refined stage inventory")
    for stage, filename in zip(replay["stages"], FILENAMES):
        expected = directory / "replay" / filename
        if stage["output"]["path"] != str(expected):
            raise ValueError("Stage partition outside its frozen output path")
        check(stage["output"])
        groups = read_partition(expected, universe)
        if type(stage["clusters"]) is not int or stage["clusters"] != len(groups):
            raise ValueError("Stage group count differs")
        coverage.append({"label": stage["label"], "groups": len(groups), "output": stage["output"]})
        check(stage["output"])
    check(native_record)
    native = read_predictions(Path(native_record["path"]), "named_groups")
    if set(membership(native)) != universe:
        raise ValueError("Native comparison partition has wrong coverage")
    native_set, replay_set = {frozenset(g) for g in native.values()}, {frozenset(g) for g in groups}
    comparison = {"partition_equal": native_set == replay_set, "native_groups": len(native_set),
                  "replay_groups": len(replay_set), "native_only_groups": len(native_set - replay_set),
                  "replay_only_groups": len(replay_set - native_set)}
    check(native_record)
    return coverage, comparison


def admit(root, plan_path, plan_sha, report_sha, job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    # Scheduler termination is required before inspecting possibly partial output.
    plan, plan_record, admission_record, names_record = corrected_evidence(plan_path, plan_sha)
    directory = root / "benchmarks/results/qfo_corrected_checked_replay_v1"
    if plan["output_root"] != str(directory):
        raise ValueError("Unexpected corrected replay output root")
    executor = root / "benchmarks/work/publication_qfo_corrected_replay_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    if revision != EXECUTOR:
        raise ValueError("Replay executor revision differs")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PLAN_SHA)
    admission = json.loads(Path(admission_record["path"]).read_text())
    checkpoint, checkpoint_sha = validate_admission(admission, primary, plan["input_fastas"])
    if plan["primary_plan"] != record(primary_path) or plan["source"] != record(executor / "benchmark_tools/prepare_qfo_corrected_replay.py"):
        raise ValueError("Replay preparation source/primary plan differs")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    expected = command_for(Path(sys.executable), launcher, directory, Path(primary["input_directory"]), checkpoint, checkpoint_sha)
    environment = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if plan["native_command"] != expected or plan["environment_overrides"] != environment:
        raise ValueError("Replay command plan differs from frozen settings")
    parent = read_frozen(directory / "results.json", report_sha)
    worker_record, replay_record = record(directory / "checked_worker.json"), record(directory / "replay.json")
    if parent["worker"] != worker_record or parent["replay"] != replay_record:
        raise ValueError("Parent-linked output changed")
    worker = json.loads((directory / "checked_worker.json").read_text())
    replay = json.loads((directory / "replay.json").read_text())
    helpers = [record(executor / "benchmark_tools" / name) for name in HELPERS]
    validate_parent(parent, worker, replay, plan, plan_record, scheduler, root, executor, helpers)
    scientific_source = record(launcher / "benchmark_tools/replay_high_sensitivity.py")
    if replay["source"] != scientific_source:
        raise ValueError("Wrong scientific replay source")
    numeric = admission["content"]["numeric_checkpoint"]
    if (any(replay["input"][k] != numeric[k] for k in ("status", "summary", "manifest"))
            or replay["input"]["accuracy_evaluated"] is not False
            or replay["input"]["auditor"] != record(launcher / "benchmark_tools/audit_accuracy_checkpoint.py")):
        raise ValueError("Replay checkpoint evidence differs from native admission")
    core, runtime_path = root / "benchmarks/work/publication_method_native_v2", root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == parent["runtime_before"] == parent["runtime_after"]:
        raise ValueError("Frozen runtime observations differ")
    stages = audit_stages(root, executor, directory, worker, replay, plan_record)
    names = Path(names_record["path"]).read_text().splitlines()
    universe = set(names)
    if len(names) != 984137 or len(universe) != len(names):
        raise ValueError("Incomplete or duplicate corrected gene universe")
    native_record = admission["content"]["native_groups"]
    coverage, comparison = validate_partitions(directory, replay, universe, native_record)
    if coverage != parent["coverage"] or comparison != parent["native_partition_comparison"]:
        raise ValueError("Independent partition comparisons disagree with parent")
    records = [record(directory / "results.json"), worker_record, replay_record, plan_record, admission_record,
               names_record, native_record, scientific_source, record(primary_path), plan["source"],
               replay["input"]["auditor"], record(directory / "replay.log"), record(directory / "time.txt"),
               *helpers, *stages["checked_records"], *plan["checked_records"], *[s["output"] for s in replay["stages"]],
               *[record(Path(__file__).with_name(name)) for name in
                 ("admit_qfo_corrected_replay.py", "audit_corrected_replay_stages.py", "validate_checked_replay_payload.py",
                  "audit_historical_profile_ablation.py", "score_ygob_groups.py")]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting replay provenance records")
        unique[item["path"]] = item
    for item in unique.values():
        check(item)
    corrected_evidence(plan_path, plan_sha)
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during admission")
    result = {"status": "corrected_checked_replay_admitted", "accuracy_evaluated": False, "publication_ready": False,
              "source": record(__file__), "scheduler": scheduler, "accounting": accounting,
              "plan": plan_record, "source_report": record(directory / "results.json"),
              "clustering": stages["clustering"], "coverage": coverage, "native_partition_comparison": comparison,
              "checked_records": list(unique.values()),
              "limitations": ["Cached shared-host replay, not dedicated end-to-end timing.",
                  "Native/replay disagreement is retained without retries or score transfer.",
                  "Checks preserved native observations, not a complete historical process or memory trace.",
                  "Candidate preparation, reconciliation and independent benchmark assessment remain required."]}
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "plan", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("plan-sha256", "report-sha256", "job"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.plan.resolve(), args.plan_sha256, args.report_sha256, args.job, args.output.resolve())
