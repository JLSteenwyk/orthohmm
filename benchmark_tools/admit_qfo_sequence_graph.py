"""Admit one terminal sequence-control graph run from retained native evidence."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.audit_accuracy_checkpoint import audit as audit_checkpoint
from benchmark_tools.audit_corrected_replay_stages import audit as audit_stages
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.compare_qfo_search_coverage import frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_sequence_graph import validate_completion
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.sequence_graph_evidence import sequence_evidence
from benchmark_tools.verify_qfo_replay_launcher import verify

EXECUTOR = "69adea4d5a1dc634783b20ff7bbe9fbe1c6464db"
HELPERS = ("sequence_graph_evidence.py", "checked_replay_interceptor.py", "checked_replay_payload_worker.py",
           "validate_checked_replay_payload.py", "repeat_qfo_saved_graph.py", "checked_python_pair_worker.py",
           "probe_leiden_boundary.py", "run_sequence_graph_control.py")


def validate_parent(parent, worker, replay, plan, plan_record, label, scheduler, root, executor):
    variant = plan["variants"][label]
    source = record(executor / "benchmark_tools/run_qfo_sequence_graph.py")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"
            or scheduler["ReqMem"] not in (str(variant["requested_memory_gib"]) + "G",
                                           str(variant["requested_memory_gib"]) + "Gn")
            or parent["status"] != "sequence_checked_graph_complete_pending_admission"
            or type(parent["exit_code"]) is not int or parent["exit_code"] != 0
            or parent["job_id"] != scheduler["JobIDRaw"] or parent["executor_commit"] != EXECUTOR
            or parent["variant"] != label or worker["variant"] != label
            or parent["accuracy_evaluated"] is not False or worker["accuracy_evaluated"] is not False
            or parent["source"] != source or worker["source"] != source
            or parent["plan"] != plan_record or worker["plan"] != plan_record
            or parent["helpers"] != [record(executor / "benchmark_tools" / name) for name in HELPERS]):
        raise ValueError("Parent/worker execution identity differs")
    times = (parent["started_epoch"], parent["finished_epoch"])
    if any(type(t) not in (int, float) or not math.isfinite(t) for t in times) or not 0 < times[0] <= times[1]:
        raise ValueError("Invalid execution timestamps")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    command = [sys.executable, str(executor / "benchmark_tools/run_qfo_sequence_graph.py"),
               "--root", str(root), "--plan", plan_record["path"], "--plan-sha256", plan_record["sha256"],
               "--variant", label, "--replay-worker"]
    environment = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if (parent["worker_command"] != command or replay["command"] != variant["native_command"]
            or replay["cwd"] != str(launcher) or plan["cwd"] != str(launcher)
            or plan["environment_overrides"] != environment
            or replay["source"] != record(launcher / "benchmark_tools/replay_high_sensitivity.py")
            or worker["replay_source"] != replay["source"]):
        raise ValueError("Scientific command/source/environment differs")
    expected = dict(accuracy_profile="high_sensitivity", cpm_resolution=.1, profile_expansion=False,
                    profile_iterations=1, jackknife_profile_thresholds=False, jackknife_single_copy_profiles=False,
                    profile_min_species=1, matrix="BLOSUM62", leiden_seed=4)
    if replay["parameters"] != expected or any("official_orthobench" in stage for stage in replay["stages"]):
        raise ValueError("Unexpected scientific settings or scoring during inference")
    validate_completion(worker, replay, variant)


def partitions(directory, replay, universe):
    if [stage["label"] for stage in replay["stages"]] != ["multipass", "multipass_refined"]:
        raise ValueError("Wrong output stage inventory")
    coverage = []
    for stage in replay["stages"]:
        path = directory / "replay" / ("orthogroups_" + stage["label"] + ".txt")
        if stage["output"]["path"] != str(path):
            raise ValueError("Stage output outside planned location")
        check(stage["output"])
        groups = read_partition(path, universe)
        if type(stage["clusters"]) is not int or stage["clusters"] != len(groups):
            raise ValueError("Stage group count differs")
        coverage.append(dict(label=stage["label"], groups=len(groups), output=stage["output"]))
        check(stage["output"])
    return coverage


def admit(root, plan_path, plan_sha, label, job, report_sha, destination):
    if destination.exists():
        raise FileExistsError(destination)
    plan, plan_record, admission_record, names_record = sequence_evidence(plan_path, plan_sha, label)
    variant = plan["variants"][label]
    scheduler = completed(job, 32, str(variant["requested_memory_gib"]) + "G")
    executor = frozen(root, "publication_qfo_sequence_graph_v1", EXECUTOR)
    if plan["source"] != record(executor / "benchmark_tools/prepare_qfo_sequence_graph.py"):
        raise ValueError("Plan preparation source differs")
    directory = Path(variant["output_root"])
    if directory != root / "benchmarks/results/qfo_sequence_graph_v1" / label:
        raise ValueError("Unexpected sequence output root")
    parent = read_frozen(directory / "results.json", report_sha)
    worker_record, replay_record = record(directory / "checked_worker.json"), record(directory / "replay.json")
    if parent["worker"] != worker_record or parent["replay"] != replay_record:
        raise ValueError("Parent-linked native records changed")
    worker = json.loads(Path(worker_record["path"]).read_text())
    replay = json.loads(Path(replay_record["path"]).read_text())
    validate_parent(parent, worker, replay, plan, plan_record, label, scheduler, root, executor)
    checkpoint = Path(variant["checkpoint_manifest"]["path"]).parent
    numeric = audit_checkpoint(checkpoint, variant["checkpoint_manifest"]["sha256"])
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if (any(replay["input"][key] != numeric[key] for key in ("status", "summary", "manifest"))
            or replay["input"]["accuracy_evaluated"] is not False
            or replay["input"]["auditor"] != record(launcher / "benchmark_tools/audit_accuracy_checkpoint.py")):
        raise ValueError("Replay numeric input differs")
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == parent["runtime_before"] == parent["runtime_after"]:
        raise ValueError("Runtime observations differ")
    clustering = audit_stages(root, executor, directory, worker, replay, plan_record, sequence_variant=label)
    names = Path(names_record["path"]).read_text().splitlines()
    universe = set(names)
    if len(names) != 984137 or len(universe) != len(names):
        raise ValueError("Incomplete or duplicate corrected gene universe")
    coverage = partitions(directory, replay, universe)
    if coverage != parent["coverage"]:
        raise ValueError("Post-run coverage differs from parent")
    records = [record(directory / "results.json"), worker_record, replay_record, plan_record, admission_record,
        names_record, replay["source"], replay["input"]["auditor"], parent["source"], *parent["helpers"],
        record(directory / "replay.log"), record(directory / "time.txt"), *plan["checked_records"],
        *clustering["checked_records"], *[stage["output"] for stage in replay["stages"]],
        *[record(Path(__file__).with_name(name)) for name in ("admit_qfo_sequence_graph.py",
            "audit_corrected_replay_stages.py", "validate_checked_replay_payload.py", "sequence_graph_evidence.py")]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting provenance records")
        unique[item["path"]] = item
    for item in unique.values():
        check(item)
    sequence_evidence(plan_path, plan_sha, label)
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during admission")
    result = dict(status="corrected_sequence_graph_admitted", variant=label, accuracy_evaluated=False,
        publication_ready=False, source=record(__file__), scheduler=scheduler, plan=plan_record,
        source_report=record(directory / "results.json"), numeric=numeric, coverage=coverage,
        clustering=clustering["clustering"], checked_records=list(unique.values()),
        prediction=coverage[-1]["output"],
        limitations=["Retained native observations verified, not a reconstruction of historical process memory.",
                     "Incremental shared-host replay is not matched end-to-end timing.",
                     "Benchmark conversion, scoring and paired uncertainty remain separate requirements."])
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "plan", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("plan-sha256", "job", "report-sha256"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--variant", required=True, choices=("all_hits", "top100"))
    args = parser.parse_args()
    admit(args.root.resolve(), args.plan.resolve(), args.plan_sha256, args.variant, args.job,
          args.report_sha256, args.output.resolve())
