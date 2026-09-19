"""Independently validate the unchanged CPM replay before changed-arm execution."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_corrected_replay_stages import audit
from benchmark_tools.admit_qfo_corrected_replay import STAGES, FILENAMES
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_corrected_replay import validate_completion
from benchmark_tools.run_qfo_cpm_control import BASELINE_SHA, compare_partitions
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "ef5fcd6f9d0ed46d8facc67631dcc7d3d3bb8f3e"
JOB = "21956"
HELPERS = ("cpm_replay_context.py", "checked_replay_payload_worker.py", "checked_replay_interceptor.py",
           "validate_checked_replay_payload.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py",
           "probe_leiden_boundary.py", "run_qfo_corrected_replay.py", "audit_historical_profile_ablation.py")


def validate_parent(parent, worker, replay, context, scheduler, root, executor, source, helpers, inputs):
    if (scheduler["JobIDRaw"] != JOB or scheduler["State"] != "COMPLETED"
            or scheduler["ExitCode"] != "0:0" or scheduler["NodeList"] != "bizon"
            or scheduler["AllocCPUS"] != "32" or parent["job_id"] != JOB):
        raise ValueError("Control scheduler identity/allocation differs")
    if (parent["status"] != "cpm_control_reproduced_pending_independent_admission"
            or type(parent["exit_code"]) is not int or parent["exit_code"] != 0
            or parent["executor_commit"] != EXECUTOR
            or parent["accuracy_evaluated"] is not False or worker["accuracy_evaluated"] is not False
            or parent["changed_arms_authorized"] is not False):
        raise ValueError("Control parent has not completed under frozen conditions")
    timestamps = [parent["started_epoch"], parent["finished_epoch"]]
    if (any(type(t) not in (int, float) or not math.isfinite(t) for t in timestamps)
            or not 0 < timestamps[0] <= timestamps[1]):
        raise ValueError("Invalid control timestamps")
    if (parent["context"] != context or worker["context"] != context
            or parent["source"] != source or worker["source"] != source
            or parent["helpers"] != helpers or parent["checked_inputs"] != inputs):
        raise ValueError("Control context/source/input binding differs")
    expected = [sys.executable, str(executor / "benchmark_tools/run_qfo_cpm_control.py"),
                "--root", str(root), "--worker"]
    if parent["worker_command"] != expected or replay["command"] != context["native_command"]:
        raise ValueError("Control command differs")
    if replay["cwd"] != context["cwd"] or replay["source"] != worker["replay_source"]:
        raise ValueError("Control replay source/cwd differs")
    parameters = {"accuracy_profile": "high_sensitivity", "cpm_resolution": .1, "profile_expansion": True,
                  "profile_iterations": 1, "jackknife_profile_thresholds": False,
                  "jackknife_single_copy_profiles": False, "profile_min_species": 1,
                  "matrix": "BLOSUM62", "leiden_seed": 4}
    if replay["parameters"] != parameters or context["expected_stages"] != STAGES:
        raise ValueError("Control scientific settings differ")
    validate_completion(worker, replay, STAGES)
    if type(replay["counts"]["species"]) is not int or replay["counts"]["species"] != 78:
        raise ValueError("Wrong control species universe")
    for stage, filename in zip(replay["stages"], FILENAMES):
        if (stage["output"]["path"] != str(Path(context["output_root"]) / "replay" / filename)
                or "official_orthobench" in stage):
            raise ValueError("Wrong control stage path or unexpected evaluation")


def admit(root, destination):
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, JOB)
    # Do not inspect outputs while the scheduled producer may still write them.
    local_helpers = [record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, admission_record, names_record = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, "control")
    directory = Path(context["output_root"])
    executor = root / "benchmarks/work/publication_qfo_cpm_control_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Control executor revision differs")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    baseline_path = root / "benchmark_tools/results/qfo_corrected_replay_admission_21757.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    if baseline["status"] != "corrected_checked_replay_admitted" or baseline["plan"] != plan_record:
        raise ValueError("Wrong corrected baseline admission")
    inputs = [record(baseline_path), names_record, *context["checked_records"], *baseline["checked_records"]]
    for item in inputs:
        check(item)
    parent_record, worker_record, replay_record = [record(directory / name) for name in
        ("results.json", "checked_worker.json", "replay.json")]
    parent, worker, replay = [read_frozen(Path(r["path"]), r["sha256"]) for r in
        (parent_record, worker_record, replay_record)]
    if parent["worker"] != worker_record or parent["replay"] != replay_record:
        raise ValueError("Control output changed after parent completion")
    source = record(executor / "benchmark_tools/run_qfo_cpm_control.py")
    helpers = [record(executor / "benchmark_tools" / name) for name in HELPERS]
    validate_parent(parent, worker, replay, context, scheduler, root, executor, source, helpers, inputs)
    launcher = Path(context["cwd"])
    scientific_source = record(launcher / "benchmark_tools/replay_high_sensitivity.py")
    if replay["source"] != scientific_source:
        raise ValueError("Wrong scientific source")
    admission = read_frozen(Path(admission_record["path"]), admission_record["sha256"])
    numeric = admission["content"]["numeric_checkpoint"]
    if (any(replay["input"][k] != numeric[k] for k in ("status", "summary", "manifest"))
            or replay["input"]["accuracy_evaluated"] is not False
            or replay["input"]["auditor"] != record(launcher / "benchmark_tools/audit_accuracy_checkpoint.py")):
        raise ValueError("Wrong control checkpoint evidence")
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == parent["runtime_before"] == parent["runtime_after"]:
        raise ValueError("Control runtime differs")
    stages = audit(root, executor, directory, worker, replay, plan_record, cpm_arm="control")
    names = Path(names_record["path"]).read_text().splitlines()
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong control gene universe")
    comparisons = compare_partitions(replay, baseline, set(names))
    if comparisons != parent["stage_comparisons"] or not all(r["partition_equal"] for r in comparisons):
        raise ValueError("Independent control partition comparison disagrees")
    if any(type(s["clusters"]) is not int or s["clusters"] != r["groups"]
           for s, r in zip(replay["stages"], comparisons)):
        raise ValueError("Control stage count differs")
    records = [parent_record, worker_record, replay_record, plan_record, admission_record, source,
        scientific_source, replay["input"]["auditor"], record(directory / "replay.log"), record(directory / "time.txt"),
        *inputs, *helpers, *local_helpers, *stages["checked_records"], *[r["output"] for r in comparisons]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting control provenance")
        unique[item["path"]] = item
        check(item)
    corrected_evidence(plan_path, REPLAY_SHA)
    evidence(root, plan, plan_record, "control")
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during control admission")
    result = {"status": "cpm_control_reproduced_and_admitted", "source": record(__file__),
        "scheduler": scheduler, "accounting": accounting, "context": context, "source_report": parent_record,
        "clustering": stages["clustering"], "stage_comparisons": comparisons,
        "checked_records": list(unique.values()), "accuracy_evaluated": False, "publication_ready": False,
        "changed_arms_authorized": ["cpm_low", "cpm_high"],
        "limitations": ["Authorizes only the two prespecified CPM replay experiments, not their results.",
            "Candidate construction, phylogeny, conversion and scoring require separate validation.",
            "Cached shared-host replay is not controlled end-to-end timing."]}
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
