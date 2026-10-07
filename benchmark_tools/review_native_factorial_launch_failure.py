"""Dispose of a verified pre-native CPU-binding failure without retry or scoring."""

import argparse
import csv
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import (
    bind_session, native_command, runtime_review, scheduler_fields,
)
from benchmark_tools.run_native_factorial_cost import (
    ROOT, SCOPE, read, validate_plan, validate_request, verify_terminal,
)
from benchmark_tools.validate_native_factorial_outputs import Evidence, require


ACCOUNT_FIELDS = ("JobIDRaw", "State", "ExitCode", "ElapsedRaw", "AllocCPUS", "ReqMem", "NodeList")
MEASUREMENT_FILES = {"command.json", "step.log", "host_processes.jsonl", "go.json", "release.json"}
ABSENT_MEASUREMENT = ("ready.json", "done.json", "native.log", "lineage_report.json")
ABSENT_NATIVE = ("native_execution.json", "metrics.json", "native")


def accounting_rows(raw):
    rows = list(csv.reader(raw.strip().splitlines(), delimiter="|"))
    require(len(rows) == 3 and all(len(row) == len(ACCOUNT_FIELDS) for row in rows),
            "Require exactly the allocation, batch and one failed launch step")
    records = [dict(zip(ACCOUNT_FIELDS, row)) for row in rows]
    require(len({row["JobIDRaw"] for row in records}) == 3, "Duplicate accounting identity")
    return records


def failure_scope(request, run, result, verification, terminal, command, step_log, abort, cleanup, rows):
    job = request["job_id"]
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields.get("ExitCode") == "1:0",
            "Require the retained failed allocation, not a successful native result")
    require(result.get("status") == "factorial_attempt_failed_retained"
        and verification.get("status") == "verified_wrapper_failed"
        and verification.get("error_type") == "TimeoutError"
        and verification.get("error") == str(Path(run["output_root"]) / "measurement/ready.json"),
        "Require a failed wrapper waiting for the unlaunched worker")
    require(command.get("cpus") == 32 and type(command.get("cpus")) is int
        and command.get("timeout_s") == 85800 and type(command.get("timeout_s")) is int
        and type(command.get("interval_s")) in (int, float)
        and command.get("interval_s") == 1., "Frozen command limits differ")
    require(set(abort) == {"abort"} and abort["abort"] is True
        and set(cleanup) == {"release"} and cleanup["release"] is True,
            "Require original abort and cleanup gates, never a native go decision")
    lines = step_log.strip().splitlines()
    require(len(lines) == 4, "Unclassified additional launch diagnostics")
    match = re.fullmatch(
        r"srun: error: CPU binding outside of job step allocation, allocated CPUs are: (0x[0-9A-Fa-f]+)\.",
        lines[0])
    require(match is not None and lines[1:] == [
        f"srun: error: Task launch for StepId={job}.0 failed on node bizon: Unable to satisfy cpu bind request",
        "srun: error: Application launch failed: Unable to satisfy cpu bind request",
        "srun: Job step aborted"], "Require the actual classified Slurm binding refusal")
    mask = int(match[1], 16)
    require(mask > 0 and mask.bit_length() <= 192, "Invalid host allocation mask")
    allocated = [cpu for cpu in range(mask.bit_length()) if mask & (1 << cpu)]
    missing = sorted(set(range(32)) - set(allocated))
    require(missing and len(allocated) == 64, "Require a real mismatch within the frozen 64-slot allocation")
    by_id = {row["JobIDRaw"]: row for row in rows}
    require(len(rows) == 3 and set(by_id) == {str(job), f"{job}.batch", f"{job}.0"},
            "Another job or step is not this failed attempt")
    for identity, row in by_id.items():
        expected = ("CANCELLED", "0:64") if identity.endswith(".0") else ("FAILED", "1:0")
        require((row["State"], row["ExitCode"]) == expected and row["AllocCPUS"] == "64"
            and row["NodeList"] == "bizon" and re.fullmatch(r"[0-9]+", row["ElapsedRaw"]) is not None,
            "Actual failed launch/accounting differs")
    require(by_id[str(job)]["ReqMem"] == "128G" and by_id[f"{job}.0"]["ElapsedRaw"] == "0",
            "Wrong allocation or evidence of native step execution")
    return dict(failure_class="scheduler_rejected_frozen_affinity_before_native_release",
        allocated_os_cpu_mask=match[1], allocated_os_cpu_ids=allocated,
        requested_native_cpu_ids=list(range(32)), missing_requested_cpu_ids=missing,
        allocation_elapsed_s=int(by_id[str(job)]["ElapsedRaw"]), native_outcome="not_started")


def absent_artifacts(root):
    measurement = root / "measurement"
    absent = [measurement / name for name in ABSENT_MEASUREMENT] + [root / name for name in ABSENT_NATIVE]
    require(all(not path.exists() and not path.is_symlink() for path in absent),
            "Native or resource artifacts contradict a pre-native refusal")
    entries = list(measurement.iterdir())
    require({p.name for p in entries} == MEASUREMENT_FILES
        and all(p.is_file() and not p.is_symlink() for p in entries)
        and (measurement / "host_processes.jsonl").stat().st_size == 0,
        "Unexpected measurement inventory or started host monitor")
    return [str(path) for path in absent]


def review(request_ref, destination):
    request = read(request_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    run = validate_plan(plan)[request["index"]]
    validate_request(request, plan_ref, request["job_id"])
    terminal = verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Actual request comment differs")
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(root)
        and not destination.is_relative_to(session) and not destination.exists()
        and not destination.is_symlink(), "Require a fresh separate failure disposition")
    evidence = Evidence()
    result = evidence.json(session / "result.json")
    verification = evidence.json(root / "verification.json")
    bind_session(result, verification, request_ref, request, plan_ref, run)
    command = evidence.json(root / "measurement/command.json")
    require(command["command"] == native_command(plan_ref, run, read(plan["baseline"])),
            "Native command differs from the frozen request")
    step_path = evidence.bind(root / "measurement/step.log")
    abort = evidence.json(root / "measurement/go.json")
    cleanup = evidence.json(root / "measurement/release.json")
    account_command = ["sacct", "-j", str(request["job_id"]), "-n", "-P", "--format=" + ",".join(ACCOUNT_FIELDS)]
    account = subprocess.run(account_command, capture_output=True, text=True, check=True, timeout=5)
    rows = accounting_rows(account.stdout)
    failure = failure_scope(request, run, result, verification, terminal, command,
                            step_path.read_text(), abort, cleanup, rows)
    absent = absent_artifacts(root)
    evidence.bind(root / "measurement/host_processes.jsonl")
    source = record(__file__)
    for ref in [request_ref, plan_ref, request["policy"], plan["baseline"], source,
                *request["evidence"], *plan["helper_sources"], *plan["evidence"]]:
        check(ref)
        evidence.bind(ref["path"])
    runtime = runtime_review(plan, run, session, verification, evidence)
    require(absent_artifacts(root) == absent, "Failure artifact inventory changed")
    checked = evidence.finish()
    check(source)
    destination.mkdir(parents=True, exist_ok=False)
    save(destination / "runtime.json", runtime)
    runtime_ref = record(destination / "runtime.json")
    report = dict(schema="native_factorial_launch_failure_review_v1",
        status="pre_native_cpu_binding_failure_reviewed_retained", job_id=request["job_id"],
        index=run["index"], dataset=run["dataset"], cell=run["cell"], repeat=run["repeat"],
        request=request_ref, plan=plan_ref, source=source, scheduler=terminal,
        scheduler_state=fields.get("JobState", fields.get("State")), scheduler_exit_code=fields["ExitCode"],
        step_accounting=dict(command=account_command, stdout=account.stdout, stderr=account.stderr, rows=rows),
        runtime=runtime_ref, evidence=checked, absent_artifacts=absent, **failure,
        terminal_reviewed=True, next_identity_authorized=True,
        next_identity_scope="Next different frozen identity only after all unchanged fresh gates; never retry this attempt.",
        native_inference_started=False, native_outputs_validated=False, accuracy_evaluated=False,
        native_command_success=False, scheduler_success=False, resources=None,
        primary_resources_replayed=False, shared_host_resources_reviewed=False,
        scientific_timings_admitted=False, eligible_for_timing_comparison=False,
        original_receipts_rewritten=False, inference_reexecuted=False, automatic_retry=False,
        execution_scope=SCOPE, uncontended_timing=False, publication_ready=False,
        limitations=["Distinct launch-failure disposition, not an ordinary successful terminal resource review.",
            "Allocation elapsed time is not native inference runtime; missing resource/accuracy values remain missing.",
            "Runtime brackets and a fresh inventory are checked, not continuous runtime integrity.",
            "No completed scientific output exists to recover or score.",
            "Prefix resolution does not solve scheduler placement or authorize an incompatible allocation.",
            "Other native/measurement failures are not covered; no retry or relaxed accounting gates."])
    save(destination / "review.json", report)
    return record(destination / "review.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Explicit request binding differs")
    print(review(ref, args.output_directory))
