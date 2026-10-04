"""Answer one parked-worker preflight request; never submit, release or retry jobs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import psutil

from benchmark_tools.observe_host_competition import membership
from benchmark_tools.audit_dgx_pressure import pressure_summary
from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.probe_host_counters import snapshot as host_snapshot
from benchmark_tools.review_threadripper_process_policy import number, review as process_review, shared_environment, SHARED_SCOPE
from benchmark_tools.review_threadripper_pressure_stream import native_pressure_role
from benchmark_tools.run_threadripper_scaling import PLAN_SHA, deployment, expect, read, review, select
from benchmark_tools.slurm_resource_snapshot import scoped_path


def job_scope(pid, job):
    scope = scoped_path(Path(f"/proc/{pid}/cgroup").read_text(), job)
    return str(Path(*scope.parts[:scope.parts.index(f"job_{job}") + 1]))


def live_identity(row):
    process = psutil.Process(row["pid"])
    if (process.create_time() != row["created"] or membership(row["pid"]) != row["cgroup"]
            or process.name() != row["name"] or not process.is_running()):
        raise ValueError("Approved user process changed during executable check")


def images(policy, snapshot):
    """Hash loaded main executables, not merely the current files at their paths."""
    observed = {r["pid"]: r for r in snapshot["processes"]}
    cache, results = {}, []
    for entry in policy["ordinary_processes"]:
        if entry["classification"] == "reviewed_kernel_thread":
            continue
        row = observed[entry["pid"]]
        image_ref = entry["image"]
        check(image_ref)
        live_identity(row)
        path = Path(f"/proc/{row['pid']}/exe")
        with path.open("rb") as handle:
            stat = os.fstat(handle.fileno())
            key = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            if key not in cache:
                digest = hashlib.sha256()
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
                cache[key] = digest.hexdigest()
            current = os.fstat(handle.fileno())
            if key != (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns, current.st_ctime_ns):
                raise ValueError("Loaded executable changed during hashing")
        current = path.stat()
        if key != (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns, current.st_ctime_ns):
            raise ValueError("Process switched executable during hashing")
        live_identity(row)
        if cache[key] != image_ref["sha256"] or stat.st_size != image_ref["bytes"]:
            raise ValueError("Loaded executable differs from reviewed image")
        results.append(dict(pid=row["pid"], created=row["created"], image=image_ref,
                            loaded_sha256=cache[key], device=stat.st_dev, inode=stat.st_ino))
    return results


def collector_ready(directory, scope, job, *, shared_host=False):
    ready_path = directory / "ready.json"
    ready_ref = record(ready_path)
    ready = read(ready_ref)
    if job_scope(ready["pid"], job) != scope:
        raise ValueError("Parked native worker belongs to another job")
    if membership(ready["pid"]) != str(scoped_path(ready["cgroup"], job)):
        raise ValueError("Parked native worker changed step membership")
    stream = directory / "host_processes.jsonl"
    with stream.open() as handle:
        first = json.loads(handle.readline())
        if not shared_host and handle.readline():
            raise ValueError("Collector already advanced beyond its initial sample")
    if first.get("index") != 0 or first.get("interval") is not None:
        raise ValueError("Missing initial collector observation")
    pid = first["observer_pid"]
    rows = [r for r in first["snapshot"]["processes"] if r["pid"] == pid]
    errors = first["snapshot"]["errors"]
    if (not isinstance(errors, list) or len(rows) != 1 or job_scope(pid, job) != scope
            or errors and (not shared_host or any(
                error.get("pid") in {pid, ready["pid"]} for error in errors))):
        raise ValueError("Collector identity, initial sample or job scope is invalid")
    live_identity(rows[0])
    native_rows = [r for r in first["snapshot"]["processes"] if r["pid"] == ready["pid"]]
    if len(native_rows) != 1:
        raise ValueError("Native worker missing from initial collector sample")
    live_identity(native_rows[0])
    if (directory / "go.json").exists() or (directory / "done.json").exists():
        raise ValueError("Native worker has already been released or finished")
    if not shared_host:
        return [ready_ref, record(stream)]
    # The periodic observer owns the append-only stream; pin only its first sample.
    initial_path = directory / "preflight_initial_process_sample.json"
    if not initial_path.exists():
        save(initial_path, first)
    initial_ref = record(initial_path)
    if read(initial_ref) != first:
        raise ValueError("Initial collector observation changed during preflight")
    return [ready_ref, initial_ref]


def respond(request_ref, policy_ref, *, root=None, sample=enriched_snapshot,
            image_check=images, collector_check=collector_ready, clock=time.time_ns,
            sleep=time.sleep, waiter=wait_file, host_sample=host_snapshot,
            capacity=psutil.virtual_memory):
    root = Path(root or Path(__file__).resolve().parent.parent)
    request = read(request_ref)
    selected = deployment(request)
    job = int(os.environ["SLURM_JOB_ID"])
    if os.uname().nodename != "bizon":
        raise ValueError("Require the authorized local host")
    run, _, _, sources = select(request, root, job)
    scope = job_scope(os.getpid(), job)
    policy = review(policy_ref, dict(decision="reviewed",
                                    host="bizon", plan_sha256=selected["plan_sha"]))
    pressure_role = native_pressure_role(policy)
    shared = shared_environment(policy)
    if shared != (request.get("execution_scope") == SHARED_SCOPE):
        raise ValueError("Request and environmental execution scopes differ")
    if selected["name"] == "private_v2_20260928" and pressure_role != "diagnostic_only":
        raise ValueError("Private timing requires explicit v2 diagnostic-only native pressure policy")
    ready = read(request["readiness_review"])
    if ready.get("environment_policy") != policy_ref:
        raise ValueError("Readiness does not bind this environmental policy")
    process_policy = read(policy["process_policy"])
    expect(process_policy, dict(schema="threadripper_process_policy_v3" if shared else "threadripper_process_policy_v2"))
    cpu_limit = number(policy["maximum_foreign_average_cores"])
    pressure_limits = policy["maximum_pressure_percent"]
    if set(pressure_limits) != {"cpu", "io", "memory"} or any(
            number(value) > 100 for value in pressure_limits.values()):
        raise ValueError("Require prospectively reviewed CPU, I/O and memory pressure bounds")
    configuration = policy["configuration_files"]
    if not isinstance(configuration, list) or not configuration and not shared:
        raise ValueError("Require reviewed service/configuration file evidence")
    minimum_available = policy.get("minimum_available_memory_bytes") if shared else None
    if shared and (type(minimum_available) is not int or minimum_available < 128 * 1024**3):
        raise ValueError("Shared-host launch must allow the full 128-GiB native memory limit")
    evidence = [request_ref, policy_ref, request["readiness_review"], request["recipe"],
                policy["process_policy"], *configuration, *policy["evidence"], *sources]
    for ref in evidence:
        check(ref)
    directory = Path(run["measurement_directory"])
    marker_path = directory / "environment_review_requested.json"
    # Preparation may take up to an hour. This wait does not authorize a rerun.
    waiter(marker_path, seconds=4000)
    marker_ref = record(marker_path)
    marker = read(marker_ref)
    expect(marker, dict(job_id=job, index=request["index"], request=request_ref,
                       review_path=request["environment_preflight_path"], wait_seconds=20,
                       native_released=False))
    output = Path(request["environment_preflight_path"])
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    detail_path = output.with_name("environment_worker_evidence.json")
    started = clock()
    requested = marker["requested_unix_ns"]
    if type(requested) is not int or not 0 < requested <= started < requested + 20_000_000_000:
        raise ValueError("Invalid or expired environmental request")
    detail = dict(schema="threadripper_environment_worker_evidence_v1", status="failed",
                  request=request_ref, policy=policy_ref, marker=marker_ref, snapshots=[],
                  scientific_timings_admitted=False, native_release_authorized=False,
                  limitations=["Policy classifications and configuration-inventory completeness require independent review.",
                      "Main-executable hashes do not attest shared libraries, interpreted code or future process behavior.",
                      "Preflight observations cannot exclude short-lived work or certify whole-run isolation.",
                      "A preflight pass does not submit a job, release native inference or admit timing results."])
    response = dict(schema="threadripper_environment_preflight_v1", decision="failed", job_id=job,
        index=request["index"], plan_sha256=selected["plan_sha"], recipe_sha256=request["recipe"]["sha256"],
        readiness_review_sha256=request["readiness_review"]["sha256"], environment_policy=policy_ref,
        whole_run_observer_ready=False, unrelated_scientific_work_present=None,
        observation_started_unix_ns=started, review_reference=policy["review_reference"],
        scientific_timings_admitted=False)
    if pressure_role == "diagnostic_only":
        response.update(native_pressure_role=pressure_role, preflight_pressure_limits_used=not shared)
    if shared:
        response.update(execution_scope=SHARED_SCOPE, background_competition_recorded=False,
                        uncontended_timing=False, foreign_cpu_used_for_eligibility=False)
    def check_collector():
        extra = dict(shared_host=True) if shared and collector_check is collector_ready else {}
        return collector_check(directory, scope, job, **extra)
    try:
        collector_refs = check_collector()
        if shared:
            detail["available_memory_bytes"] = [capacity().available]
        detail["host_snapshots"] = [host_sample()]
        before = sample()
        detail["snapshots"].append(before)
        sleep(3.)
        after = sample()
        detail["snapshots"].append(after)
        boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        result = process_review(process_policy, before, after, boot_id=boot,
                                job_scope=scope, observer_pid=os.getpid())
        detail["process_review"] = result
        response["boot_id"] = boot
        if not result["process_policy_matched"]:
            raise ValueError("Unresolved outside process inventory")
        if not shared and result["cpu_diagnostic"]["sum_observed_foreign_average_cores"] > cpu_limit:
            raise ValueError("Reviewed prospective background CPU bound exceeded")
        detail["loaded_images"] = [] if shared else image_check(process_policy, after)
        detail["host_snapshots"].append(host_sample())
        if any(s["errors"] for s in detail["host_snapshots"]):
            raise ValueError("Host pressure/counter observation errors")
        if (any(s["raw"]["boot_id"].strip() != boot for s in detail["host_snapshots"])
                or job_scope(os.getpid(), job) != scope):
            raise ValueError("Host observations or review worker changed identity")
        detail["pressure"] = {resource: pressure_summary(detail["host_snapshots"], resource)
                              for resource in pressure_limits}
        if not shared and any(detail["pressure"][resource]["midpoint_percent"]["some"] > limit
               for resource, limit in pressure_limits.items()):
            raise ValueError("Reviewed prospective host pressure bound exceeded")
        if shared:
            detail["available_memory_bytes"].append(capacity().available)
            if any(type(value) is not int or value < minimum_available
                   for value in detail["available_memory_bytes"]):
                raise ValueError("Insufficient safe available memory for the shared-host launch")
        for ref in [*evidence, marker_ref, *collector_refs]:
            check(ref)
        check_collector()
        if Path("/proc/sys/kernel/random/boot_id").read_text().strip() != boot:
            raise ValueError("Boot changed during environmental review")
        if not requested <= clock() < requested + 20_000_000_000:
            raise TimeoutError("Environmental review missed its release deadline")
        response.update(decision="passed", whole_run_observer_ready=True)
        if shared:
            response.update(background_competition_recorded=True,
                observed_foreign_average_cores=result["cpu_diagnostic"]["sum_observed_foreign_average_cores"],
                available_memory_bytes=detail["available_memory_bytes"])
            detail["limitations"].append("Shared-host observation: background images/services are not certified ordinary or isolated.")
        else:
            response["unrelated_scientific_work_present"] = False
        detail.update(status="reviewed_preflight_observation", collector_evidence=collector_refs)
    except BaseException as error:
        response["decision"] = "failed"
        detail.update(error_type=type(error).__name__, error=str(error))
        if not isinstance(error, Exception):
            raise
    finally:
        response["observation_finished_unix_ns"] = clock()
        if not started <= response["observation_finished_unix_ns"] < requested + 20_000_000_000:
            response["decision"] = "failed"
            detail["status"] = "deadline_expired"
        save(detail_path, detail)
        response["evidence"] = [record(detail_path), *evidence, marker_ref]
        response["observation_finished_unix_ns"] = clock()
        if not started <= response["observation_finished_unix_ns"] < requested + 20_000_000_000:
            response.update(decision="failed", publication_deadline_expired=True)
        save(output, response)
    return response


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--policy-sha256", required=True)
    args = parser.parse_args()
    request_ref, policy_ref = record(args.request), record(args.policy)
    if request_ref["sha256"] != args.request_sha256 or policy_ref["sha256"] != args.policy_sha256:
        raise ValueError("Request or policy digest differs")
    raise SystemExit(0 if respond(request_ref, policy_ref)["decision"] == "passed" else 1)
