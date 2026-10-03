"""Replay typed host process samples without treating them as timing admission.

Only two snapshots are retained in memory. The caller must bind the input bytes,
reviewed policy and native interval to the actual run and verify their provenance.
"""

from collections import Counter
import json
from pathlib import Path

from benchmark_tools.review_threadripper_process_policy import number, review, SHARED_SCOPE, shared_environment
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.review_threadripper_pressure_stream import evaluate as evaluate_pressure
from benchmark_tools.review_threadripper_pressure_stream import native_pressure_role
from benchmark_tools.verify_lineage_native_provenance import same


def evaluate(lines, policy, *, boot_id, job_scope, observer_pid, launch, end,
             maximum_foreign_average_cores, maximum_sample_period_s):
    """Check every consecutive observation, never just endpoints or mean load.

    Invalid records break continuity and cannot be bridged into a passing stream.
    Thresholds must be supplied prospectively, not selected from observed results.
    """
    number(launch, positive=True)
    number(end, positive=True)
    shared = policy.get("schema") == "threadripper_process_policy_v3"
    if shared and policy.get("execution_scope") != SHARED_SCOPE:
        raise ValueError("Require explicit shared-host scope")
    if end <= launch or policy.get("schema") not in {"threadripper_process_policy_v2", "threadripper_process_policy_v3"}:
        raise ValueError("Require a positive native interval and typed process policy")
    number(maximum_foreign_average_cores)
    number(maximum_sample_period_s, positive=True)
    failures = Counter()
    previous = None
    first_finished = last_started = None
    expected_index = records = intervals = matched = 0
    maximum_cores = maximum_period = 0.
    for line in lines:
        records += 1
        try:
            row = json.loads(line)
            if not isinstance(row, dict) or "observation_error" in row:
                raise ValueError("Missing process observation")
            if (type(row.get("index")) is not int or row["index"] != expected_index
                    or type(row.get("observer_pid")) is not int
                    or row["observer_pid"] != observer_pid):
                raise ValueError("Process stream index or observer differs")
            expected_index += 1
            sample = row["snapshot"]
            start = number(sample["started_monotonic_s"])
            finish = number(sample["finished_monotonic_s"])
            if finish < start:
                raise ValueError("Invalid snapshot duration")
            if first_finished is None:
                first_finished = finish
            last_started = start
            if previous is not None:
                verdict = review(policy, previous, sample, boot_id=boot_id,
                                 job_scope=job_scope, observer_pid=observer_pid)
                intervals += 1
                period = start - previous["started_monotonic_s"]
                maximum_period = max(maximum_period, period)
                cores = verdict["cpu_diagnostic"]["sum_observed_foreign_average_cores"]
                maximum_cores = max(maximum_cores, cores)
                if not verdict["process_policy_matched"]:
                    failures["process_policy_mismatch"] += 1
                if not shared and cores > maximum_foreign_average_cores:
                    failures["foreign_cpu_bound_exceeded"] += 1
                if period > maximum_sample_period_s or finish - start > maximum_sample_period_s:
                    failures["sample_period_bound_exceeded"] += 1
                if verdict["process_policy_matched"]:
                    matched += 1
            elif finish - start > maximum_sample_period_s:
                failures["sample_period_bound_exceeded"] += 1
            previous = sample
        except (ValueError, TypeError, KeyError, OverflowError) as error:
            failures["invalid_record_" + type(error).__name__] += 1
            previous = None
    bracketed = (first_finished is not None and last_started is not None
                 and first_finished <= launch and last_started >= end)
    if not bracketed:
        failures["native_interval_not_bracketed"] += 1
    if intervals == 0 or intervals != records - 1:
        failures["incomplete_interval_chain"] += 1
    passed = not failures
    result = dict(schema="threadripper_process_stream_review_v1",
        status="sampled_process_policy_satisfied" if passed else "unresolved_process_stream",
        sampled_process_policy_satisfied=passed, scientific_timings_admitted=False,
        controlled_workload_verified=False, records=records, intervals=intervals,
        policy_matched_intervals=matched, failures=dict(failures),
        boot_id=boot_id, job_scope=job_scope, observer_pid=observer_pid,
        command_bracketed_by_samples=bracketed, launch=launch, end=end,
        maximum_observed_foreign_average_cores=maximum_cores,
        maximum_observed_sample_period_s=maximum_period,
        bounds=dict(maximum_foreign_average_cores=maximum_foreign_average_cores,
                    maximum_sample_period_s=maximum_sample_period_s),
        limitations=["Input and policy provenance must be independently bound to the run.",
            "Periodic snapshots miss short-lived work and cannot certify exclusivity.",
            "This checks process identities and CPU only, not executable or configuration drift, PSI or device contention.",
            "All intervals are checked, including bracketing intervals outside native execution.",
            "No production timing admission, automatic retry or next-run authorization."])
    if shared:
        result.update(execution_scope=SHARED_SCOPE, foreign_cpu_used_for_eligibility=False)
        result["limitations"].append("Background activity and churn do not reject shared-host observations; distortion is unknown.")
    return result


def audit(directory, policy_ref, preflight_ref, *, job_id, index, collector_arm="periodic"):
    """Bind a completed collector stream to its reviewed policy and native clock."""
    if not isinstance(collector_arm, str) or collector_arm not in {"periodic", "boundary"}:
        raise ValueError("Require explicit supported native collector arm")
    directory = Path(directory).resolve()
    output = directory / "process_stream_review.json"
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    references = [policy_ref, preflight_ref]

    def read(ref):
        path = Path(ref["path"])
        if not path.is_absolute() or path.is_symlink() or path.resolve() != path:
            raise ValueError("Require direct absolute evidence path")
        check(ref)
        value = json.loads(path.read_text())
        check(ref)
        return value

    policy = read(policy_ref)
    shared = shared_environment(policy)
    pressure_role = native_pressure_role(policy)
    if collector_arm == "boundary" and pressure_role != "diagnostic_only":
        raise ValueError("Boundary environment review requires diagnostic-only native pressure")
    preflight = read(preflight_ref)
    if (policy.get("decision") != "reviewed" or policy.get("host") != "bizon"
            or preflight.get("schema") != "threadripper_environment_preflight_v1"
            or preflight.get("decision") != "passed"
            or type(job_id) is not int or job_id <= 0 or type(index) is not int or index < 0
            or type(preflight.get("job_id")) is not int or preflight["job_id"] != job_id
            or type(preflight.get("index")) is not int or preflight["index"] != index
            or preflight.get("environment_policy") != policy_ref
            or not policy.get("plan_sha256") or preflight.get("plan_sha256") != policy["plan_sha256"]):
        raise ValueError("Policy and preflight do not identify the same reviewed attempt")
    process_ref = policy["process_policy"]
    process_policy = read(process_ref)
    if shared != (process_policy.get("schema") == "threadripper_process_policy_v3"):
        raise ValueError("Environment and process contention scopes differ")
    if shared and (preflight.get("execution_scope") != SHARED_SCOPE
                   or preflight.get("background_competition_recorded") is not True):
        raise ValueError("Shared-host preflight scope or annotation differs")
    configuration = policy.get("configuration_files")
    if not isinstance(configuration, list) or not configuration and not shared:
        raise ValueError("Require reviewed configuration file inventory")
    references.extend([process_ref, *configuration, *policy["evidence"], *preflight["evidence"]])
    refs = {name: record(directory / name) for name in
            ("ready.json", "done.json", "host_processes.jsonl")}
    references.extend(refs.values())
    for ref in references:
        check(ref)
    ready = read(refs["ready.json"])
    done = read(refs["done.json"])
    if collector_arm == "boundary":
        boundary_ref = record(directory / "boundary_report.json")
        boundary = read(boundary_ref)
        if (boundary.get("schema") != "threadripper_boundary_control_v1"
                or boundary.get("collector_arm") != "boundary"
                or type(boundary.get("job_id")) is not int or boundary["job_id"] != job_id
                or not same(boundary.get("native"), done)
                or not same(boundary.get("policy"), dict(native_points=2, periodic_native_sampling=False,
                    completion_poll_interval_s=1., common_host_interval_s=30.))):
            raise ValueError("Boundary report does not identify the expected native collector")
        references.append(boundary_ref)
    for key in ("started_ns", "finished_ns"):
        if type(done.get(key)) is not int or done[key] <= 0:
            raise ValueError("Require integer native monotonic timestamps")
    scope = scoped_path(ready["cgroup"], job_id)
    job_scope = next(p for p in scope.parents if p.name == f"job_{job_id}")
    with Path(refs["host_processes.jsonl"]["path"]).open() as handle:
        first = json.loads(handle.readline())
        observer = first["observer_pid"]
        handle.seek(0)
        result = evaluate(handle, process_policy, boot_id=preflight["boot_id"],
            job_scope=str(job_scope), observer_pid=observer,
            launch=done["started_ns"] / 1e9, end=done["finished_ns"] / 1e9,
            maximum_foreign_average_cores=policy["maximum_foreign_average_cores"],
            maximum_sample_period_s=policy["maximum_sample_period_s"])
    point_paths = sorted(directory.glob("point_*.json"))
    if point_paths != [directory / f"point_{i:06d}.json" for i in range(len(point_paths))]:
        raise ValueError("Require contiguous numbered pressure observation files")
    if collector_arm == "boundary" and len(point_paths) != 2:
        raise ValueError("Boundary review requires exactly two native points")
    def points():
        for path in point_paths:
            ref = record(path)
            references.append(ref)
            yield read(ref)
    extra = dict(boundary_only=True) if collector_arm == "boundary" else {}
    pressure = evaluate_pressure(points(), boot_id=preflight["boot_id"], job_scope=str(job_scope),
        launch_ns=done["started_ns"], end_ns=done["finished_ns"],
        limits=policy["maximum_pressure_percent"],
        maximum_period_s=policy["maximum_pressure_sample_period_s"], pressure_role=pressure_role, **extra)
    if sorted(directory.glob("point_*.json")) != point_paths:
        raise ValueError("Pressure observation inventory changed during review")
    result["pressure_review"] = pressure
    pressure_key = ("sampled_pressure_evidence_satisfied" if pressure_role == "diagnostic_only"
                    else "sampled_pressure_policy_satisfied")
    result["sampled_environment_policy_satisfied"] = bool(result["sampled_process_policy_satisfied"]
        and pressure[pressure_key])
    if pressure_role == "diagnostic_only":
        result.update(schema="threadripper_process_stream_review_v2", native_pressure_role=pressure_role,
                      pressure_thresholds_used_for_eligibility=False)
    if shared:
        result.update(execution_scope=SHARED_SCOPE, uncontended_timing=False,
                      background_cpu_used_for_eligibility=False)
    if collector_arm == "boundary":
        result.update(schema="threadripper_boundary_environment_review_v1", collector_arm=collector_arm,
                      native_pressure_observation="boundary_only", periodic_pressure_cadence_checked=False)
    for ref in references:
        check(ref)
    result["configuration_endpoint_hashes_verified"] = bool(configuration)
    if shared:
        result["background_configuration_review_required"] = False
    result["limitations"].append(
        "Configuration bytes are checked before and after post-run review; transient changes during native execution may be missed.")
    result.update(job_id=job_id, index=index, evidence=references,
                  source=record(Path(__file__).resolve()))
    save(output, result)
    return record(output), result
