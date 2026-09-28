"""Replay typed host process samples without treating them as timing admission.

Only two snapshots are retained in memory. The caller must bind the input bytes,
reviewed policy and native interval to the actual run and verify their provenance.
"""

from collections import Counter
import json
from pathlib import Path

from benchmark_tools.review_threadripper_process_policy import number, review
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.slurm_resource_snapshot import scoped_path


def evaluate(lines, policy, *, boot_id, job_scope, observer_pid, launch, end,
             maximum_foreign_average_cores, maximum_sample_period_s):
    """Check every consecutive observation, never just endpoints or mean load.

    Invalid records break continuity and cannot be bridged into a passing stream.
    Thresholds must be supplied prospectively, not selected from observed results.
    """
    number(launch, positive=True)
    number(end, positive=True)
    if end <= launch or policy.get("schema") != "threadripper_process_policy_v2":
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
                if cores > maximum_foreign_average_cores:
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
    return dict(schema="threadripper_process_stream_review_v1",
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


def audit(directory, policy_ref, preflight_ref, *, job_id, index):
    """Bind a completed collector stream to its reviewed policy and native clock."""
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
    preflight = read(preflight_ref)
    if (policy.get("schema") != "threadripper_environment_policy_v1"
            or policy.get("decision") != "reviewed" or policy.get("host") != "bizon"
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
    references.extend([process_ref, *policy["evidence"], *preflight["evidence"]])
    refs = {name: record(directory / name) for name in
            ("ready.json", "done.json", "host_processes.jsonl")}
    references.extend(refs.values())
    for ref in references:
        check(ref)
    ready = read(refs["ready.json"])
    done = read(refs["done.json"])
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
    for ref in references:
        check(ref)
    result.update(job_id=job_id, index=index, evidence=references,
                  source=record(Path(__file__).resolve()))
    save(output, result)
    return record(output), result
