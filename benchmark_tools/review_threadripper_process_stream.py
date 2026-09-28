"""Replay typed host process samples without treating them as timing admission.

Only two snapshots are retained in memory. The caller must bind the input bytes,
reviewed policy and native interval to the actual run and verify their provenance.
"""

from collections import Counter
import json

from benchmark_tools.review_threadripper_process_policy import number, review


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
