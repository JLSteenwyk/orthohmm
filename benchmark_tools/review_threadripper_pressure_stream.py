"""Check retained per-point PSI evidence without attributing stalls to competitors."""

from collections import Counter

from benchmark_tools.audit_dgx_pressure import pressure_summary
from benchmark_tools.probe_host_counters import parse_group
from benchmark_tools.review_threadripper_process_policy import group, number


def native_pressure_role(policy):
    schema = policy.get("schema")
    if schema == "threadripper_environment_policy_v1" and "native_pressure_role" not in policy:
        return "eligibility"
    if (schema == "threadripper_environment_policy_v2"
            and policy.get("native_pressure_role") == "diagnostic_only"):
        return "diagnostic_only"
    raise ValueError("Require an explicit supported environmental pressure policy")


def evaluate(points, *, boot_id, job_scope, launch_ns, end_ns, limits, maximum_period_s,
             pressure_role="eligibility"):
    if (set(limits) != {"cpu", "memory", "io"}
            or any(number(v) > 100 for v in limits.values())
            or type(launch_ns) is not int or type(end_ns) is not int
            or not 0 < launch_ns < end_ns):
        raise ValueError("Require native timestamps and explicit three-resource PSI limits")
    if not isinstance(pressure_role, str) or pressure_role not in {"eligibility", "diagnostic_only"}:
        raise ValueError("Invalid native pressure role")
    number(maximum_period_s, positive=True)
    scope = group(job_scope)
    failures = Counter()
    exceedances = Counter()
    maxima = dict.fromkeys(limits, 0.)
    previous = None
    first_finish = last_start = None
    count = intervals = 0
    maximum_period = 0.
    for point in points:
        count += 1
        try:
            samples = point["host"]
            if not isinstance(samples, list) or len(samples) != 2:
                raise ValueError("Require two host brackets per point")
            for sample in samples:
                if (sample["errors"] or sample["raw"]["boot_id"].strip() != boot_id
                        or not group(parse_group(sample["raw"]["cgroup_membership"])).is_relative_to(scope)):
                    raise ValueError("Host read errors, boot or job membership differs")
            for resource in limits:
                pressure_summary(samples, resource)
            if first_finish is None:
                first_finish = samples[-1]["finished_monotonic_ns"]
            last_start = samples[-1]["started_monotonic_ns"]
            duration = (samples[-1]["finished_monotonic_ns"] - samples[0]["started_monotonic_ns"]) / 1e9
            if duration > maximum_period_s:
                failures["point_duration_exceeded"] += 1
            if previous is not None:
                # Validate continuity through both brackets, but use cadence-scale
                # intervals for policy rates, not very short within-point spans.
                for resource in limits:
                    pressure_summary([previous, *samples], resource)
                    value = pressure_summary([previous, samples[-1]], resource)
                    percent = value["midpoint_percent"]["some"]
                    maxima[resource] = max(maxima[resource], percent)
                    if percent > limits[resource]:
                        exceedances[resource + "_pressure_bound_exceeded"] += 1
                        if pressure_role == "eligibility":
                            failures[resource + "_pressure_bound_exceeded"] += 1
                intervals += 1
                period = (last_start - previous["started_monotonic_ns"]) / 1e9
                maximum_period = max(maximum_period, period)
                if period > maximum_period_s:
                    failures["pressure_sample_period_exceeded"] += 1
            previous = samples[-1]
        except (ValueError, TypeError, KeyError, OverflowError) as error:
            failures["invalid_point_" + type(error).__name__] += 1
            previous = None
    bracketed = first_finish is not None and first_finish <= launch_ns and last_start >= end_ns
    if not bracketed:
        failures["native_interval_not_bracketed"] += 1
    if not intervals or intervals != count - 1:
        failures["incomplete_pressure_chain"] += 1
    result = dict(schema="threadripper_pressure_stream_review_v1",
        sampled_pressure_policy_satisfied=not failures, points=count, intervals=intervals,
        failures=dict(failures), maximum_observed_some_percent=maxima,
        maximum_observed_period_s=maximum_period, bounds=dict(some_percent=limits,
            maximum_period_s=maximum_period_s), native_interval_bracketed=bracketed,
        scientific_timings_admitted=False, controlled_workload_verified=False,
        limitations=["System PSI includes native and foreign work; a failed bound does not identify its cause.",
            "Non-atomic reads and delayed accounting make midpoint percentages approximate; values are not clipped.",
            "CPU full is undefined at system scope and is not interpreted.",
            "Checks do not establish absence of cache, memory-bandwidth or device interference."])
    if pressure_role == "diagnostic_only":
        # Native-induced stalls are method behavior, not evidence of outside work.
        result.pop("sampled_pressure_policy_satisfied")
        result.update(schema="threadripper_pressure_stream_review_v2",
            native_pressure_role=pressure_role, pressure_thresholds_used_for_eligibility=False,
            sampled_pressure_evidence_satisfied=not failures,
            diagnostic_thresholds_satisfied=not failures and not exceedances,
            diagnostic_threshold_exceedances=dict(exceedances))
        result["limitations"].append(
            "Valid native-interval PSI magnitudes are diagnostic only; they neither exclude a slow method nor certify an uncontended host.")
    return result
