"""Recompute host CPU observations; this does not certify an isolated host."""

from collections import Counter
import json
import math
from pathlib import Path

from benchmark_tools.observe_host_competition import analyze
from benchmark_tools.verify_lineage_native_provenance import same


def replay(path, summary, scope, launch, end):
    observer = summary["observer_pid"]
    if type(observer) is not int or observer <= 0:
        raise ValueError("Invalid observer identity")
    previous = None
    first = last = None
    boundary = -math.inf
    count = errors = 0
    intervals = Counter()
    maximum = gap = 0.
    with Path(path).open() as handle:
        for line in handle:
            row = json.loads(line)
            if "observation_error" in row:
                if (set(row) != {"observation_error", "at_monotonic_s"}
                        or not isinstance(row["observation_error"], str) or not row["observation_error"]):
                    raise ValueError("Invalid observation error")
                start = finish = row["at_monotonic_s"]
                errors += 1
                previous = None
            else:
                if (set(row) != {"index", "observer_pid", "snapshot", "interval"}
                        or type(row["index"]) is not int or row["index"] != count
                        or type(row["observer_pid"]) is not int or row["observer_pid"] != observer):
                    raise ValueError("Process stream identity/order differs")
                sample = row["snapshot"]
                start, finish = sample["started_monotonic_s"], sample["finished_monotonic_s"]
                interval = None if previous is None else analyze(previous, sample, scope, observer)
                if not same(interval, row["interval"]):
                    raise ValueError("Process interval does not reproduce")
                if interval is not None:
                    intervals[interval["status"]] += 1
                    maximum = max(maximum, interval["sum_observed_foreign_average_cores"])
                    gap = max(gap, start - previous["finished_monotonic_s"])
                if first is None:
                    first = finish
                last = start
                previous = sample
                count += 1
            if (any(type(t) not in (int, float) or not math.isfinite(t) for t in (start, finish))
                    or not 0 <= start <= finish or start < boundary):
                raise ValueError("Invalid process observation time/order")
            boundary = finish
    bracketed = count >= 2 and first <= launch and last >= end
    status = ("competing_cpu_observed" if intervals["competing_cpu_observed"] else
              "inconclusive" if errors or intervals["inconclusive"] or not bracketed or not intervals else
              "no_large_persistent_competitor_observed")
    expected = dict(status=status, controlled_workload_verified=False,
        command_bracketed_by_samples=bracketed, observer_pid=observer,
        command_launch_started_monotonic_s=launch, command_wait_finished_monotonic_s=end,
        successful_snapshots=count, observation_errors=errors, interval_counts=dict(intervals),
        maximum_observed_foreign_average_cores=maximum, maximum_between_snapshot_gap_s=gap,
        threshold_average_cores=.25)
    if not same({k: summary[k] for k in expected}, expected):
        raise ValueError("Host process summary does not reproduce")
    return expected
