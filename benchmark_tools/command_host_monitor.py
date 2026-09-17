"""Streaming host CPU observations bracketing a measured command, without exclusivity claims."""

from collections import Counter
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observe_host_competition import analyze, snapshot
from monitor_slurm_resources import source_record


class HostMonitor:
    def __init__(self, handle, scope, sample_fn=None):
        self.handle = handle
        self.scope = scope
        self.observer_pid = os.getpid()
        self.sample_fn = sample_fn or snapshot
        self.previous = None
        self.first_finished = None
        self.last_started = None
        self.count = 0
        self.errors = 0
        self.intervals = Counter()
        self.maximum_cores = 0.
        self.maximum_gap = 0.

    def observe(self):
        try:
            sample = self.sample_fn()
            interval = None
            if self.previous is not None:
                interval = analyze(self.previous, sample, self.scope, self.observer_pid)
                self.intervals[interval["status"]] += 1
                self.maximum_cores = max(self.maximum_cores, interval["sum_observed_foreign_average_cores"])
                self.maximum_gap = max(self.maximum_gap, sample["started_monotonic_s"] - self.previous["finished_monotonic_s"])
            self.handle.write(json.dumps({"index": self.count, "observer_pid": self.observer_pid,
                                          "snapshot": sample, "interval": interval}, sort_keys=True) + "\n")
            self.handle.flush()
            if self.first_finished is None:
                self.first_finished = sample["finished_monotonic_s"]
            self.last_started = sample["started_monotonic_s"]
            self.previous = sample
            self.count += 1
        except Exception as error:
            # Native command execution continues; missing workload evidence is never quiet evidence.
            self.errors += 1
            self.previous = None
            try:
                self.handle.write(json.dumps({"observation_error": type(error).__name__, "at_monotonic_s": time.monotonic()}) + "\n")
                self.handle.flush()
            except Exception:
                pass

    def summary(self, launch, end):
        bracketed = (self.count >= 2 and self.first_finished <= launch and self.last_started >= end)
        state = ("competing_cpu_observed" if self.intervals["competing_cpu_observed"] else
                 "inconclusive" if self.errors or self.intervals["inconclusive"] or not bracketed or not self.intervals else
                 "no_large_persistent_competitor_observed")
        return {"status": state, "controlled_workload_verified": False, "command_bracketed_by_samples": bracketed,
                "observer_pid": self.observer_pid, "command_launch_started_monotonic_s": launch,
                "command_wait_finished_monotonic_s": end,
                "successful_snapshots": self.count, "observation_errors": self.errors,
                "interval_counts": dict(self.intervals), "maximum_observed_foreign_average_cores": self.maximum_cores,
                "maximum_between_snapshot_gap_s": self.maximum_gap, "threshold_average_cores": .25,
                "source": source_record(__file__), "observer_source": source_record(Path(__file__).with_name("observe_host_competition.py")),
                "limitations": ["Snapshots miss short-lived processes and non-CPU contention.",
                                "Per-process CPU intervals are not simultaneous; no exclusivity or quiet-host certification.",
                                "Sampling overhead is included in measured resources; raw process inventory is local evidence."]}
