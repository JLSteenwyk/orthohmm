import json
import os
import sys
import time

import pytest

from benchmark_tools.audit_slurm_measurement import audit, source_record
from benchmark_tools.measure_slurm_command import measure


@pytest.fixture
def evidence(tmp_path):
    def snapshot(pid, job):
        return {"pid": pid, "job_id": job, "scope": "/job_123/step_batch/user/task_0",
                "processes": [{"pid": pid, "rss_bytes": 100}], "sampling_errors": [],
                "sampled_sum_process_rss_bytes": 100,
                "metrics": {"effective_cpus": [0], "cpu": {"usage_usec": 1, "user_usec": 1, "system_usec": 0},
                            "memory_peak_since_creation_or_reset_bytes": 100},
                "ancestor_limits": [{"memory_max": "1024"}]}

    def host():
        now = time.monotonic()
        return {"started_monotonic_s": now, "finished_monotonic_s": now, "processes": [], "errors": []}

    directory = tmp_path / "measurement"
    measure([sys.executable, "-c", "pass"], directory, 123, 1, 1024, 5.,
            snapshot_fn=snapshot, monitor_host=True, host_snapshot_fn=host)
    return directory


def test_saved_evidence_replays_without_quiet_host_claim(evidence):
    result = audit(evidence, source_record(evidence / "results.json")["sha256"])
    assert result["status"] == "collector_evidence_replayed"
    assert result["benchmark_admitted"] is False
    assert result["controlled_workload_verified"] is False
    assert result["resource_observations"] == result["host_observations"] == 2


def test_wrong_report_hash(evidence):
    with pytest.raises(ValueError, match="report hash"):
        audit(evidence, "wrong")


@pytest.mark.parametrize("change", ["raw", "duration", "summary", "cpuset", "cap", "elapsed", "pid", "interval"])
def test_changed_evidence_rejected_even_if_report_rehashed(evidence, change):
    path = evidence / "results.json"
    report = json.loads(path.read_text())
    if change == "raw":
        (evidence / "command.log").write_text("changed")
    elif change == "duration":
        report["command_wall_s"] += 1
    elif change == "summary":
        report["summary"]["maximum_sampled_sum_rss_bytes"] += 1
    else:
        name = "host_samples.jsonl" if change in {"pid", "interval"} else "samples.jsonl"
        rows = [json.loads(line) for line in (evidence / name).read_text().splitlines()]
        if change == "cpuset":
            rows[-1]["snapshot"]["metrics"]["effective_cpus"] = [1]
        elif change == "cap":
            rows[-1]["snapshot"]["ancestor_limits"] = [{"memory_max": "2048"}]
        elif change == "elapsed":
            rows[-1]["elapsed_s"] += 1
        elif change == "pid":
            rows[-1]["observer_pid"] = os.getpid() + 1
        else:
            rows[-1]["interval"]["sum_observed_foreign_average_cores"] = 999
        (evidence / name).write_text("".join(json.dumps(row) + "\n" for row in rows))
        report[name] = source_record(evidence / name)
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        audit(evidence, source_record(path)["sha256"])
