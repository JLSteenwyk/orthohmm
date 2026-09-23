import io
import json

import pytest

from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.replay_host_process_observation import replay
from tests.unit.test_command_host_monitor import sample


def fixture(tmp_path, cpu=0., failure=False):
    values = iter([sample(1.), None, sample(3.), sample(4., cpu)] if failure else
                  [sample(1.), sample(4., cpu)])
    def observe():
        value = next(values)
        if value is None:
            raise OSError("fixture")
        return value
    handle = io.StringIO()
    monitor = HostMonitor(handle, "/slurm/job_1", observe)
    for _ in range(4 if failure else 2):
        monitor.observe()
    rows = [json.loads(line) for line in handle.getvalue().splitlines()]
    if failure:
        rows[1]["at_monotonic_s"] = 2.
    path = tmp_path / "observations.jsonl"
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    return path, monitor.summary(1.5, 3.5), rows


@pytest.mark.parametrize("cpu,failure,status", [(0.,False,"no_large_persistent_competitor_observed"),
    (2.,False,"competing_cpu_observed"), (0.,True,"inconclusive")])
def test_recompute_stream_and_retain_uncertainty(tmp_path, cpu, failure, status):
    path, summary, _ = fixture(tmp_path, cpu, failure)
    result = replay(path, summary, "/slurm/job_1", 1.5, 3.5)
    assert result["status"] == status
    assert result["controlled_workload_verified"] is False


@pytest.mark.parametrize("fault", ["index", "observer", "interval", "time", "truncate", "summary", "scope", "boundary"])
def test_tampered_observations_rejected(tmp_path, fault):
    path, summary, rows = fixture(tmp_path, 2.)
    scope, end = "/slurm/job_1", 3.5
    if fault == "index": rows[1]["index"] = 3
    elif fault == "observer": rows[1]["observer_pid"] += 1
    elif fault == "interval": rows[1]["interval"]["sum_observed_foreign_average_cores"] = 0.
    elif fault == "time": rows[1]["snapshot"]["finished_monotonic_s"] = 0.
    elif fault == "truncate": rows.pop()
    elif fault == "summary": summary["controlled_workload_verified"] = True
    elif fault == "scope": scope = "/other"
    elif fault == "boundary": end = 5.
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    with pytest.raises(ValueError):
        replay(path, summary, scope, 1.5, end)
