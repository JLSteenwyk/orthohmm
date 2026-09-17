from pathlib import PurePosixPath

import pytest

from benchmark_tools.slurm_resource_snapshot import scoped_path, counters, cpus, metrics


def test_scope_requires_requested_job_and_step():
    text = "0::/system.slice/slurmstepd.scope/job_123/step_batch/user/task_0\n"
    assert scoped_path(text, 123) == PurePosixPath(text.strip().split("::")[1])
    with pytest.raises(ValueError):
        scoped_path(text, 12)


@pytest.mark.parametrize("text", ["0::/user.slice/other", "0::/job_123", "0::/job_123/notstep", "0::/job_123/step_batch/../x",
                                  "0::relative/job_123/step_batch", "0::/job_123/step_batch\n0::/job_123/step_batch"])
def test_reject_unscoped_or_ambiguous_path(text):
    with pytest.raises(ValueError):
        scoped_path(text, 123)


@pytest.mark.parametrize("text", ["a 1\na 2", "a -1", "a x", "a 1 extra"])
def test_reject_invalid_counters(text):
    with pytest.raises(ValueError):
        counters(text)


def test_counter_order_and_cpuset():
    assert counters("system_usec 2\nusage_usec 5\nuser_usec 3\n") == {"usage_usec": 5, "user_usec": 3, "system_usec": 2}
    assert cpus("1-3,8,10-11\n") == [1, 2, 3, 8, 10, 11]


@pytest.mark.parametrize("text", ["", "x", "4-2", "1-2-3", "1,", "0-999999"])
def test_invalid_cpuset(text):
    with pytest.raises(ValueError):
        cpus(text)


def test_metrics_preserve_raw_and_distinct_accounting(tmp_path):
    raw = {"cpu.stat": "usage_usec 5\nuser_usec 3\nsystem_usec 2\n", "memory.current": "100\n", "memory.peak": "200\n",
           "memory.stat": "anon 60\nfile 40\n", "memory.events": "oom 0\noom_kill 0\n",
           "cpuset.cpus.effective": "2-3\n", "memory.max": "max\n"}
    for name, value in raw.items():
        (tmp_path / name).write_text(value)
    result = metrics(tmp_path)
    assert result["raw"] == raw
    assert result["memory_current_bytes"] == 100
    assert result["memory_peak_since_creation_or_reset_bytes"] == 200
    assert result["memory_max_bytes"] is None
    assert result["effective_cpus"] == [2, 3]
