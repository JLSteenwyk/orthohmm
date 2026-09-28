from copy import deepcopy

import pytest

from benchmark_tools.replay_threadripper_scaling import validate_job_memory


def fixture():
    before = dict(scope="/slurm/job_1", errors=[], started_ns=1, finished_ns=2,
        raw={"memory.current": "100", "memory.peak": "200", "memory.max": str(128*1024**3),
             "memory.swap.max": "max", "memory.stat": "shmem 64\n",
             "memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\n"})
    after = deepcopy(before)
    after.update(started_ns=8, finished_ns=9)
    after["raw"]["memory.peak"] = "300"
    return before, after, "/slurm/job_1", dict(started_ns=3, finished_ns=5), dict(finished_ns=7)


def test_job_peak_not_added_or_subtracted():
    result = validate_job_memory(*fixture())
    assert result["peak_bytes_since_job_creation"] == 300
    assert result["includes_preparation_and_observer"] is True
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", ["scope", "error", "ram", "reset", "window", "events", "shmem", "swap", "gauge"])
def test_reject_wrong_memory_evidence(change):
    args = fixture()
    after = args[1]
    if change == "scope":
        after["scope"] = "/slurm/job_2"
    elif change == "error":
        after["errors"] = [dict(error="read failed")]
    elif change == "ram":
        after["raw"]["memory.max"] = str(2*1024**3)
    elif change == "reset":
        after["raw"]["memory.peak"] = "150"
    elif change == "window":
        after["started_ns"] = 6
    elif change == "events":
        after["raw"]["memory.events"] = "low 0\n"
    elif change == "shmem":
        after["raw"]["memory.stat"] = "anon 100\n"
    elif change == "swap":
        after["raw"]["memory.swap.max"] = "-1"
    else:
        after["raw"]["memory.current"] = "400"
    with pytest.raises(ValueError):
        validate_job_memory(*args)
