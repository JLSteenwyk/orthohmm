import pytest

from benchmark_tools.probe_tmpfs_memory_charge import evaluate, SIZE


def observation(scope, start, shmem):
    return dict(scope=scope, started_ns=start, finished_ns=start+1, raw={
        "memory.current": str(shmem+100), "memory.peak": str(shmem+200),
        "memory.stat": f"shmem {shmem}\n"})


def fixture():
    before = observation("/cgroup/job_3", 1, 0)
    prepared = observation("/cgroup/job_3", 3, SIZE)
    native = dict(before=observation("/cgroup/job_3/step_0", 5, 0),
                  after=observation("/cgroup/job_3/step_0", 7, 0), bytes_read=SIZE, sha256="expected")
    return before, prepared, native


def test_charge_observed():
    result = evaluate(*fixture(), "expected")
    assert result["prepared_job_shmem_delta_bytes"] == SIZE
    assert result["native_step_shmem_after_bytes"] == 0
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", ["job", "step", "bytes", "hash", "charge", "native_charge", "window"])
def test_reject_invalid_control(change):
    before, prepared, native = fixture()
    if change == "job":
        prepared["scope"] = "/cgroup/job_4"
    elif change == "step":
        native["before"]["scope"] = "/cgroup/job_3/step_batch"
    elif change == "bytes":
        native["bytes_read"] -= 1
    elif change == "hash":
        native["sha256"] = "wrong"
    elif change == "charge":
        prepared["raw"]["memory.stat"] = "shmem 0\n"
    elif change == "native_charge":
        native["after"]["raw"]["memory.stat"] = f"shmem {SIZE}\n"
    else:
        native["before"]["started_ns"] = 2
    with pytest.raises(ValueError):
        evaluate(before, prepared, native, "expected")
