import copy
import json

import pytest

from benchmark_tools.monitor_slurm_resources import collect, summarize


def row(index):
    return {"elapsed_s": float(index), "anchor_created": 1., "snapshot": {"job_id": 123, "pid": 456, "scope": "/job_123/step_batch",
        "metrics": {"cpu": {"usage_usec": 1000000 * index, "user_usec": 900000 * index, "system_usec": 100000 * index},
                    "memory_peak_since_creation_or_reset_bytes": 1000}, "sampled_sum_process_rss_bytes": 100 + index,
        "sampling_errors": []}}


def test_summary_deltas_not_lifetime_peak_subtraction():
    result = summarize([row(5), row(6), row(7)])
    assert result["cpu_delta_usec"]["usage_usec"] == 2000000
    assert result["mean_cpu_cores_over_observed_span"] == 1.
    assert result["maximum_reported_cgroup_peak_bytes"] == 1000
    assert result["maximum_sampled_sum_rss_bytes"] == 107
    assert summarize([row(5)])["mean_cpu_cores_over_observed_span"] is None


@pytest.mark.parametrize("problem", ["scope", "pid", "job_id", "identity", "counter", "time"])
def test_reject_incomparable_observations(problem):
    samples = [row(5), row(6)]
    if problem in ("scope", "pid", "job_id"):
        samples[1]["snapshot"][problem] = "different"
    elif problem == "identity":
        samples[1]["anchor_created"] = 2.
    elif problem == "counter":
        samples[1]["snapshot"]["metrics"]["cpu"]["usage_usec"] = 0
    else:
        samples[1]["elapsed_s"] = 5.
    with pytest.raises(ValueError):
        summarize(samples)


def test_success_and_partial_failure_preserve_samples(tmp_path):
    ticks = iter(float(i) for i in range(100))
    sample_index = iter(range(3))
    def sample(pid, job):
        index = next(sample_index)
        if index == 2:
            raise ProcessLookupError("fixture anchor gone")
        return copy.deepcopy(row(index)["snapshot"])
    output = tmp_path / "failure"
    with pytest.raises(ProcessLookupError):
        collect(456, 123, output, 3, 1., sample_fn=sample, anchor_fn=lambda p: 1.,
                host_fn=lambda: {}, clock=lambda: next(ticks), sleep=lambda _: None)
    report = json.loads((output / "results.json").read_text())
    assert report["status"] == "observation_failed" and report["observations_retained"] == 2
    assert len((output / "samples.jsonl").read_text().splitlines()) == 2
    ticks = iter(float(i) for i in range(100))
    result = collect(456, 123, tmp_path / "success", 2, 1., sample_fn=lambda p, j: row(1)["snapshot"],
                     anchor_fn=lambda p: 1., host_fn=lambda: {}, clock=lambda: next(ticks), sleep=lambda _: None)
    assert result["status"] == "bounded_observations_complete"


@pytest.mark.parametrize("count,interval", [(1, 1.), (2, 0.), (2, -1.), (2, float("nan")), (2, float("inf"))])
def test_invalid_collection_plan(tmp_path, count, interval):
    with pytest.raises(ValueError):
        collect(456, 123, tmp_path / "unused", count, interval)


def test_existing_output_is_not_reused(tmp_path):
    with pytest.raises(FileExistsError):
        collect(456, 123, tmp_path, 2, 1.)
