import copy

import pytest

from benchmark_tools.probe_cgroup_lineage import scopes
from benchmark_tools.run_lineage_read_crossing_control import assess, run


def point(times, counters):
    target = "/system.slice/slurm.scope/job_1"
    names = scopes(target)
    identities = {name: [1, i + 1] for i, name in enumerate(names)}
    return dict(status="aggregate_lineage_snapshot", target=target,
        boot_before="boot", boot_after="boot", identities_before=identities,
        identities_after=copy.deepcopy(identities), rows=[dict(scope=name,
            started_ns=t, finished_ns=t+1, raw=f"usage_usec {v}\n")
            for name, t, v in zip(names, times, counters)])


@pytest.fixture
def evidence():
    return [point([0, 10, 20, 30], [100, 80, 60, 40]),
            point([100, 200, 210, 220], [110, 100, 80, 60]),
            point([300, 310, 320, 330], [750200, 120, 100, 80]),
            dict(started_ns=110, finished_ns=190),
            dict(status="lifecycle_control_evaluated", result=dict(
                descendant_cpu_retention_response=True, outside_target_response=True))]


def test_crossing_retains_negative_partial_span_without_promoting_timing(evidence):
    result = assess(*evidence)
    assert result["control_met"] is True
    assert result["spans"]["before_to_crossing"]["root_minus_target_cpu_usec"] == -10
    assert result["spans"]["full"]["root_minus_target_cpu_usec"] == 750060
    assert result["scientific_timings_admitted"] is False
    assert result["environmental_validity_established"] is False


@pytest.mark.parametrize("fault", ["early", "late", "reversed", "bool", "nested", "boot", "identity"])
def test_invalid_crossing_rejected(evidence, fault):
    event = evidence[3]
    if fault == "early":
        event["started_ns"] = 100
    elif fault == "late":
        event["finished_ns"] = 201
    elif fault == "reversed":
        event["finished_ns"] = 109
    elif fault == "bool":
        event["started_ns"] = True
    elif fault == "nested":
        evidence[4]["result"]["outside_target_response"] = False
    elif fault == "boot":
        evidence[2]["boot_before"] = "other"
    else:
        evidence[2]["identities_after"]["/"][1] += 1
    with pytest.raises(ValueError):
        assess(*evidence)


def test_insufficient_full_span_is_retained_failed_control(evidence):
    evidence[2]["rows"][0]["raw"] = "usage_usec 1000\n"
    assert assess(*evidence)["control_met"] is False


def test_local_host_cannot_launch_owned_service(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    with pytest.raises(ValueError, match="Require exclusive DGX"):
        run(tmp_path / "output")
    assert not (tmp_path / "output").exists()
