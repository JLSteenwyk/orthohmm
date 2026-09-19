import pytest

from benchmark_tools import run_native_pressure_controls as module
from tests.unit.test_probe_native_pressure import point


def row(mode="contended"):
    points = [point(0), point(3)]
    native = dict(pid=42, cgroup=points[0]["native_membership"], affinity=[1],
                  started_ns=1000, finished_ns=2_000_000_000, cpu_seconds=.76)
    competitor = dict(native, pid=43, cgroup=points[0]["host"][0]["raw"]["cgroup_membership"],
                      started_ns=1000000)
    return dict(block=0, mode=mode, points=points, ready=dict(cpu=1, allowed=[1], pid=42,
                membership=native["cgroup"]), observer_allowed=[1, 2], observer_cpu=2,
                worker_exit_code=0, native=native, competitor=competitor if mode == "contended" else None,
                pressure=module.compare(*points, 123))


@pytest.mark.parametrize("mode", ["quiet", "native-only", "contended"])
def test_injection_witness(mode):
    result = module.validate_trial(row(mode), 123)
    assert result["status"] == "control_injection_validated"
    assert (result["overlap_s"] is not None) == (mode == "contended")


@pytest.mark.parametrize("fault", ["cpu", "observer", "dose", "pid", "membership", "overlap",
                                    "competitor_scope", "competitor_dose", "exit", "pressure", "bracket"])
def test_reject_invalid_witness(fault):
    value = row()
    if fault == "cpu":
        value["native"]["affinity"] = [2]
    elif fault == "observer":
        value["observer_cpu"] = 1
    elif fault == "dose":
        value["native"]["cpu_seconds"] = .5
    elif fault == "pid":
        value["native"]["pid"] = 100
    elif fault == "membership":
        value["ready"]["membership"] = "0::/other\n"
    elif fault == "overlap":
        value["competitor"]["started_ns"] = 1_900_000_000
    elif fault == "competitor_scope":
        value["competitor"]["cgroup"] = value["native"]["cgroup"]
    elif fault == "competitor_dose":
        value["competitor"]["cpu_seconds"] = 1
    elif fault == "exit":
        value["worker_exit_code"] = 1
    elif fault == "pressure":
        value["pressure"]["native_stall_usec"]["cpu"]["some"] += 1
    else:
        value["native"]["started_ns"] = 0
    with pytest.raises(ValueError):
        module.validate_trial(value, 123)


def rows():
    result = []
    for block, modes in enumerate(module.ORDER):
        for mode in modes:
            value = row(mode)
            value.update(block=block, status="validated")
            value["pressure"]["native_stall_usec"]["cpu"]["some"] = 200001 if mode == "contended" else 1
            result.append(value)
    return result


def test_all_blocks_must_pass_no_scientific_admission():
    result = module.summarize(rows())
    assert result["all_response_checks_passed"] is True
    assert result["scientific_timings_admitted"] is False
    assert [b["difference_usec"] for b in result["blocks"]] == [200000]*3


def test_failed_control_not_zero_filled_or_dropped():
    values = rows()
    values[0]["status"] = "failed"
    result = module.summarize(values)
    assert result["all_response_checks_passed"] is False
    assert result["blocks"][0]["difference_usec"] is None


def test_single_failed_response_blocks_summary():
    values = rows()
    values[2]["pressure"]["native_stall_usec"]["cpu"]["some"] = 100000
    assert module.summarize(values)["all_response_checks_passed"] is False


@pytest.mark.parametrize("change", ["missing", "reordered"])
def test_frozen_order_required(change):
    values = rows()
    if change == "missing":
        values.pop()
    else:
        values[0], values[1] = values[1], values[0]
    with pytest.raises(ValueError, match="inventory/order"):
        module.summarize(values)


def test_protocol_hash_still_frozen():
    import hashlib
    path = module.Path(module.__file__).parent / "results/NATIVE_PRESSURE_CONTROL_PROTOCOL_20260918.md"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == module.PROTOCOL_SHA
