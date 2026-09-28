import pytest

from benchmark_tools.probe_threadripper_threads import validate_control


def fixture(mode):
    identities = [dict(pid=10+i*10, child=bool(i), mode=mode, leader_affinity=list(range(32)),
        threads=[dict(tid=11+i*10, affinity=[0, 96] if i and mode == "widen" else [0]),
                 dict(tid=12+i*10, affinity=list(range(32)))]) for i in range(2)]
    threads = []
    for row in identities:
        threads.append(dict(tid=row["pid"], affinity=row["leader_affinity"]))
        threads.extend(row["threads"])
    points = [dict(thread_affinity=dict(threads=threads, errors=[],
        status="violation" if mode == "widen" else "observed_within_affinity",
        violating_tids=[21] if mode == "widen" else []))]
    return identities, points


@pytest.mark.parametrize("mode", ["clean", "widen"])
def test_controls(mode):
    identities, points = fixture(mode)
    result = validate_control(mode, identities, points)
    assert result["complete_matching_points"] == [0]
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", ["missed_thread", "wrong_mask", "hidden_violation", "read_error", "wrong_child"])
def test_rejects_invalid_control(change):
    identities, points = fixture("widen")
    value = points[0]["thread_affinity"]
    if change == "missed_thread":
        value["threads"] = value["threads"][:-1]
    elif change == "wrong_mask":
        value["threads"][0]["affinity"] = [0]
    elif change == "hidden_violation":
        value["violating_tids"] = []
    elif change == "read_error":
        value["errors"] = [dict(error="read failed")]
    else:
        identities[1]["child"] = False
    with pytest.raises(ValueError):
        validate_control("widen", identities, points)
