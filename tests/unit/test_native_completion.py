from copy import deepcopy
import json

import pytest

from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.replay_threadripper_scaling import replay_completion


def observation(start, end):
    return dict(status="observed_within_affinity", errors=[], violating_tids=[],
                initial_tids=[9], final_tids=[9], scope="/job_1/step_0/user",
                threads=[dict(tid=9, start_ticks=100)], started_ns=start, finished_ns=end)


def test_anchor_only():
    result = completion_evidence(observation(1, 2), observation(4, 5), 9, 3)
    assert result["errors"] == []
    assert result["status"] == "anchor_only_at_boundaries"
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", [
    {"initial_tids": [9, 10]}, {"final_tids": [9, 10]}, {"initial_tids": []},
    {"errors": ["unreadable"]}, {"violating_tids": [9]}, {"status": "incomplete"},
    {"threads": []}, {"threads": [dict(tid=9, start_ticks=101)]},
    {"threads": [dict(tid=9, start_ticks=True)]}, {"scope": "/other"},
    {"started_ns": 2}, {"finished_ns": 3}, {"started_ns": float("nan")},
])
def test_bad_final_boundary(change):
    after = observation(4, 5)
    after.update(change)
    assert completion_evidence(observation(1, 2), after, 9, 3)["errors"]


def test_bad_initial_inventory():
    before = observation(1, 2)
    before["initial_tids"].append(10)
    assert completion_evidence(before, observation(4, 5), 9, 3)["errors"]


def test_input_not_changed():
    before, after = observation(1, 2), observation(4, 5)
    original = deepcopy((before, after))
    completion_evidence(before, after, 9, 3)
    assert (before, after) == original


@pytest.mark.parametrize("pid", [0, -1, True, 9.0])
def test_bad_pid(pid):
    with pytest.raises(ValueError):
        completion_evidence(observation(1, 2), observation(4, 5), pid, 3)


@pytest.mark.parametrize("corrupt", [None, "missing", "embedded", "file", "descendant", "symlink"])
def test_replay_completion(tmp_path, corrupt):
    before, after = observation(1, 2), observation(4, 5)
    receipt = completion_evidence(before, after, 9, 3)
    measured = dict(native_completion=deepcopy(receipt))
    path = tmp_path / "native_completion.json"
    if corrupt == "file":
        receipt["anchor_pid"] = 10
    if corrupt == "embedded":
        measured["native_completion"]["anchor_pid"] = 10
    if corrupt == "descendant":
        after["final_tids"].append(10)
    if corrupt == "symlink":
        target = tmp_path / "other.json"
        target.write_text(json.dumps(receipt))
        path.symlink_to(target)
    elif corrupt != "missing":
        path.write_text(json.dumps(receipt))
    args = (tmp_path, measured, [{"thread_affinity": before}, {"thread_affinity": after}],
            {"pid": 9}, {"finished_ns": 3})
    if corrupt:
        with pytest.raises(ValueError):
            replay_completion(*args)
    else:
        replayed, record = replay_completion(*args)
        assert replayed == receipt and record["sha256"]
