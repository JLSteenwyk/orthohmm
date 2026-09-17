from copy import deepcopy

import pytest

from benchmark_tools import admit_simulation_tree_panel as audit


def accounting(change=None):
    rows = [[f"{21405 if i == 0 else 21406}_{i}", str(30000 + i), "COMPLETED", "0:0", "00:01:00"]
            for i in range(210)]
    if change:
        change(rows)
    return "JobID|JobIDRaw|State|ExitCode|Elapsed\n" + "\n".join("|".join(r) for r in rows)


def test_exact_split_array_inventory():
    rows = audit.terminal_panel(accounting())
    assert len(rows) == 210
    assert rows[0]["JobID"] == "21405_0"
    assert rows[-1]["JobID"] == "21406_209"


@pytest.mark.parametrize("change", [
    lambda r: r.pop(),
    lambda r: r.append(r[0]),
    lambda r: r[0].__setitem__(0, "21406_0"),
    lambda r: r[10].__setitem__(2, "RUNNING"),
    lambda r: r[10].__setitem__(2, "PENDING"),
    lambda r: r[10].__setitem__(3, "1:0"),
])
def test_partial_duplicate_wrong_array_or_contradiction_rejected(change):
    with pytest.raises(ValueError):
        audit.terminal_panel(accounting(change))


def test_terminal_failures_remain_in_inventory():
    def change(rows):
        rows[5][2:4] = ["TIMEOUT", "0:15"]
    rows = audit.terminal_panel(accounting(change))
    assert len(rows) == 210 and rows[5]["State"] == "TIMEOUT"


def identity_fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, "record", lambda p: {"path": str(p)})
    task = {"JobID": "21406_1", "JobIDRaw": "999"}
    tree, config = {"variant": "nni1"}, {"methods": {"x": {"argv": ["frozen"]}}}
    provenance = {"index": 1, "tree": tree, "configured": config, "executor_commit": audit.EXECUTOR,
                  "job_id": "999", "array_job_id": "21406", "array_task_id": "1",
                  "source": audit.record(tmp_path / "benchmark_tools/run_simulation_tree_experiment.py"),
                  "helpers": [audit.record(tmp_path / "benchmark_tools" / n) for n in audit.HELPERS]}
    return provenance, task, tree, config


def test_identity_matches(tmp_path, monkeypatch):
    p, task, tree, config = identity_fixture(tmp_path, monkeypatch)
    audit.check_identity(p, 1, task, tree, config, tmp_path)


@pytest.mark.parametrize("key,value", [("index", 2), ("job_id", "other"), ("array_job_id", "21405"),
    ("array_task_id", "2"), ("executor_commit", "other"), ("helpers", []), ("source", {}),
    ("tree", {}), ("configured", {})])
def test_changed_identity_rejected(tmp_path, monkeypatch, key, value):
    p, task, tree, config = identity_fixture(tmp_path, monkeypatch)
    p = deepcopy(p)
    p[key] = value
    with pytest.raises(ValueError):
        audit.check_identity(p, 1, task, tree, config, tmp_path)


@pytest.mark.parametrize("state", ["FAILED", "TIMEOUT", "COMPLETED"])
def test_missing_preflight_never_becomes_admitted(tmp_path, state):
    tree = {"condition": "baseline", "seed": 1, "variant": "generating", "tree": {}}
    manifest = {"datasets": [{"label": "baseline_1"}]}
    args = (tmp_path, 0, {"State": state}, tree, manifest, {}, tmp_path, {}, {})
    if state == "COMPLETED":
        with pytest.raises(ValueError, match="missing preflight"):
            audit.validate_cell(*args)
    else:
        rows = audit.validate_cell(*args)
        assert {r["method"] for r in rows} == set(audit.METHODS)
        assert all(r["status"] == "failed" for r in rows)
