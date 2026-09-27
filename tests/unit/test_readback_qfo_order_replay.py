import copy
import json

import pytest

from benchmark_tools.readback_qfo_order_replay import (
    PLAN_SHA, constraint_difference, partition_difference, read_partition,
    load_membership_constraints, verify_arm,
)


def test_partition_labels_order_and_changed_genes():
    left = [{"a", "b"}, {"c"}, {"d"}]
    assert partition_difference(left, list(reversed(left)))["equal"]
    result = partition_difference(left, [{"a"}, {"b", "c"}, {"d"}])
    assert result == dict(equal=False, left_groups=3, right_groups=3,
                         shared_groups=1, left_only_groups=2, right_only_groups=2,
                         genes_in_changed_groups=3)


def test_constraints_ignore_diagnostics_but_keep_direction_and_order():
    row = dict(source_genes=["a", "b"], target_genes=["c"], score=1)
    same = dict(source_genes=["b", "a"], target_genes=["c"], score=2)
    reverse = dict(source_genes=["c"], target_genes=["a", "b"])
    assert constraint_difference([row], [same])["equal"]
    assert not constraint_difference([row], [reverse])["equal"]
    assert constraint_difference([row, reverse], [reverse, row])["unequal_positions"] == 2
    assert constraint_difference([row], [row, reverse])["first_unequal_position"] == 1
    assert constraint_difference([], [])["first_unequal_position"] is None


@pytest.mark.parametrize("contents", ["a a\nb c\n", "a b\nb c\n", "a b\n", "a b c z\n"])
def test_invalid_complete_partition(tmp_path, contents):
    path = tmp_path / "partition"
    path.write_text(contents)
    with pytest.raises(ValueError):
        read_partition(path, {"a", "b", "c"})


@pytest.mark.parametrize("source,target", [(["a"], ["c"]), (["a"], ["a"]),
    (["a", "a"], ["b"]), ([], ["b"]), (["z"], ["b"])])
def test_invalid_merge_membership(tmp_path, source, target):
    partition, trace = tmp_path / "partition", tmp_path / "trace"
    partition.write_text("a b\nc\n")
    trace.write_text(json.dumps([dict(source_genes=source, target_genes=target)]))
    with pytest.raises(ValueError):
        load_membership_constraints(trace, partition)


def arm_fixture():
    source = dict(path="/repo/benchmark_tools/run_qfo_order_replay.py", sha256="source", bytes=1)
    pipeline = dict(path="/venv/pipeline.py", sha256="pipeline", bytes=2)
    plan = dict(directory="/run", repo="/repo", python="/venv/python", checked_records=[source, pipeline])
    identity = dict(path="/run/plan.json", sha256=PLAN_SHA, bytes=3)
    arm = "retained_order"
    execution = dict(arm=arm, returncode=0, command=["/usr/bin/time", "-v", "-o",
        "/run/retained_order.time.txt", "/venv/python", "-I", source["path"], "worker",
        "--plan", identity["path"], "--sha256", PLAN_SHA, "--arm", arm])
    started = dict(arm=arm, plan=identity, attempts=1, executable=plan["python"], source=source, pipeline=pipeline)
    complete = dict(arm=arm, plan=identity, status="candidate_arm_complete_pending_readback", accuracy_evaluated=False)
    return plan, identity, arm, execution, started, complete


def test_accept_bound_arm():
    verify_arm(*arm_fixture())


@pytest.mark.parametrize("index,key,value", [
    (3, "returncode", 1), (3, "command", []), (3, "arm", "canonical_order"),
    (4, "attempts", 2), (4, "executable", "/other/python"), (4, "plan", {}),
    (4, "pipeline", dict(path="/other", bytes=2, sha256="pipeline")),
    (5, "accuracy_evaluated", True), (5, "status", "partial"), (5, "plan", {}),
])
def test_reject_changed_arm(index, key, value):
    args = copy.deepcopy(arm_fixture())
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_arm(*args)
