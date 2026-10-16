from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools import fresh_protected_union_search as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json


def result():
    return SimpleNamespace(query_indices=np.array([0, 0], dtype=np.int32),
        target_indices=np.array([0, 1], dtype=np.int32), scores=np.array([0., 1.]),
        evalues=np.array([1e-5, 2e-5]), candidate_count=2)


def test_exact_complete_candidate_values_required():
    assert runner.require_same(result(), result(), ["q"], ["a", "b"])


@pytest.mark.parametrize("field,value", [("scores", [0., 2.]), ("scores", [-0., 1.]),
    ("evalues", [1e-4, 2e-5]), ("target_indices", [1, 0])])
def test_candidate_order_and_bitwise_values_must_match(field, value):
    actual = result()
    setattr(actual, field, np.array(value, dtype=getattr(actual, field).dtype))
    with pytest.raises(ValueError, match="candidate values differ"):
        runner.require_same(actual, result(), ["q"], ["a", "b"])


@pytest.mark.parametrize("field,dtype", [("scores", "float32"), ("evalues", "float32"),
    ("target_indices", "int64"), ("query_indices", "int64")])
def test_even_equal_values_need_native_layout(field, dtype):
    actual = result()
    setattr(actual, field, getattr(actual, field).astype(dtype))
    with pytest.raises(ValueError, match="candidate layout differs"):
        runner.require_same(actual, result(), ["q"], ["a", "b"])


def test_partial_candidate_snapshots_refused():
    actual = result()
    actual.candidate_count = 3
    with pytest.raises(ValueError, match="complete unfiltered"):
        runner.require_same(actual, result(), ["q"], ["a", "b"])


def test_observed_routes_must_cover_every_fresh_candidate(tmp_path):
    library = tmp_path / "kernel.so"
    library.write_bytes(b"trusted test library")
    pin = record(library)
    path = tmp_path / "routes_1.json"
    write_json(path, {"library": pin, "routes": {"initial_search": {"calls": 288, "scored_pairs": 4560971}}})
    assert runner.require_routes(tmp_path, pin, 4560971, 288) == record(path)
    with pytest.raises(ValueError, match="Actual kernel routes"):
        runner.require_routes(tmp_path, pin, 4560972, 288)
    with pytest.raises(ValueError, match="Actual kernel routes"):
        runner.require_routes(tmp_path, pin, 4560971, 287)
    write_json(tmp_path / "routes_2.json", {})
    with pytest.raises(ValueError, match="one observed"):
        runner.require_routes(tmp_path, pin, 4560971, 288)


@pytest.fixture
def preparation(tmp_path, monkeypatch):
    repo, source = tmp_path / "repo", tmp_path / "executor"
    tools, original = source / "benchmark_tools", repo / "original"
    tools.mkdir(parents=True)
    original.mkdir(parents=True)
    copies, old = [], []
    for name in runner.LIBRARIES:
        path = original / name
        path.write_bytes((name + " fixture").encode())
        copy = source / "orthohmm/search/csrc" / name
        copy.parent.mkdir(parents=True, exist_ok=True)
        copy.write_bytes(path.read_bytes())
        copies.append(record(copy))
        old.append(record(path))
    (original / "plan.json").write_text("{}\n")
    plan_pin = record(original / "plan.json")
    reports = [{"label": label, "execution_plan": plan_pin} for label in ("control", "shared_hit1")]
    for report in reports:
        path = original / report["label"] / "report.json"
        path.parent.mkdir()
        write_json(path, report)
    cached = repo / "cached.json"
    write_json(cached, {"checked": []})
    receipt = repo / "receipt.json"
    write_json(receipt, {"status": "cached_verified", "advance_gate_passed": True, "producer": record(cached)})
    protocol = {"prerequisite": {**record(receipt), "path": "receipt.json", "status": "cached_verified"}}
    for name in ("fresh_protected_union_search.py", runner.SLURM, runner.PROTOCOL,
                 "diagnose_protected_candidate_union.py", "strict_native_cli.py"):
        (tools / name).write_text("# frozen source fixture\n")
    write_json(tools / runner.PROTOCOL, protocol)
    execution = {"sources": old, "output": str(original), "runtime": {"files": []}}
    monkeypatch.setattr(runner, "context", lambda *args: ({}, deepcopy(execution), deepcopy(reports), [plan_pin, *old], "frozen"))
    monkeypatch.setattr(runner, "source_records", lambda *args: copies)
    def git(args, **kwargs):
        assert args[-2] == "show"
        return (source / args[-1].split(":", 1)[1]).read_bytes()
    monkeypatch.setattr(runner.subprocess, "check_output", git)
    return repo, source, tmp_path / "output", old, copies


def test_preparation_pins_original_libraries_and_both_retained_arms(preparation):
    repo, source, output, old, copies = preparation
    plan = runner.prepare(repo, source, output)
    assert plan["libraries"] == [{"original": a, "copy": b} for a, b in zip(old, copies)]
    assert len(plan["original_reports"]) == 2 and plan["threads"] == 8
    assert plan == json.loads((output / "execution_plan.json").read_text())
    assert plan["defaults_promoted"] is plan["goal_completed"] is False


def test_preparation_rejects_changed_control_library(preparation):
    repo, source, output, _, _ = preparation
    (source / "orthohmm/search/csrc/hmm_viterbi.so").write_bytes(b"different library")
    with pytest.raises(ValueError, match="CPU library differs"):
        runner.prepare(repo, source, output)
    assert not output.exists()


def test_preparation_rejects_changed_committed_executor(preparation, monkeypatch):
    repo, source, output, _, _ = preparation
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *args, **kwargs: b"other committed content")
    with pytest.raises(ValueError, match="committed frozen source"):
        runner.prepare(repo, source, output)


def test_preparation_does_not_replace_an_existing_attempt(preparation):
    repo, source, output, _, _ = preparation
    output.mkdir()
    with pytest.raises(ValueError, match="fresh diagnostic namespace"):
        runner.prepare(repo, source, output)
