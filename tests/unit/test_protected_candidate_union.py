from copy import deepcopy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.diagnose_protected_candidate_union import advance_gate, protected_union
from benchmark_tools import diagnose_protected_candidate_union as diagnostic
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json


def result(pairs, scores, evalues):
    return SimpleNamespace(query_indices=np.array([p[0] for p in pairs], dtype=np.int32),
        target_indices=np.array([p[1] for p in pairs], dtype=np.int32), candidate_count=len(pairs),
        scores=np.array(scores, dtype=np.float64), evalues=np.array(evalues, dtype=np.float64))


@pytest.fixture
def arms():
    return (result([(0, 1), (0, 0), (1, 1)], [3, 2, 1], [1e-5, .5, 1e-5]),
            result([(0, 2), (0, 1), (1, 0)], [4, 3, 5], [1e-5, 1e-5, .5]))


def test_union_preserves_baseline_and_relaxed_candidates_and_values(arms):
    merged, stats = protected_union(*arms, ["q1", "q2"], ["t1", "t2", "t3"])
    assert list(zip(merged.query_indices, merged.target_indices)) == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    assert merged.scores.tolist() == [2, 3, 4, 5, 1]
    assert stats == {"baseline": 3, "relaxed": 3, "common": 1, "baseline_only": 2,
        "relaxed_only": 2, "union": 5, "per_query_union_histogram": {"2": 1, "3": 1},
        "common_values_bitwise_equal": True}


def test_union_does_not_mutate_inputs_and_is_order_independent(arms):
    saved = deepcopy(arms)
    forward, _ = protected_union(*arms, ["q1", "q2"], ["t1", "t2", "t3"])
    reverse, _ = protected_union(*reversed(arms), ["q1", "q2"], ["t1", "t2", "t3"])
    for a, b in zip(arms, saved):
        for field in ("query_indices", "target_indices", "scores", "evalues"):
            assert np.array_equal(getattr(a, field), getattr(b, field))
            assert np.array_equal(getattr(forward, field), getattr(reverse, field))


@pytest.mark.parametrize("field", ["scores", "evalues"])
def test_common_candidate_difference_refused(arms, field):
    getattr(arms[1], field)[1] += .001
    with pytest.raises(ValueError, match="bit-for-bit identical"):
        protected_union(*arms, ["q1", "q2"], ["t1", "t2", "t3"])


def test_common_signed_zero_difference_refused():
    a, b = result([(0, 0)], [0.], [1.]), result([(0, 0)], [-0.], [1.])
    with pytest.raises(ValueError, match="bit-for-bit identical"):
        protected_union(a, b, ["q"], ["t"])


@pytest.mark.parametrize("field,value", [("scores", float("nan")), ("scores", float("inf")),
    ("evalues", -1), ("evalues", float("inf"))])
def test_invalid_values_refused(arms, field, value):
    getattr(arms[0], field)[0] = value
    with pytest.raises(ValueError, match="Invalid numerical"):
        protected_union(*arms, ["q1", "q2"], ["t1", "t2", "t3"])


def test_original_cap_and_duplicate_candidates_refused():
    a = result([(0, i) for i in range(101)], [1.] * 101, [1.] * 101)
    with pytest.raises(ValueError, match="100-candidate"):
        protected_union(a, a, ["q"], [str(i) for i in range(101)])
    a = result([(0, 0), (0, 0)], [1., 1.], [1., 1.])
    with pytest.raises(ValueError, match="Duplicate directed candidate"):
        protected_union(a, a, ["q"], ["t"])


def test_empty_and_exact_200_candidate_bound():
    empty = result([], [], [])
    union, stats = protected_union(empty, empty, ["q"], ["t"])
    assert union.candidate_count == 0 and stats["per_query_union_histogram"] == {"0": 1}
    a = result([(0, i) for i in range(100)], [1.] * 100, [1.] * 100)
    b = result([(0, i) for i in range(100, 200)], [1.] * 100, [1.] * 100)
    union, stats = protected_union(a, b, ["q"], [str(i) for i in range(200)])
    assert union.candidate_count == 200 and stats["per_query_union_histogram"] == {"200": 1}


@pytest.mark.parametrize("union,gained,lost,expected", [(125, 1, 0, True), (126, 1, 0, False),
    (125, 0, 0, False), (125, 1, 1, False)])
def test_prespecified_gate_is_not_retuned(union, gained, lost, expected):
    assert advance_gate({"baseline": 100, "union": union}, [
        {"versus_control": {"gained": gained, "lost": lost}}]) is expected


@pytest.fixture
def frozen_context(tmp_path, monkeypatch):
    source, repo = tmp_path / "executor", tmp_path / "repo"
    tools = source / "benchmark_tools"
    tools.mkdir(parents=True)
    repo.mkdir()
    for name in diagnostic.SOURCES:
        (tools / name).write_text("# frozen fixture\n")
    inputs = repo / "input.json"
    write_json(inputs, {})
    pin = record(inputs)
    root = repo / "original"
    root.mkdir()
    plan = {"source_revision": "original", "output": str(root), "inputs": [pin],
            "sources": [pin], "runtime": {"files": [pin]}}
    path = root / "execution_plan.json"
    write_json(path, plan)
    plan_pin = record(path)
    reports, arms = [], []
    for label, minimum in (("control", 4), ("shared_hit1", 1)):
        folder = root / label
        folder.mkdir()
        report = {"label": label, "execution_plan": plan_pin, "source_revision": "original",
            "status": "fresh_search_completed_pending_independent_audit", "explicit_control_fixture_exact": True,
            "settings": {"min_total_hits": minimum, "fixed": True}}
        write_json(folder / "report.json", report)
        reports.append(report)
        arms.append({"label": label, "status": "independently_verified", "retained_report": record(folder / "report.json")})
    audit = repo / "audit.json"
    write_json(audit, {"execution_plan": plan_pin, "status": "all_search_arms_independently_verified", "arms": arms})
    spec = {"prior_evidence": [{**record(audit), "path": "audit.json"}], "fixed_settings": {"fixed": True},
        "inputs": {"execution_plan": "original/execution_plan.json", "execution_plan_sha256": plan_pin["sha256"],
            "scientific_revision": "original", "arms": ["control", "shared_hit1"],
            "control_report_sha256": arms[0]["retained_report"]["sha256"],
            "relaxed_report_sha256": arms[1]["retained_report"]["sha256"]}}
    write_json(tools / diagnostic.PROTOCOL, spec)
    def git(args, **kwargs):
        if args[-2:] == ["rev-parse", "HEAD"]:
            return "executor\n"
        assert args[-2] == "show"
        return (source / args[-1].split(":", 1)[1]).read_bytes()
    monkeypatch.setattr(diagnostic.subprocess, "check_output", git)
    monkeypatch.setattr(diagnostic, "runtime_records", lambda: {"files": [pin]})
    return source, repo, spec, plan, reports


def test_context_authenticates_original_plan_audit_arm_reports_and_source(frozen_context):
    source, repo, spec, plan, reports = frozen_context
    actual, execution, arms, pins, revision = diagnostic.context(repo, source)
    assert actual == spec and execution == plan and arms == reports
    assert revision == "executor" and all(p["sha256"] for p in pins)


@pytest.mark.parametrize("change", ["audit", "plan", "arm", "runtime", "source"])
def test_changed_context_refused(frozen_context, monkeypatch, change):
    source, repo, spec, _, _ = frozen_context
    if change == "runtime":
        monkeypatch.setattr(diagnostic, "runtime_records", lambda: {})
    elif change == "source":
        monkeypatch.setattr(diagnostic.subprocess, "check_output", lambda args, **kwargs:
            "executor\n" if args[-2:] == ["rev-parse", "HEAD"] else b"wrong committed bytes")
    else:
        target = {"audit": repo / "audit.json", "plan": repo / "original/execution_plan.json",
                  "arm": repo / "original/control/report.json"}[change]
        data = json.loads(target.read_text())
        data["changed"] = True
        write_json(target, data)
    with pytest.raises(ValueError):
        diagnostic.context(repo, source)
