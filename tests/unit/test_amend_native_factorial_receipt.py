"""Cross-plan scope and refusal tests; no new native inference or submission."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import amend_native_factorial_receipt as amend
from benchmark_tools import run_native_factorial_cost as cost
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def plans(monkeypatch):
    old = json.loads((amend.ROOT / "benchmark_tools/results/native_factorial_cost_plan_repaired_20261004.json").read_text())
    new = deepcopy(old)
    new["source_commit"] = "prospective"
    new["helper_sources"] = [record(amend.EXECUTOR) if r["path"] == str(amend.EXECUTOR) else deepcopy(r)
                             for r in old["helper_sources"]]
    new["evidence"].append(dict(path="/prospective", bytes=1, sha256="a" * 64))
    monkeypatch.setattr(amend, "historical_check", lambda ref: None)
    monkeypatch.setattr(amend, "check", lambda ref: None)
    return old, new


def test_exact_scientific_invocation_preserved():
    amend.scientific_source_unchanged()


def test_prospective_control_only_plan(plans):
    old, new = plans
    amend.compare_plans(old, new)
    assert old["runs"] == new["runs"] and old["resources"] == new["resources"]


@pytest.mark.parametrize("key", ["baseline", "runtime_lookup", "core_commit", "resources",
    "panel_root", "filename_enumeration_probes", "replacement_reason", "automatic_retry"])
def test_no_science_resource_or_history_refresh(plans, key):
    old, new = plans
    new[key] = "changed"
    with pytest.raises((ValueError, TypeError)):
        amend.compare_plans(old, new)


@pytest.mark.parametrize("key", ["inputs", "native_order", "input_creation_order", "repeat", "output_root"])
def test_no_identity_or_input_mutation(plans, key):
    old, new = plans
    new["runs"][1][key] = "changed"
    with pytest.raises((ValueError, TypeError)):
        amend.compare_plans(old, new)


def test_no_unrelated_helper_or_evidence_change(plans):
    old, new = plans
    other = next(r for r in new["helper_sources"] if r["path"] != str(amend.EXECUTOR))
    other["sha256"] = "a" * 64
    with pytest.raises(ValueError, match="helper"):
        amend.compare_plans(old, new)
    new["helper_sources"] = deepcopy(old["helper_sources"])
    with pytest.raises(ValueError, match="helper"):
        amend.compare_plans(old, new)
    new["helper_sources"] = [record(amend.EXECUTOR) if r["path"] == str(amend.EXECUTOR) else r
                             for r in old["helper_sources"]]
    new["evidence"].pop(0)
    with pytest.raises(ValueError, match="evidence"):
        amend.compare_plans(old, new)


@pytest.fixture
def adoption(plans, monkeypatch):
    old, new = plans
    recovery_ref = record(amend.ROOT / "benchmark_tools/results/native_factorial_recovery_22427.json")
    recovery = json.loads(Path(recovery_ref["path"]).read_text())
    prior_ref = recovery["terminal_review"]
    prior = json.loads(Path(prior_ref["path"]).read_text())
    new["history_adoption"] = dict(schema="native_factorial_receipt_history_adoption_v1",
        prior_plan=prior["plan"], prior_review=prior_ref, recovery=recovery_ref, executed_source=record(amend.ARCHIVE))
    values = {prior_ref["path"]: prior, prior["plan"]["path"]: old, recovery_ref["path"]: recovery,
              recovery["outputs"]["path"]: dict(input_genes=251378)}
    monkeypatch.setattr(amend, "read", lambda ref: values[ref["path"]])
    return prior_ref, prior, new, recovery


def test_cross_plan_adopts_only_reviewed_failed_zero(adoption):
    ref, prior, plan, recovery = adoption
    result = amend.adopt(ref, prior, plan)
    assert result["status"] == "reviewed_failed_index_zero_adopted_without_retry"
    assert result["scheduler_success"] is result["timing_success_established"] is False
    assert result["inference_reexecuted"] is False and result["retained_evidence_checked"] > 100
    assert prior["plan"] == plan["history_adoption"]["prior_plan"]


@pytest.mark.parametrize("key,value", [("index", 1), ("job_id", 999), ("cell", "p0_c0_r1"),
    ("status", "native_success"), ("scheduler_state", "COMPLETED"), ("scheduler_exit_code", "0:0"),
    ("terminal_reviewed", False), ("next_identity_authorized", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False)])
def test_unreviewed_or_other_attempt_refused(adoption, key, value):
    ref, prior, plan, recovery = adoption
    prior[key] = value
    with pytest.raises(ValueError):
        amend.adopt(ref, prior, plan)


@pytest.mark.parametrize("key,value", [("scheduler_success", True), ("timing_success_established", True),
    ("original_receipts_rewritten", True), ("inference_reexecuted", True), ("automatic_retry", True),
    ("native_outputs_validated", False), ("accuracy_evaluated", False), ("cell", "p0_c1_r0")])
def test_recovery_scope_changes_refused(adoption, key, value):
    ref, prior, plan, recovery = adoption
    recovery[key] = value
    with pytest.raises(ValueError, match="Recovery"):
        amend.adopt(ref, prior, plan)


def test_conflicting_evidence_refused(adoption):
    ref, prior, plan, recovery = adoption
    recovery["evidence"].append(dict(prior["source"], sha256="a" * 64))
    with pytest.raises(ValueError, match="Conflicting"):
        amend.adopt(ref, prior, plan)


def test_no_implicit_cross_plan_acceptance(adoption):
    ref, prior, plan, recovery = adoption
    plan.pop("history_adoption")
    with pytest.raises(ValueError, match="Unapproved"):
        amend.adopt(ref, prior, plan)


def test_existing_destination_preserved(tmp_path):
    with pytest.raises(ValueError, match="results destination"):
        amend.prepare(tmp_path)
    destination = amend.ROOT / "benchmark_tools/results"
    with pytest.raises(FileExistsError):
        amend.prepare(destination)


def test_no_arbitrary_added_helper(plans):
    old, new = plans
    new["helper_sources"].append(dict(path="/unapproved.py", bytes=1, sha256="a" * 64))
    with pytest.raises(ValueError, match="helper"):
        amend.compare_plans(old, new)


def test_controller_refuses_retry_of_adopted_zero():
    with pytest.raises(ValueError, match="must not be retried"):
        cost.reviewed_history(dict(index=0, history=[]), {}, dict(history_adoption={}))


@pytest.mark.parametrize("cross_plan", [False, True])
@pytest.mark.parametrize("state,code,success", [("FAILED", "1:0", True),
    ("COMPLETED", "0:0", False), ("FAILED", "2:0", False)])
def test_controller_queries_actual_previous_terminal(monkeypatch, cross_plan, state, code, success):
    old_ref, new_ref = dict(path="/old"), dict(path="/new")
    prior_ref = dict(path="/review")
    prior = dict(index=0, job_id=22427, plan=old_ref if cross_plan else new_ref,
        terminal_reviewed=True, next_identity_authorized=True, scheduler_state="FAILED", scheduler_exit_code="1:0")
    monkeypatch.setattr(cost, "read", lambda ref: prior)
    jobs, adopted = [], []
    def terminal(job):
        jobs.append(job)
        return dict(verified=dict(State=state, ExitCode=code))
    def adopt(ref, review, plan):
        adopted.append((ref, review))
        return dict(status="explicit_test_adoption")
    monkeypatch.setattr(cost, "verify_terminal", terminal)
    monkeypatch.setattr(amend, "adopt", adopt)
    if success:
        result = cost.reviewed_history(dict(index=1, history=[prior_ref]), new_ref, {})
        assert result[0]["adoption"] == (dict(status="explicit_test_adoption") if cross_plan else None)
    else:
        with pytest.raises(ValueError, match="outcome changed"):
            cost.reviewed_history(dict(index=1, history=[prior_ref]), new_ref, {})
    assert jobs == [22427] and len(adopted) == int(cross_plan)


def test_controller_refuses_unresolved_history_before_query(monkeypatch):
    monkeypatch.setattr(cost, "read", lambda ref: dict(index=0, terminal_reviewed=False))
    monkeypatch.setattr(cost, "verify_terminal", lambda job: pytest.fail("Must not query unresolved identity"))
    with pytest.raises(ValueError, match="unresolved"):
        cost.reviewed_history(dict(index=1, history=[{}]), {}, {})
