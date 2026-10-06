from copy import deepcopy
import json

import pytest

from benchmark_tools import review_native_factorial_measurement_failure as module
from benchmark_tools.prepare_native_factorial_qfo_pairs import admit_conversion
from tests.unit.test_review_native_factorial_attempt import runtime  # noqa: F401


def fixture():
    ref = dict(path="/request.json", bytes=1, sha256="a")
    plan = dict(path="/plan.json", bytes=2, sha256="b")
    request = dict(job_id=22437, index=7, plan=plan)
    run = dict(index=7, cell="p0_c0_r1")
    pressure = dict(native_pressure_role="diagnostic_only", pressure_thresholds_used_for_eligibility=False,
        sampled_pressure_evidence_satisfied=False, failures={"pressure_sample_period_exceeded": 1})
    audit = dict(schema="native_factorial_cadence_failure_audit_v1",
        status="retained_measurement_cadence_failure_reproduced", request=ref, plan=plan,
        job_id=22437, index=7, cell="p0_c0_r1", primary_resources=None,
        census=dict(failures=[dict(error=module.ERROR, wall_s=1.66)],
            failure_counts={module.ERROR: 1}, unchanged_cadence_bounds_s=[.5, 1.5]),
        retained_pressure_review=deepcopy(pressure))
    for key in ("full_resource_replay", "scientific_timings_admitted", "native_outputs_validated",
        "next_identity_authorized", "original_receipts_rewritten", "inference_reexecuted", "automatic_retry"):
        audit[key] = False
    terminal = dict(verified=dict(State="FAILED", ExitCode="1:0"))
    wrapper = dict(status="verified_wrapper_failed", error=module.ERROR)
    done = dict(exit_code=0, timed_out=False)
    environment = dict(job_id=22437, index=7, execution_scope=module.SCOPE,
        uncontended_timing=False, sampled_process_policy_satisfied=True,
        sampled_environment_policy_satisfied=False, background_cpu_used_for_eligibility=False,
        failures={}, pressure_review=pressure)
    return [audit, ref, request, run, terminal, wrapper, done, environment]


def test_exact_classified_failure_scope():
    assert module.failure_scope(*fixture()) is None


@pytest.mark.parametrize("change", ["scheduler", "exit", "audit_schema", "job", "index", "cell",
    "request", "resource_imputation", "resource_admission", "retry", "native_exit", "native_exit_bool",
    "timeout", "wrapper", "other_census_error", "relaxed_bounds", "normal_cadence", "process_failure",
    "pressure_unknown", "pressure_different", "background_exclusion", "environment_success"])
def test_rejects_unclassified_or_relabelled_failure(change):
    args = fixture()
    audit, _, _, _, terminal, wrapper, done, environment = args
    if change == "scheduler": terminal["verified"]["State"] = "COMPLETED"
    elif change == "exit": terminal["verified"]["ExitCode"] = "0:0"
    elif change == "audit_schema": audit["schema"] = "native_factorial_terminal_review_v1"
    elif change in {"job", "index", "cell"}: audit[{"job": "job_id"}.get(change, change)] = "different"
    elif change == "request": audit["request"] = {}
    elif change == "resource_imputation": audit["primary_resources"] = dict(cpu_s=0)
    elif change == "resource_admission": audit["scientific_timings_admitted"] = True
    elif change == "retry": audit["automatic_retry"] = True
    elif change == "native_exit": done["exit_code"] = 1
    elif change == "native_exit_bool": done["exit_code"] = False
    elif change == "timeout": done["timed_out"] = True
    elif change == "wrapper": wrapper["error"] = "unknown error"
    elif change == "other_census_error": audit["census"]["failure_counts"]["decreasing counter"] = 1
    elif change == "relaxed_bounds": audit["census"]["unchanged_cadence_bounds_s"] = [.1, 2.]
    elif change == "normal_cadence": audit["census"]["failures"][0]["wall_s"] = 1.
    elif change == "process_failure": environment["sampled_process_policy_satisfied"] = False
    elif change == "pressure_unknown": environment["pressure_review"]["failures"]["unknown"] = 1
    elif change == "pressure_different": audit["retained_pressure_review"] = {}
    elif change == "background_exclusion": environment["background_cpu_used_for_eligibility"] = True
    else: environment["sampled_environment_policy_satisfied"] = True
    with pytest.raises(ValueError):
        module.failure_scope(*args)


def original_runtime(data, tmp_path):
    evidence = module.Evidence()
    value = module.runtime_review(data["plan"], data["run"], data["session"], data["verification"], evidence)
    path = tmp_path / "original_runtime_review.json"
    path.write_text(json.dumps(value))
    return module.record(path)


def test_reused_runtime_is_not_a_new_fresh_tree_audit(runtime, tmp_path):
    ref = original_runtime(runtime, tmp_path)
    evidence = module.Evidence()
    result = module.reuse_runtime(runtime["plan"], runtime["run"], runtime["session"],
        runtime["verification"], ref, evidence)
    evidence.finish()
    assert result["original_runtime_review"] == ref
    assert result["original_terminal_inventory_reused"] is True
    assert result["fresh_runtime_tree_recheck"] is False
    assert result["continuous_runtime_integrity_established"] is False


@pytest.mark.parametrize("change", ["input", "lookup", "runtime_record", "status"])
def test_runtime_reuse_rejects_changed_bindings(runtime, tmp_path, change):
    ref = original_runtime(runtime, tmp_path)
    if change == "input":
        (runtime["source"].parent / "run/input/source.fa").write_text(">gene\nDIFFERENT\n")
    elif change == "lookup":
        path = runtime["session"] / "lookup_checks/check_02/orthohmm.json"
        value = json.loads(path.read_text())
        value["modules"] = {"orthohmm": "/foreign/core.py"}
        path.write_text(json.dumps(value))
    else:
        path = tmp_path / "original_runtime_review.json"
        value = json.loads(path.read_text())
        if change == "status": value["status"] = "failed"
        else: value["fresh_terminal_inventory"][0]["records"] = 999
        path.write_text(json.dumps(value))
        ref = module.record(path)
    with pytest.raises(ValueError):
        module.reuse_runtime(runtime["plan"], runtime["run"], runtime["session"],
            runtime["verification"], ref, module.Evidence())


def test_existing_success_only_converter_still_rejects_failure_recovery():
    audit, request_ref, request, run, *_ = fixture()
    run.update(dataset="qfo_corrected", repeat=0)
    recovered = dict(schema="native_factorial_measurement_failure_review_v1", request=request_ref,
        plan=request["plan"], job_id=request["job_id"], index=7, cell=run["cell"], dataset=run["dataset"],
        repeat=0, status="native_scientific_outputs_recovered_measurement_failure_retained",
        native_outputs_validated=True, terminal_reviewed=True, primary_resources_replayed=False)
    with pytest.raises(ValueError, match="successful"):
        admit_conversion(recovered, request_ref, request, run)
