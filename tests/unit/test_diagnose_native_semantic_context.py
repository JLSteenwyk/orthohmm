"""Invented failed-call contexts, without running selected scientific data."""

import json

import pytest

from benchmark_tools import diagnose_native_semantic_context as diagnostic


def test_passing_kernel_is_nested_and_does_not_authorize_admission():
    context = dict(example="invented")
    calls = []

    def validator(value):
        calls.append(value)
        return dict(native_outputs_validated=True, terminal_reviewed=False, accuracy_evaluated=False)

    result = diagnostic.evaluate(context, validator)
    assert calls == [context]
    assert result["status"] == "semantic_probe_passed"
    assert result["semantic_result"]["native_outputs_validated"] is True
    assert all(result[key] is False for key in (
        "native_outputs_validated", "accuracy_evaluated", "terminal_reviewed", "next_identity_authorized",
        "automatic_retry", "resources_admitted", "historical_failure_cause_established"))


@pytest.mark.parametrize("owners,value,missing,extra", [
    ({"a": "s1", "b": "s2"}, {"g": ["a"]}, 1, 0),
    ({"a": "s1"}, {"g": ["a", "x"]}, 0, 1),
    ({"a": "s1", "b": "s1"}, {"g": ["x"]}, 2, 1),
])
def test_exception_frames_capture_the_actual_failed_coverage(owners, value, missing, extra):
    def validator(context):
        local_owners = owners
        local_value = value

        def fail(owners, value):
            raise ValueError("Native partition loses or adds input genes")

        return fail(local_owners, local_value)

    result = diagnostic.evaluate({}, validator)
    assert result["status"] == "semantic_probe_failed"
    assert result["semantic_result"] is None
    assert result["error_type"] == "ValueError"
    failed = result["failure_frames"][-1]
    assert failed["owner_count"] == len(owners)
    assert failed["missing_count"] == missing
    assert failed["extra_count"] == extra
    assert result["historical_failure_cause_established"] is False


def test_noncoverage_error_preserved_without_invented_owners():
    def validator(context):
        raise RuntimeError("invented other stage failure")

    result = diagnostic.evaluate({}, validator)
    assert result["error_type"] == "RuntimeError"
    assert result["error"] == "invented other stage failure"
    assert "owner_count" not in result["failure_frames"][-1]
    assert "RuntimeError" in result["traceback"]


@pytest.mark.parametrize("error", [KeyboardInterrupt, SystemExit])
def test_user_or_platform_interruptions_are_not_turned_into_completed_probes(error):
    def validator(context):
        raise error()

    with pytest.raises(error):
        diagnostic.evaluate({}, validator)


def probe_fixture(tmp_path):
    def saved(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        binding = diagnostic.Bindings()
        binding.bind(path)
        return binding.files[str(path)]

    baseline = saved("baseline.json", dict(core_root="/invented/core", tool_entrypoints={
        name: dict(absolute_path="/invented/" + name) for name in ("orthohmm_python", "mafft", "FastTree")}))
    run = dict(index=0, output_root=str(tmp_path / "run"), cell="p1_c0_r1", genes=3,
               proteomes=2, inputs=[], native_order=["beta.fa", "alpha.fa"])
    plan = saved("plan.json", dict(runs=[run], baseline=baseline))
    amendment = saved("amendment.json", dict(historical_plan=plan))
    request = saved("request.json", dict(index=0, job_id=999, amendment=amendment, plan=plan))
    failure = saved("failure.json", dict(request=request))
    diagnosis = saved("diagnosis.json", dict(schema="native_partition_diagnosis_v1",
        status="diagnosis_completed", index=0, job_id=999, plan=plan,
        original_failure=failure, evidence=[], terminal_reviewed=False, native_outputs_validated=False))
    return saved, run, diagnosis


def test_probe_reconstructs_original_caller_context_without_scheduler_or_resource_work(tmp_path, monkeypatch):
    from benchmark_tools import validate_native_factorial_outputs as kernel
    from benchmark_tools.native_factorial_allocated_execution import native_command

    _, run, diagnosis = probe_fixture(tmp_path)
    observed = []
    monkeypatch.setattr(kernel, "validate_semantics", lambda context: observed.append(context) or {})
    result = diagnostic.probe(diagnosis)
    context = observed[0]
    assert context["cpu"] == 32
    assert context["threads_per_worker"] == 4
    assert context["input_directory"] == str(tmp_path / "run/input")
    assert context["cwd"] == "/invented/core"
    assert context["aligner"] == "/invented/mafft"
    assert context["tree_builder"] == "/invented/FastTree"
    assert context["cell"] == run["cell"]
    assert context["native_order"] == run["native_order"]
    binding = diagnostic.Bindings()
    original = binding.read(result["request"]["path"])
    baseline = binding.read(str(tmp_path / "baseline.json"))
    assert context["command"] == native_command(original["amendment"], run, baseline, metrics=True)
    assert result["status"] == "semantic_probe_passed"
    assert result["native_outputs_validated"] is False


@pytest.mark.parametrize("key,value", [("index", 1), ("job_id", 1000), ("plan", {})])
def test_probe_rejects_wrong_original_request_before_kernel(tmp_path, monkeypatch, key, value):
    from benchmark_tools import validate_native_factorial_outputs as kernel

    saved, _, diagnosis = probe_fixture(tmp_path)
    request = json.loads((tmp_path / "request.json").read_text())
    request[key] = value
    request = saved("request.json", request)
    failure = saved("failure.json", dict(request=request))
    data = json.loads((tmp_path / "diagnosis.json").read_text())
    data["original_failure"] = failure
    diagnosis = saved("diagnosis.json", data)
    calls = []
    monkeypatch.setattr(kernel, "validate_semantics", lambda context: calls.append(context))
    with pytest.raises(ValueError):
        diagnostic.probe(diagnosis)
    assert calls == []
