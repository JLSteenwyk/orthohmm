"""Recovered accuracy conversion must never turn measurement failure into success."""

import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import prepare_measurement_failed_native_qfo_pairs as module
from benchmark_tools import prepare_native_factorial_qfo_pairs as original
from tests.unit.test_prepare_native_factorial_qfo_pairs import (
    admission_fixture, conversion_fixture, joined_conversion,  # noqa: F401
)


def recovered_fixture(index=7):
    review, ref, request, run = admission_fixture(index)
    review.pop("reviews")
    review.update(schema="native_factorial_measurement_failure_review_v1",
        status="native_scientific_outputs_recovered_measurement_failure_retained",
        scheduler_state="FAILED", scheduler_exit_code="1:0", resources=None,
        native_command_success=True, primary_resources_replayed=False, shared_host_resources_reviewed=False)
    for key in ("scientific_timings_admitted", "eligible_for_timing_comparison", "scheduler_success",
        "original_receipts_rewritten", "inference_reexecuted", "accuracy_evaluated"):
        review[key] = False
    return review, ref, request, run


@pytest.mark.parametrize("index", range(6, 13))
def test_recovery_has_explicit_separate_admission(index):
    args = recovered_fixture(index)
    assert module.admit_recovery(*args) == ("native" if index in {7, 9, 10, 12} else "group")
    with pytest.raises(ValueError, match="successful"):
        original.admit_conversion(*args)


@pytest.mark.parametrize("key,value", [
    ("schema", "native_factorial_terminal_review_v1"), ("status", "native_success"),
    ("scheduler_state", "COMPLETED"), ("scheduler_exit_code", "0:0"),
    ("job_id", 111), ("index", 6), ("cell", "p0_c0_r0"), ("dataset", "orthobench"),
    ("repeat", 3), ("request", {}), ("plan", {}), ("terminal_reviewed", False),
    ("native_outputs_validated", False), ("native_command_success", False), ("resources", {}),
    ("primary_resources_replayed", True), ("shared_host_resources_reviewed", True),
    ("scientific_timings_admitted", True), ("eligible_for_timing_comparison", True),
    ("scheduler_success", True), ("original_receipts_rewritten", True), ("inference_reexecuted", True),
    ("automatic_retry", True), ("accuracy_evaluated", True), ("uncontended_timing", True),
    ("execution_scope", "isolated"),
])
def test_recovery_rejects_relabelling_or_changed_identity(key, value):
    review, ref, request, run = recovered_fixture()
    review[key] = value
    with pytest.raises(ValueError, match="failed timing retained"):
        module.admit_recovery(review, ref, request, run)


def test_missing_null_resource_field_is_not_accepted():
    args = recovered_fixture()
    del args[0]["resources"]
    with pytest.raises(ValueError):
        module.admit_recovery(*args)


@pytest.mark.parametrize("kind,expected", [("native", 1), ("group", 5)])
def test_reuses_original_conversion_and_all_input_coverage(conversion_fixture, tmp_path, kind, expected):
    inputs, owners, group, native, mapping = conversion_fixture
    normalized = module.normalize_owners(owners, mapping)
    a, b, count, _, coverage = module.materialize(kind, native if kind == "native" else group,
        tmp_path, inputs, owners, mapping, normalized, 1 if kind == "native" else None)
    assert count == expected and a.read_bytes() == b.read_bytes()
    assert coverage["input_accessions"] == 5
    assert coverage["accessions_in_any_pair"] == (2 if kind == "native" else 4)


def test_mapping_loss_retains_failed_partial(conversion_fixture, tmp_path):
    inputs, owners, _, native, mapping = conversion_fixture
    normalized = module.normalize_owners(owners, mapping)
    with pytest.raises(ValueError, match="mapping loses"):
        module.materialize("native", native, tmp_path, inputs, owners, {"A1": 1}, normalized, 1)
    assert (tmp_path / "pairs.partial.tsv").read_text() == "A1\tB1\n"
    assert (tmp_path / "pairs.qfo.partial.tsv").read_text() == ""


@pytest.fixture
def recovered_joined(joined_conversion, monkeypatch):
    request_ref, review_ref, destination, index = joined_conversion
    root = destination.parent
    old = json.loads(Path(review_ref["path"]).read_text())
    new, _, _, _ = recovered_fixture(index)
    shutil.copyfile(Path(module.__file__).with_name("review_native_factorial_measurement_failure.py"),
        root / "benchmark_tools/review_native_factorial_measurement_failure.py")
    shutil.copyfile(Path(module.__file__).with_name("prepare_native_factorial_qfo_pairs.py"),
        root / "benchmark_tools/prepare_native_factorial_qfo_pairs.py")
    for key in ("request", "plan"):
        new[key] = old[key]
    output_ref = old["reviews"]["outputs_or_failure"]
    outputs = module.read(output_ref)
    for key in ("request", "plan", "job_id", "index", "evidence"):
        outputs.pop(key)
    outputs.update(schema="native_factorial_output_review_v1", status="native_outputs_validated",
        source=module.record(root / "benchmark_tools/validate_native_factorial_outputs.py"), input_genes=5,
        factors=module.expected_factors(new["cell"]), accuracy_evaluated=False,
        resource_measurements_admitted=False, next_identity_authorized=False)
    Path(output_ref["path"]).write_text(json.dumps(outputs))
    new.update(source=module.record(root / "benchmark_tools/review_native_factorial_measurement_failure.py"),
        outputs=module.record(output_ref["path"]), original_failed_review=old["reviews"]["runtime"],
        cadence_diagnosis=old["reviews"]["runtime"], environment_report=old["reviews"]["runtime"],
        runtime_readback=old["reviews"]["runtime"], evidence=[])
    Path(review_ref["path"]).write_text(json.dumps(new))
    monkeypatch.setattr(module, "ROOT", root)
    monkeypatch.setattr(module, "FIXED_INPUTS", original.FIXED_INPUTS)
    monkeypatch.setattr(module, "ENV_SHA", original.ENV_SHA)
    monkeypatch.setattr(module, "validate_plan", original.validate_plan)
    monkeypatch.setattr(module, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
        verified=dict(State="FAILED", ExitCode="1:0")))
    return request_ref, module.record(review_ref["path"]), destination, index


def test_joined_recovery_converts_without_success_or_overwrite(recovered_joined):
    request, review, destination, index = recovered_joined
    result = module.read(module.prepare(request, review, destination))
    assert result["status"] == "measurement_failed_native_qfo_pairs_prepared_unscored"
    assert result["total_pairs"] == result["retained_pairs"] == (1 if index == 7 else 5)
    assert result["native_scheduler"]["verified"]["State"] == "FAILED"
    assert result["resources"] is None
    for key in ("scientific_timings_admitted", "eligible_for_timing_comparison", "original_native_scheduler_success",
        "accuracy_evaluated", "automatic_retry", "native_inference_reexecuted", "next_identity_authorized"):
        assert result[key] is False
    assert result["participant"].startswith("ohmm_qfo_recovered_native_")
    with pytest.raises(FileExistsError):
        module.prepare(request, review, destination)


def test_joined_recovery_preserves_failure_after_materialization(recovered_joined, monkeypatch):
    request, review, destination, _ = recovered_joined
    def fail(*args):
        (destination / "pairs.partial.tsv").write_text("retained partial\n")
        raise ValueError("deliberate conversion failure")
    monkeypatch.setattr(module, "materialize", fail)
    with pytest.raises(ValueError, match="deliberate"):
        module.prepare(request, review, destination)
    result = module.read(module.record(destination / "results.json"))
    assert result["status"] == "measurement_failed_native_qfo_conversion_failed_retained"
    assert result["scientific_timings_admitted"] is result["accuracy_evaluated"] is False
    assert (destination / "pairs.partial.tsv").exists()


def test_changed_recovery_binding_refused_before_destination(recovered_joined):
    request, review, destination, _ = recovered_joined
    value = module.read(review)
    Path(value["outputs"]["path"]).write_text("{}")
    with pytest.raises(ValueError):
        module.prepare(request, review, destination)
    assert not destination.exists()


@pytest.mark.parametrize("change", ["scheduler", "source", "cell", "factors", "genes", "accuracy", "conflict"])
def test_rebound_but_wrong_recovered_semantics_refused(recovered_joined, monkeypatch, change):
    request, review_ref, destination, _ = recovered_joined
    review = module.read(review_ref)
    output = module.read(review["outputs"])
    if change == "scheduler":
        monkeypatch.setattr(module, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
            verified=dict(State="COMPLETED", ExitCode="0:0")))
    elif change == "source": review["source"] = review["outputs"]
    elif change == "conflict":
        review["evidence"] = [dict(review["outputs"], sha256="0" * 64)]
    else:
        if change == "cell": output["cell"] = "p1_c0_r0"
        elif change == "factors": output["factors"] = {}
        elif change == "genes": output["input_genes"] = 1
        else: output["accuracy_evaluated"] = True
        Path(review["outputs"]["path"]).write_text(json.dumps(output))
        review["outputs"] = module.record(review["outputs"]["path"])
    Path(review_ref["path"]).write_text(json.dumps(review))
    with pytest.raises(ValueError):
        module.prepare(request, module.record(review_ref["path"]), destination)
    assert not destination.exists()
