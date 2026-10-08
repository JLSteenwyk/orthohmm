import pytest

from benchmark_tools import native11_composed_review_binding as current
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import admit_conversion


def fixture():
    request_ref = {"request": 1}
    request = dict(plan={"plan": 1}, amendment={"amendment": 1}, job_id=23985)
    run = dict(index=11, cell="p1_c1_r0", dataset="qfo_corrected", repeat=0)
    review = dict(schema=current.composed.SCHEMA, status="native_success", request=request_ref,
        plan=request["plan"], amendment=request["amendment"], job_id=23985, review_job_id=24034,
        index=11, cell=run["cell"], dataset=run["dataset"], repeat=0,
        scheduler_state="COMPLETED", scheduler_exit_code="0:0", execution_scope=current.SCOPE,
        resource_scopes=current.SCOPES, terminal_reviewed=True, composed_full_review_complete=True,
        native_outputs_validated=True, primary_resources_replayed=True, shared_host_resources_reviewed=True,
        current_original_inventory_equality=False, original_ordinary_full_review_success=False,
        next_identity_authorized=False, downstream_adoption_complete=False, accuracy_evaluated=False,
        scientific_timings_admitted=False, uncontended_timing=False, automatic_retry=False, publication_ready=False,
        historical_failures_retained=[23986, 24033],
        reviews={k: {} for k in ("runtime", "resources", "environment", "outputs_or_failure")})
    producer = dict(JobIDRaw="24034", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="128G")
    return review, request_ref, request, run, producer


def test_exact_new_type_can_bind_conversion_without_authorizing_history():
    values = fixture()
    assert current.admit_composed_review(*values) == "group"
    review, request_ref, request, run, producer = values
    assert review["schema"] == current.composed.SCHEMA
    assert review["next_identity_authorized"] is False
    with pytest.raises(ValueError):
        admit_conversion(review, request_ref, request, run)


@pytest.mark.parametrize("field", ["schema", "status", "request", "plan", "amendment", "job_id",
    "review_job_id", "index", "cell", "dataset", "repeat", "scheduler_state", "scheduler_exit_code",
    "execution_scope", "resource_scopes", "terminal_reviewed", "composed_full_review_complete",
    "native_outputs_validated", "primary_resources_replayed", "shared_host_resources_reviewed",
    "current_original_inventory_equality", "original_ordinary_full_review_success", "next_identity_authorized",
    "downstream_adoption_complete", "accuracy_evaluated", "scientific_timings_admitted", "uncontended_timing",
    "automatic_retry", "publication_ready", "historical_failures_retained", "reviews"])
def test_changed_identity_proof_or_admission_claim_is_rejected(field):
    values = fixture()
    review = values[0]
    review[field] = not review[field] if type(review[field]) is bool else "wrong"
    with pytest.raises(ValueError):
        current.admit_composed_review(*values)


@pytest.mark.parametrize("field", ["JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem"])
def test_pending_failed_or_wrong_resource_review_producer_cannot_admit(field):
    values = fixture()
    values[-1][field] = "wrong"
    with pytest.raises(ValueError, match="producer"):
        current.admit_composed_review(*values)


def test_bool_repeat_is_not_an_integer_repeat():
    values = fixture()
    values[0]["repeat"] = False
    with pytest.raises(ValueError):
        current.admit_composed_review(*values)


def test_runtime_component_alone_cannot_admit_as_full_review():
    values = fixture()
    values[0]["schema"] = current.composed.runtime_v2.SCHEMA
    with pytest.raises(ValueError):
        current.admit_composed_review(*values)
