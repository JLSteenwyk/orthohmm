from copy import deepcopy

import pytest

from benchmark_tools.measure_root_context_scaling import classify_measurement


def evidence(code=0, timed_out=False):
    wall = 85800. if timed_out else 2.
    done = dict(exit_code=code, timed_out=timed_out, started_ns=1_000_000_000,
                finished_ns=1_000_000_000+int(wall*1e9))
    measured = dict(status="command_exited_zero" if code == 0 else "command_failed",
                    native=done, native_wall_s=wall, job_id=123)
    wrapper = dict(status=measured["status"], measurement=measured, before={}, after={}, scientific_results_admitted=False)
    replay = dict(status="scaling_root_context_measurement_replayed", measured=deepcopy(measured),
        native_outcome="timed_out" if timed_out else "exited_zero" if code == 0 else "exited_nonzero",
        native_exit_code=code, native_wall_s=wall, scientific_timings_admitted=False,
        environmental_validity_established=False, native_outputs_validated=False, publication_ready=False)
    return wrapper, replay


@pytest.mark.parametrize("code,timed_out,expected", [(0,False,"exited_zero"), (7,False,"exited_nonzero"),
    (-9,False,"exited_nonzero"), (124,False,"exited_nonzero"), (124,True,"timed_out")])
def test_corroborated_native_outcomes_do_not_authorize_followup(code, timed_out, expected):
    wrapper, replay = evidence(code, timed_out)
    before = deepcopy((wrapper, replay))
    result = classify_measurement(wrapper, replay, 123)
    assert result["status"] == "native_" + expected
    assert result["corroborated_native_outcome"] == expected
    assert result["reported_native"] == wrapper["measurement"]["native"]
    for key in ("next_submission_authorized", "automatic_retry", "native_outputs_validated",
                "runtime_identity_verified", "scientific_timings_admitted", "environmental_validity_established"):
        assert result[key] is False
    assert (wrapper, replay) == before


@pytest.mark.parametrize("status", ["verified_wrapper_failed", "runtime_changed_or_unverifiable"])
def test_infrastructure_failure_takes_precedence_even_after_native_success(status):
    wrapper, replay = evidence()
    wrapper["status"] = status
    result = classify_measurement(wrapper, replay, 123)
    assert result["status"] == "infrastructure_or_provenance_failure"
    assert result["reported_native"]["exit_code"] == 0
    assert result["corroborated_native_outcome"] is None


def test_missing_replay_is_not_a_native_only_failure():
    wrapper, _ = evidence(7)
    result = classify_measurement(wrapper, None, 123)
    assert result["status"] == "measurement_evidence_incomplete"
    assert result["corroborated_native_outcome"] is None
    assert result["reported_native"]["exit_code"] == 7


@pytest.mark.parametrize("fault", ["missing_before", "missing_after", "admission", "replay_admission", "replay_status",
    "raw_mismatch", "job", "bool_job", "outcome", "code", "wall", "status", "exit_bool"])
def test_contradictions_require_pause(fault):
    wrapper, replay = evidence()
    if fault == "missing_before": del wrapper["before"]
    elif fault == "missing_after": del wrapper["after"]
    elif fault == "admission": wrapper["scientific_results_admitted"] = True
    elif fault == "replay_admission": replay["scientific_timings_admitted"] = True
    elif fault == "replay_status": replay["status"] = "other"
    elif fault == "raw_mismatch": replay["measured"]["native"]["exit_code"] = 7
    elif fault in {"job", "bool_job"}:
        wrapper["measurement"]["job_id"] = replay["measured"]["job_id"] = True if fault == "bool_job" else 124
    elif fault == "outcome": replay["native_outcome"] = "exited_nonzero"
    elif fault == "code": replay["native_exit_code"] = False
    elif fault == "wall": replay["native_wall_s"] = 3.
    elif fault == "status": wrapper["status"] = "command_failed"
    elif fault == "exit_bool":
        wrapper["measurement"]["native"]["exit_code"] = replay["measured"]["native"]["exit_code"] = False
    result = classify_measurement(wrapper, replay, 123)
    assert result["status"] == "infrastructure_or_provenance_failure"
    assert result["next_submission_authorized"] is False


@pytest.mark.parametrize("job", [True, 0, -1, "123"])
def test_expected_job_must_be_explicit_positive_integer(job):
    with pytest.raises(ValueError): classify_measurement({}, None, job)


@pytest.mark.parametrize("measurement", [None, [], "invalid"])
def test_malformed_measurement_cannot_authorize_submission(measurement):
    result = classify_measurement({"measurement": measurement}, None, 123)
    assert result["status"] == "measurement_evidence_incomplete"
    assert result["next_submission_authorized"] is False


def test_nonrecord_wrapper_rejected():
    with pytest.raises(ValueError): classify_measurement(None, None, 123)
