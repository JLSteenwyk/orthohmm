import pytest

from benchmark_tools import validate_species_tree_control as validator


@pytest.mark.parametrize("change", [None, "running", "job", "method", "scored", "source_changed"])
def test_control_status_and_postflight_guards(change):
    cell = {"label": "p1_c1_r1_supplied_control"}
    status = {"status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False,
              "dataset": cell["label"], "methods": {cell["label"]: {}},
              "provenance": {"job_id": "21298", "cell": cell}}
    postflight = {"status": "complete_pending_native_equivalence", "accuracy_evaluated": False,
                  "baseline_source_unchanged": True, "cell": cell}
    if change == "running":
        status["status"] = "running"
    elif change == "job":
        status["provenance"]["job_id"] = "21290"
    elif change == "method":
        status["methods"]["extra"] = {}
    elif change == "scored":
        status["accuracy_evaluated"] = True
    elif change == "source_changed":
        postflight["baseline_source_unchanged"] = False
    if change:
        with pytest.raises(ValueError):
            validator.check_status(status, postflight, cell)
    else:
        validator.check_status(status, postflight, cell)


def test_live_job_cannot_reach_native_admission(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(validator.subprocess, "check_output", lambda *a, **k:
                        "JobIDRaw|State|ExitCode|Elapsed\n21298|RUNNING|0:0|00:01:00\n")
    monkeypatch.setattr(validator, "read_frozen", lambda *a: pytest.fail("Passed live scheduler gate"))
    with pytest.raises(ValueError):
        validator.validate(tmp_path)
