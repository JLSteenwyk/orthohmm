from copy import deepcopy

import pytest

from benchmark_tools.admit_qfo_candidate_neighborhood import validate_report
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS
from benchmark_tools.prepare_qfo_candidate_neighborhood import ENVIRONMENT


def fixture():
    parameters = {"min_norm": .03, "min_margin": 1.5, "max_iterations": 2}
    baseline = {"expansion": {"parameters": parameters},
                "candidate_partition": {"sha256": "partition", "bytes": 10},
                "membership_constraints": {"sha256": "trace", "bytes": 20}}
    rows = []
    for label, delta in ARMS:
        wrapper = {"parameters": deepcopy(parameters)}
        rows.append({"label": label, "status": "candidate_prepared_unscored", "delta": dict(delta),
                     "applied_parameters": {**parameters, **delta}, "engine_calls": 1,
                     "engine_fixed_profile_report": wrapper, "incremental_seconds": 1.5,
                     "candidate_arm": {"expansion": deepcopy(wrapper),
                                       "candidate_partition": dict(baseline["candidate_partition"]),
                                       "membership_constraints": dict(baseline["membership_constraints"])}})
    rows[0]["baseline_byte_equivalent"] = True
    report = {"status": "corrected_qfo_five_candidate_neighborhood_arms_prepared_unscored",
              "job_id": "21927", "accuracy_evaluated": False, "publication_ready": False,
              "environment": dict(ENVIRONMENT), "arms": rows}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "21927"}
    return report, baseline, scheduler


def test_valid_report():
    validate_report(*fixture())


def test_running_job_rejected_before_artifact_reads(tmp_path, monkeypatch):
    from benchmark_tools.admit_qfo_candidate_neighborhood import admit
    monkeypatch.setattr("benchmark_tools.admit_qfo_candidate_neighborhood.subprocess.check_output",
        lambda *a, **k: "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21927|RUNNING|0:0|00:01:00|bizon|2\n")
    with pytest.raises(ValueError):
        admit(tmp_path, "21927", "not-yet-known", tmp_path / "out.json")
    assert not (tmp_path / "out.json").exists()


@pytest.mark.parametrize("problem", ["state", "exit", "node", "cpu", "job", "status", "scored", "ready",
    "environment", "missing", "order", "delta", "applied", "nominal", "wrapper", "calls", "bool_calls",
    "arm_state", "time_nan", "time_negative", "time_bool", "control", "partition", "trace"])
def test_rejects_incomplete_or_changed_evidence(problem):
    report, baseline, scheduler = fixture()
    row = report["arms"][1]
    if problem in ("state", "exit", "node", "cpu"):
        key = {"state": "State", "exit": "ExitCode", "node": "NodeList", "cpu": "AllocCPUS"}[problem]
        scheduler[key] = "wrong"
    elif problem == "job":
        report["job_id"] = "different"
    elif problem == "status":
        report["status"] = "preparing_unscored"
    elif problem == "scored":
        report["accuracy_evaluated"] = True
    elif problem == "ready":
        report["publication_ready"] = True
    elif problem == "environment":
        report["environment"]["OMP_NUM_THREADS"] = "2"
    elif problem == "missing":
        report["arms"].pop()
    elif problem == "order":
        report["arms"].reverse()
    elif problem == "delta":
        row["delta"] = {}
    elif problem == "applied":
        row["applied_parameters"]["min_norm"] = .9
    elif problem == "nominal":
        row["engine_fixed_profile_report"]["parameters"]["min_norm"] = .9
    elif problem == "wrapper":
        row["candidate_arm"]["expansion"]["extra"] = "changed"
    elif problem in ("calls", "bool_calls"):
        row["engine_calls"] = True if problem == "bool_calls" else 2
    elif problem == "arm_state":
        row["status"] = "preparing"
    elif problem.startswith("time_"):
        row["incremental_seconds"] = {"time_nan": float("nan"), "time_negative": -1, "time_bool": True}[problem]
    elif problem == "control":
        report["arms"][0]["baseline_byte_equivalent"] = False
    else:
        key = "candidate_partition" if problem == "partition" else "membership_constraints"
        report["arms"][0]["candidate_arm"][key]["sha256"] = "changed"
    with pytest.raises(ValueError):
        validate_report(report, baseline, scheduler)
