from copy import deepcopy

import pytest

from benchmark_tools.admit_ob_candidate_neighborhood import validate_plan
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS


def fixture_report():
    baseline = {"min_norm": .03, "min_margin": 1.5, "max_iterations": 2}
    report = {"status": "five_candidate_arms_prepared_unscored", "job_id": "21314",
        "accuracy_evaluated": False, "publication_ready": False,
        "arms": [{"label": label, "delta": deepcopy(delta), "applied_parameters": {**baseline, **delta},
                  "engine_calls": 1, "engine_fixed_profile_report": {"parameters": dict(baseline)}}
                 for label, delta in ARMS]}
    return report, baseline


def test_accepts_exact_prespecified_plan():
    validate_plan(*fixture_report())


@pytest.mark.parametrize("field,value", [("status", "running"), ("job_id", "21315"),
    ("accuracy_evaluated", True), ("publication_ready", True)])
def test_rejects_unadmitted_status(field, value):
    report, baseline = fixture_report()
    report[field] = value
    with pytest.raises(ValueError, match="identity"):
        validate_plan(report, baseline)


@pytest.mark.parametrize("change", ["order", "missing", "extra", "delta", "applied", "nominal", "calls"])
def test_rejects_changed_plan(change):
    report, baseline = fixture_report()
    if change == "order":
        report["arms"].reverse()
    elif change == "missing":
        report["arms"].pop()
    elif change == "extra":
        report["arms"].append(deepcopy(report["arms"][0]))
    elif change == "delta":
        report["arms"][1]["delta"] = {}
    elif change == "applied":
        report["arms"][1]["applied_parameters"]["min_norm"] = .04
    elif change == "nominal":
        report["arms"][1]["engine_fixed_profile_report"]["parameters"]["min_norm"] = .024
    else:
        report["arms"][1]["engine_calls"] = 2
    with pytest.raises(ValueError):
        validate_plan(report, baseline)
