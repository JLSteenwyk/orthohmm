import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_native_candidate_order_scores import prepare, run
from benchmark_tools.readback_native_candidate_factorial import collect
from tests.unit.test_probe_native_candidate_order_scores import prepared


@pytest.fixture
def executed(prepared):
    write, prep, score, output = prepared
    plan = prepare(write("preparation.json", prep), write("score.json", score), output)
    return plan, run(plan), output


def test_whole_five_arm_readback(executed):
    plan, report, _ = executed
    result = collect(plan, report)
    assert len(result["rows"]) == 5
    assert len(result["comparisons"]) == 10
    assert len(result["accepted_trace_comparisons_against_historical"]) == 4
    assert all(row["whole_partition_reconstructed"] for row in result["rows"])
    assert not result["candidate_scores_recalculated"]
    assert not result["accuracy_scored"]


@pytest.mark.parametrize("problem", ["controls", "comparisons", "labels", "count", "accuracy"])
def test_inconsistent_aggregate_rejected(executed, problem):
    plan, report_ref, _ = executed
    path = Path(report_ref["path"])
    report = json.loads(path.read_text())
    if problem == "controls":
        report["controls"]["historical"] = False
    elif problem == "comparisons":
        next(iter(report["comparisons"].values()))["partition_equal"] = False
    elif problem == "labels":
        report["rows"].reverse()
    elif problem == "count":
        report["rows"][0]["groups"] += 1
    else:
        report["accuracy_scored"] = True
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        collect(plan, record(path))


def test_failure_not_promoted_to_success(executed):
    plan, report, output = executed
    (output / "failure.json").write_text("{}")
    with pytest.raises(ValueError, match="Retained failure"):
        collect(plan, report)


def test_tampered_accepted_membership_rejected(executed):
    plan, report, output = executed
    trace_path = output / "historical_order_historical_scores_merges.json"
    trace = json.loads(trace_path.read_text())
    trace[0]["source_genes"] = ["unindexed"]
    trace_path.write_text(json.dumps(trace))
    with pytest.raises(ValueError):
        collect(plan, report)


def test_actual_readback_reproduces_and_order_only_differences():
    root = Path(__file__).resolve().parents[2]
    path = root / "benchmark_tools/results/native_candidate_factorial_readback_22429.json"
    if not path.exists():
        pytest.skip("Actual local native candidate diagnostic not installed")
    saved = json.loads(path.read_text())
    if any(not Path(ref["path"]).is_file() for ref in saved["checked_records"]):
        pytest.skip("Actual retained raw assets not installed")
    current = collect(saved["plan"], saved["report"])
    assert json.loads(json.dumps(current, allow_nan=False)) == saved
    comparisons = saved["comparisons"]
    assert comparisons["historical_order_historical_scores__vs__historical_order_fresh_scores"]["partition_equal"]
    assert comparisons["fresh_order_historical_scores__vs__fresh_order_fresh_scores"]["partition_equal"]
    assert comparisons["fresh_order_fresh_scores__vs__fresh_full_self_control"]["partition_equal"]
    assert comparisons["historical_order_historical_scores__vs__fresh_order_historical_scores"]["genes_in_changed_groups"] == 58
    features = saved["accepted_trace_comparisons_against_historical"]["fresh_order_historical_scores"]["common_feature_differences"]
    for key in ("forward_hits", "reverse_hits", "forward_maximum", "reverse_maximum", "forward_coverage", "reverse_coverage"):
        assert features[key]["changed_values"] == 0
    assert features["support"]["changed_values"] == 1725
