import copy
from itertools import combinations
import json

import pytest

from benchmark_tools.compare_qfo_scored_pairs import compare
from benchmark_tools.plot_qfo_scored_pair_decomposition import BASELINE, METHODS, counts, endpoints, run
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def panel():
    keys = [*METHODS, BASELINE]
    # Fixed shared scores; different eligible sets and zero-valued scores are retained.
    scores = {key: {("a", "b"): 500000, ("a", str(i)): i * 100000} for i, key in enumerate(keys)}
    result = dict(status="corrected_all_method_scored_pair_panel", uncertainty_admitted=False,
                  publication_ready=False, checked_records=[],
                  methods=[dict(key=key) for key in keys], endpoints=[], comparisons=[])
    for metric in ("GO", "EC"):
        for key in keys:
            mean = sum(scores[key].values()) / (len(scores[key]) * 1000000)
            result["endpoints"].append(dict(metric=metric, method=key, assessed_pairs=len(scores[key]),
                rounded_mean=mean, admitted_mean=mean, historical_execution_output_bound=True))
        for left, right in combinations(keys, 2):
            result["comparisons"].append(dict(metric=metric, left=left, right=right,
                                               result=compare(scores[left], scores[right])))
    return result


def test_complete_inventory_and_exact_terms():
    rows = endpoints(panel())
    assert len(rows) == 14
    assert [(r["metric"], r["method"]) for r in rows] == [(metric, key) for metric in ("GO", "EC") for key in METHODS]
    row = rows[0]
    assert row["reference"] == BASELINE
    assert row["shared_pairs"] == 1
    assert row["method_shared_fraction"] == .5
    assert row["method_only_mean"] == 0
    assert row["difference_points"] == -35
    assert row["shared_denominator_points"] == 0
    assert row["method_only_points"] == 0
    assert row["negative_reference_only_points"] == -35
    for row in rows:
        assert row["difference_points"] == pytest.approx(sum(row[k] for k in (
            "shared_denominator_points", "method_only_points", "negative_reference_only_points")))


def test_reverse_baseline_orientation():
    value = panel()
    expected = endpoints(value)
    for row in value["comparisons"]:
        r = row["result"]
        left = {("x", "y"): 500000, ("a", row["left"]): r["left_score_sum_millionths"] - 500000}
        right = {("x", "y"): 500000, ("a", row["right"]): r["right_score_sum_millionths"] - 500000}
        row["left"], row["right"] = row["right"], row["left"]
        row["result"] = compare(right, left)
    assert endpoints(value) == expected


@pytest.mark.parametrize("case", ["missing", "duplicate", "unknown", "method", "summary", "bound", "admitted", "nan"])
def test_invalid_panel(case):
    value = panel()
    if case == "missing":
        value["comparisons"].pop()
    elif case == "duplicate":
        value["comparisons"].append(copy.deepcopy(value["comparisons"][0]))
    elif case == "unknown":
        value["comparisons"][0]["metric"] = "FAS"
    elif case == "method":
        value["methods"][0] = value["methods"][1]
    elif case == "summary":
        value["endpoints"][0]["assessed_pairs"] += 1
    elif case == "bound":
        value["endpoints"][0]["historical_execution_output_bound"] = False
    elif case == "admitted":
        value["endpoints"][0]["admitted_mean"] += .01
    else:
        value["comparisons"][0]["result"]["left_mean"] = float("nan")
    with pytest.raises(ValueError):
        endpoints(value)


@pytest.mark.parametrize("field,value", [
    ("left_pairs", True), ("shared_pairs", 3), ("left_only_pairs", 9),
    ("left_score_sum_millionths", 9000000), ("shared_left_sum_millionths", 9000000),
    ("shared_pairs_with_different_serialized_scores", 1), ("original_mean_difference", .01)])
def test_invalid_counts(field, value):
    result = compare({("a", "b"): 500000}, {("a", "b"): 500000, ("c", "d"): 0})
    result[field] = value
    with pytest.raises(ValueError):
        counts(result)


def test_empty_and_complete_intersections():
    a = {("a", "b"): 0}
    b = {("c", "d"): 1000000}
    assert counts(compare(a, b))[:3] == (1, 1, 0)
    assert counts(compare(a, a))[:3] == (1, 1, 1)
    broken = compare(a, b)
    broken["shared_conditional_mean_difference"] = 0
    with pytest.raises(ValueError):
        counts(broken)


def test_component_tampering():
    value = panel()
    value["comparisons"][-1]["result"]["original_mean_difference_components"]["left_only"] += .01
    with pytest.raises(ValueError):
        endpoints(value)


def test_export_and_refuse_overwrite(tmp_path):
    path = tmp_path / "panel.json"
    path.write_text(json.dumps(panel()))
    output = tmp_path / "out"
    receipt = run(path, record(path)["sha256"], output)
    assert receipt["contrasts"] == 14
    assert receipt["checked_comparisons"] == 56
    assert receipt["new_scoring"] is False
    assert receipt["new_uncertainty"] is False
    assert receipt["native_endpoints_replaced"] is False
    assert receipt["visual_review_complete"] is False
    assert len((output / "endpoints.tsv").read_text().splitlines()) == 15
    for suffix in ("png", "pdf", "svg"):
        assert (output / f"qfo_pair_decomposition.{suffix}").stat().st_size > 1000
    for pin in receipt["outputs"]:
        assert pin == record(pin["path"])
    with pytest.raises(FileExistsError):
        run(path, record(path)["sha256"], output)


def test_source_and_bound_input_integrity(tmp_path):
    path = tmp_path / "panel.json"
    bound = tmp_path / "bound.txt"
    bound.write_text("original")
    value = panel()
    value["checked_records"] = [record(bound)]
    path.write_text(json.dumps(value))
    output = tmp_path / "out"
    with pytest.raises(ValueError):
        run(path, "0" * 64, output)
    bound.write_text("modified")
    with pytest.raises(ValueError):
        run(path, record(path)["sha256"], output)
    assert not output.exists()
