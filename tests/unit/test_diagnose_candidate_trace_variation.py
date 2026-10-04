"""Qualified trace evidence and frozen candidate-cap sensitivity diagnostics."""

from collections import Counter
import json
import math
from pathlib import Path

import pytest

from benchmark_tools.diagnose_candidate_trace_variation import (
    FEATURES, controlled_fixture, load_engine, render, trace_comparison, validate_trace,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def row(source="s", target=None, **updates):
    target = target or ["a", "b"]
    base = {"iteration": 0, "source_genes": [source], "target_genes": target,
            "source_size": 1, "target_size": len(target), "source_cluster": 1, "target_cluster": 0,
            **{key: 1.0 for key in FEATURES}}
    base.update(updates)
    return base


def test_membership_and_order_canonicalization():
    old = [row(), row("t", source_cluster=2)]
    new = [row("t", target=["b", "a"], source_cluster=2), row(target=["b", "a"])]
    result = trace_comparison(old, new, 4)
    assert result["common_semantic_merges"] == 2
    assert result["original_only_merges_by_iteration"] == {}
    assert result["common_anchors_with_changed_sources"] == []
    assert all(x["maximum_absolute_delta"] == 0 for x in result["common_feature_differences"].values())


def test_numeric_label_changes_do_not_change_semantic_identity():
    result = trace_comparison([row()], [row(source_cluster=9, target_cluster=8)], 4)
    assert result["common_semantic_merges"] == 1
    assert result["common_source_cluster_ids_changed"] == 1
    assert result["common_target_cluster_ids_changed"] == 1


def test_near_tie_and_cap_recorded_without_causal_claim():
    old = [row("shared"), row("old", source_cluster=2)]
    new = [row("shared"), row("new", source_cluster=3, support=math.nextafter(1., math.inf))]
    result = trace_comparison(old, new, 2)
    anchor = result["common_anchors_with_changed_sources"][0]
    assert anchor["both_at_attachment_cap"]
    assert anchor["original_only_sources"] == [("old",)]
    assert anchor["native_only_sources"] == [("new",)]
    assert anchor["changed_selection_support_spread"] == math.ulp(1.)
    assert "cause" not in result


def test_equal_infinite_margins_do_not_generate_nan():
    result = trace_comparison([row(margin=math.inf)], [row(margin=math.inf)], 4)
    assert result["common_feature_differences"]["margin"]["maximum_absolute_delta"] == 0
    json.dumps(result, allow_nan=False)


def test_unequal_infinite_margins_are_explicit():
    result = trace_comparison([row(margin=math.inf)], [row()], 4)
    assert result["common_feature_differences"]["margin"]["maximum_absolute_delta"] == "positive_infinity"
    json.dumps(result, allow_nan=False)


def test_duplicate_semantic_merge_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        validate_trace([row(), row(source_cluster=9)])


@pytest.mark.parametrize("updates", [
    {"source_genes": []}, {"source_genes": ["s", "s"], "source_size": 2},
    {"source_size": 2}, {"source_cluster": True}, {"source_genes": ["a"]},
    {"target_genes": [" a", "b"]}, {"iteration": 2}, {"iteration": False},
    {"support": float("nan")}, {"support": math.inf}, {"margin": -math.inf},
    {"forward_hits": True}, {"support": -1},
])
def test_invalid_trace_fields_rejected(updates):
    with pytest.raises(ValueError):
        validate_trace([row(**updates)])


def test_empty_trace_and_no_common_merges_rejected():
    with pytest.raises(ValueError, match="Empty"):
        validate_trace([])
    with pytest.raises(ValueError, match="comparable"):
        trace_comparison([row("x")], [row("y")], 4)


def test_different_rounds_are_not_collapsed():
    old = [row(), row("x", iteration=1)]
    new = [row(), row("y", iteration=1)]
    result = trace_comparison(old, new, 4)
    assert result["original_only_merges_by_iteration"] == {"1": 1}
    assert result["native_only_merges_by_iteration"] == {"1": 1}


def actual_report():
    root = Path(__file__).resolve().parents[2]
    path = root / "benchmark_tools/results/candidate_trace_variation_20261004/diagnostic.json"
    if not path.exists():
        pytest.skip("Actual diagnostic report not yet available")
    return path, json.loads(path.read_text())


def test_actual_trace_projection_readback_and_pins():
    path, report = actual_report()
    if any(not Path(ref["path"]).is_file() for ref in report["inputs"]):
        pytest.skip("Retained local trace assets not installed; fixture unit tests remain portable")
    original_ref = next(ref for ref in report["inputs"] if "publication_ob_factorial_v1/candidates" in ref["path"])
    original = json.loads(Path(original_ref["path"]).read_text())
    for point in report["points"]:
        trace = json.loads(Path(point["trace"]["path"]).read_text())
        assert point["comparison"] == json.loads(json.dumps(trace_comparison(original, trace, 4)))
        assert point["comparison"]["common_source_cluster_ids_changed"] == 0
        assert point["comparison"]["common_target_cluster_ids_changed"] == 0
        assert point["comparison"]["common_feature_differences"]["forward_hits"]["maximum_absolute_delta"] == 0
        assert point["comparison"]["common_feature_differences"]["reverse_hits"]["maximum_absolute_delta"] == 0
        assert all(anchor["both_at_attachment_cap"] for anchor in point["comparison"]["common_anchors_with_changed_sources"])
    assert [p["comparison"]["common_semantic_merges"] for p in report["points"]] == [8428, 8439, 8426]
    assert [len(p["comparison"]["common_anchors_with_changed_sources"]) for p in report["points"]] == [2, 1, 4]
    assert [p["comparison"]["original_only_anchor_identities"] for p in report["points"]] == [2, 0, 2]
    for ref in [*report["inputs"], report["source"]]:
        assert record(ref["path"]) == ref
    assert path.with_suffix(".md").read_text() == render(report)
    assert report["native_inference_or_scoring_repeated"] is False
    assert report["frozen_method_modified"] is False
    assert report["publication_ready"] is False


def test_actual_frozen_engine_fixture_reproduces():
    _, report = actual_report()
    frozen_ref = next(ref for ref in report["inputs"] if ref["path"].endswith("/orthohmm/refinement.py"))
    if not Path(frozen_ref["path"]).is_file():
        pytest.skip("Frozen local scientific source not installed")
    assert record(frozen_ref["path"]) == frozen_ref
    fixture = controlled_fixture(load_engine(Path(frozen_ref["path"])), report["fixture"]["parameters"])
    assert fixture == report["fixture"]
    cases = fixture["cases"]
    assert [case["unattached_satellites"] for case in cases] == [[8], [0], [7]]
    assert [case["partition_equal_to_baseline"] for case in cases] == [True, False, False]
    assert all(case["merges"] == 8 and case["iterations"] == 2 for case in cases)
    assert cases[1]["input_score_maximum_change"] == 0
    assert cases[2]["input_score_maximum_change"] == math.ulp(1.)
    assert all(Counter(s["iteration"] for s in case["selections"]) == {0: 4, 1: 4} for case in cases)


def test_actual_native_interpreter_report_matches_except_runtime():
    path, report = actual_report()
    native_path = path.parent.parent / "candidate_trace_variation_native_runtime_20261004/diagnostic.json"
    if not native_path.exists():
        pytest.skip("Native-interpreter check not yet available")
    native = json.loads(native_path.read_text())
    assert [key for key in report if report[key] != native[key]] == ["diagnostic_runtime"]
    assert native["diagnostic_runtime"] == {"python": "3.10.13", "numpy": "2.2.6"}
    assert native_path.with_suffix(".md").read_text() == render(native)
