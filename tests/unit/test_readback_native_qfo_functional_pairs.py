import json
from pathlib import Path

import pytest

from benchmark_tools import compare_native_qfo_functional_pairs as primary
from benchmark_tools import readback_native_qfo_functional_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_compare_native_qfo_functional_pairs import fixture
from tests.unit.test_compare_qfo_scored_pairs import raw
from tests.unit.test_export_native_qfo_factorial_scores import write

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def test_independent_sql_join_preserves_integer_scores_and_direction(tmp_path):
    left = raw(tmp_path / "left.gz", "B\tA\t0.500000\nA\tC\t0.900000\n")
    right = raw(tmp_path / "right.gz", "A\tB\t0.400000\nB\tD\t0.100000\n")
    result = module.joined(left, right, "GO")
    assert result == dict(left_pairs=2, right_pairs=2, left_sum=1400000, right_sum=500000,
        shared_pairs=1, differing_shared_pairs=1, maximum_shared_difference=100000,
        shared_left_sum=500000, shared_right_sum=400000)


def test_disjoint_sql_join_keeps_shared_count_and_sums_zero(tmp_path):
    left = raw(tmp_path / "left.gz", "B\tA\t0.500000\n")
    right = raw(tmp_path / "right.gz", "C\tD\t0.400000\n")
    result = module.joined(left, right, "GO")
    assert result["shared_pairs"] == result["shared_left_sum"] == result["shared_right_sum"] == 0
    assert result["maximum_shared_difference"] is None


def test_fixture_sql_readback_validates_six_raw_tables(tmp_path, monkeypatch):
    _, snapshot = fixture(tmp_path, monkeypatch)
    path = tmp_path / "composition.json"
    primary.run(Path(snapshot["path"]), snapshot["sha256"], path)
    result = module.readback(path, record(path)["sha256"])
    assert result["checked_raw_rows"] == 12 and len(result["comparisons"]) == 3
    assert result["independent_parsing_and_join"] is True
    assert result["uncertainty_admitted"] is result["new_scoring_or_admission"] is result["publication_ready"] is False


@pytest.mark.parametrize("change", ["source", "direction", "missing_endpoint", "raw_hash", "sum",
                                   "fas_shared", "endpoint_mean", "publication_ready", "uncertainty_admitted"])
def test_readback_refuses_contradictions_in_primary_report(tmp_path, monkeypatch, change):
    _, snapshot = fixture(tmp_path, monkeypatch)
    path = tmp_path / "composition.json"
    composition = primary.run(Path(snapshot["path"]), snapshot["sha256"], path)
    if change == "source": composition["source"]["sha256"] = "0" * 64
    elif change == "direction": composition["comparisons"][0]["left"] = "p0_c0_r0"
    elif change == "missing_endpoint": composition["endpoints"].pop()
    elif change == "raw_hash": composition["endpoints"][0]["raw"]["sha256"] = "0" * 64
    elif change == "sum": composition["comparisons"][0]["result"]["left_score_sum_millionths"] += 1
    elif change == "fas_shared": composition["comparisons"][2]["result"]["shared_sample_pairs"] -= 1
    elif change == "endpoint_mean": composition["endpoints"][0]["raw_mean"] += .01
    else: composition[change] = True
    ref = write(tmp_path / "changed.json", composition)
    with pytest.raises(ValueError): module.readback(Path(ref["path"]), ref["sha256"])


def test_actual_native_composition_keeps_subsets_fas_sampling_and_no_intervals():
    path = RESULTS / "native_qfo_functional_pair_composition_20261006_v1.json"
    result = json.loads(path.read_text())
    assert record(path)["sha256"] == "1807b3269518590eb339637cf618f07c99abdd7e547e268fc0fea932b4ddf853"
    for row, left, right in zip(result["comparisons"][:2], (78607, 116929), (145142, 186098)):
        point = row["result"]
        assert point["left_pairs"] == point["shared_pairs"] == left
        assert point["right_pairs"] == right and point["left_only_pairs"] == 0
        assert point["shared_pairs_with_different_serialized_scores"] == 0
        assert point["shared_conditional_mean_difference"] == 0
        assert sum(point["original_mean_difference_components"].values()) == pytest.approx(point["original_mean_difference"])
    fas = result["comparisons"][2]["result"]
    assert fas["left_sample_pairs"] == 252451 and fas["right_sample_pairs"] == 38205
    assert fas["shared_sample_pairs"] == 1007
    assert fas["shared_pairs_with_different_serialized_scores"] == 0
    assert fas["shared_fraction_of_left"] < .004 and fas["shared_fraction_of_right"] < .027
    assert sum(fas["original_sample_mean_difference_components"].values()) == pytest.approx(fas["original_sample_mean_difference"])
    assert result["uncertainty_admitted"] is result["scientific_timings_admitted"] is result["publication_ready"] is False


def test_actual_sql_receipt_source_and_raw_totals_match_retained_result():
    primary_path = RESULTS / "native_qfo_functional_pair_composition_20261006_v1.json"
    receipt = json.loads((RESULTS / "native_qfo_functional_pair_sql_readback_20261006.json").read_text())
    assert receipt["composition"]["sha256"] == record(primary_path)["sha256"]
    assert receipt["source"]["sha256"] == record(module.__file__)["sha256"]
    assert receipt["checked_raw_rows"] == 817432
    assert [r["independent_sql_join"]["shared_pairs"] for r in receipt["comparisons"]] == [78607, 116929, 1007]
    assert all(r["independent_sql_join"]["differing_shared_pairs"] == 0 for r in receipt["comparisons"])
    assert receipt["independent_parsing_and_join"] is True
    assert receipt["new_scoring_or_admission"] is receipt["uncertainty_admitted"] is receipt["publication_ready"] is False


def test_manuscript_and_claims_preserve_composition_not_uncertainty():
    manuscript = (RESULTS / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    section = manuscript.split("### Native Functional-Pair Composition", 1)[1].split(
        "The original-release QfO factorial", 1)[0]
    for text in ("78,607", "116,929", "145,142", "186,098", "66,535", "69,169",
                 "0.4506870677", "0.8719528955", "1,007", "0.3989%", "2.6358%",
                 "817,432", "arithmetic composition", "no paired functional-score interval",
                 "not native inference", "unknown contention effects"):
        assert text in section
    claims = (RESULTS / "PUBLICATION_CLAIMS_20260916.md").read_text()
    assert "no functional-score paired CI or independent" in claims
    assert "Conditioning\non them changes the endpoint" in claims
    assert "| The package is publication-ready | All sections below | Not achieved |" in claims
    guide = (ROOT / "benchmark_tools/PUBLICATION_REPRODUCTION.md").read_text()
    assert "benchmark_tools.compare_native_qfo_functional_pairs" in guide
    assert "benchmark_tools.readback_native_qfo_functional_pairs" in guide
    assert "six\nsmall raw GO/EC/FAS tables" in guide
    for name in ("NATIVE_QFO_FUNCTIONAL_PAIR_RESULT_20261006.md",
                 "native_qfo_functional_pair_composition_20261006_v1.json",
                 "native_qfo_functional_pair_sql_readback_20261006.json"):
        assert "(" + name + ")" in section
        assert (RESULTS / name).is_file()
