import pytest

from benchmark_tools.report_ygob_validation import METHODS, assemble_report, markdown, read_checkpoint


def test_frozen_report_keeps_diagnostic_out_of_contrasts():
    references = {"P": ["a", "b"], "Q": ["c", "d"], "S": ["e"]}
    predictions = {m: references for m in METHODS}
    predictions["orthofinder_sequence_only"] = {"merged": list("abcde")}
    report = assemble_report(predictions, references, "abcde")
    uncertainty = report["uncertainty"]
    assert uncertainty["replicates"] == 20000
    assert uncertainty["seed"] == 20260917
    assert uncertainty["multiplicity_count"] == 6
    assert set(uncertainty["comparisons"]) == {"orthohmm_satellite_v2", "orthohmm_high_sensitivity"}
    assert report["scores"]["orthofinder_sequence_only"]["counts"] == {"tp": 2, "fp": 8, "fn": 0}
    for comparison in uncertainty["comparisons"].values():
        for metric in comparison.values():
            assert metric["difference_percentage_points"] == 0
            assert metric["bonferroni_ci"] == [0, 0]
    rendered = markdown(report)
    assert "33.333333" in rendered
    assert "Primary: OrthoHMM satellite_v2" in rendered
    assert "not resolved pairwise orthology" in rendered
    assert report["completion_gates_verified_by_this_module"] is False


def test_report_rejects_missing_or_extra_methods():
    for predictions in ({}, {**dict.fromkeys(METHODS, {}), "extra": {}}):
        with pytest.raises(ValueError, match="four frozen"):
            assemble_report(predictions, {"P": ["a"]}, "a")


@pytest.mark.parametrize("case", ["valid", "duplicate_original", "missing_input", "unknown", "duplicate_member", "missing_member"])
def test_checkpoint_universe_and_id_validation(tmp_path, case):
    mapping = tmp_path / "SequenceIDs.txt"
    clusters = tmp_path / "clusters.txt"
    mapping.write_text("0_0: a description\n0_1: b\n")
    clusters.write_text("(mclmatrix\nbegin\n0 0_0 0_1 $\n)\n")
    if case == "duplicate_original":
        mapping.write_text("0_0: a\n0_1: a\n")
    elif case == "missing_input":
        mapping.write_text("0_0: a\n")
    elif case == "unknown":
        clusters.write_text("begin\n0 0_0 unknown $\n)\n")
    elif case == "duplicate_member":
        clusters.write_text("begin\n0 0_0 0_1 0_0 $\n)\n")
    elif case == "missing_member":
        clusters.write_text("begin\n0 0_0 $\n)\n")
    if case == "valid":
        assert read_checkpoint(clusters, mapping, "ab") == {"MCL0": ["a", "b"]}
    else:
        with pytest.raises(ValueError):
            read_checkpoint(clusters, mapping, "ab")
