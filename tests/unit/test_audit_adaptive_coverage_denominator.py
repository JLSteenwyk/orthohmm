import pytest

from benchmark_tools.audit_adaptive_coverage_denominator import probe_case, run


def test_exact_500_query_sample_preserves_boundary_and_cap():
    row = probe_case(500,400)
    assert row["actual_probe_calls"] == row["nominal_sample_size"] == 500
    assert row["reported_coverage"] == row["observed_fraction"] == .8
    assert row["reported_cap"] == row["observed_fraction_cap"] == 60
    assert row["legacy_formula_coverage"] == .8 and row["legacy_formula_cap"] == 60


def test_corrected_native_coverage_is_distinct_from_inflated_legacy_formula():
    row = probe_case(512,512)
    assert row["actual_probe_calls"] == 512
    assert row["reported_coverage"] == row["observed_fraction"] == 1.
    assert row["legacy_formula_coverage"] == 1.024
    assert not row["reported_coverage_exceeds_one"]


def test_corrected_candidate_tier_is_distinct_from_legacy_formula():
    row = probe_case(750,600)
    assert row["reported_coverage"] == row["observed_fraction"] == .8
    assert row["reported_cap"] == row["observed_fraction_cap"] == 60
    assert row["legacy_formula_coverage"] == 1.2 and row["legacy_formula_cap"] == 20
    assert not row["cap_would_change"]


def test_historical_reproducer_rejects_corrected_engine_without_writing(tmp_path):
    output = tmp_path / "historical_audit.json"
    with pytest.raises(ValueError, match="pinned historical engine"):
        run(output)
    assert not output.exists()


@pytest.mark.parametrize("counts", [(0,0),(500,-1),(500,501),(1.5,1)])
def test_invalid_stub_inputs_are_rejected(counts):
    with pytest.raises(ValueError):
        probe_case(*counts)
