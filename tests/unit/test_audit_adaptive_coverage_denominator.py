import pytest

from benchmark_tools.audit_adaptive_coverage_denominator import probe_case


def test_exact_500_query_sample_preserves_boundary_and_cap():
    row = probe_case(500,400)
    assert row["actual_probe_calls"] == row["nominal_sample_size"] == 500
    assert row["legacy_reported_coverage"] == row["observed_fraction"] == .8
    assert row["legacy_cap"] == row["observed_fraction_cap"] == 60


def test_unchanged_native_loop_can_report_coverage_above_one():
    row = probe_case(512,512)
    assert row["actual_probe_calls"] == 512
    assert row["legacy_reported_coverage"] == 1.024 and row["observed_fraction"] == 1.
    assert row["reported_coverage_exceeds_one"]


def test_denominator_changes_candidate_tier_at_same_explicit_hits():
    row = probe_case(750,600)
    assert row["legacy_reported_coverage"] == 1.2 and row["observed_fraction"] == .8
    assert row["legacy_cap"] == 20 and row["observed_fraction_cap"] == 60


@pytest.mark.parametrize("counts", [(0,0),(500,-1),(500,501),(1.5,1)])
def test_invalid_stub_inputs_are_rejected(counts):
    with pytest.raises(ValueError):
        probe_case(*counts)
