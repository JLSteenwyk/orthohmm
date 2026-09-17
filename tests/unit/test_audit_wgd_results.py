import pytest

from benchmark_tools.audit_wgd_results import recalculate, verify_entry


def inputs():
    owners = {"a": "Scerevisiae", "b": "Scerevisiae", "c": "other", "d": "other", "x": "other", "u": "other"}
    reference = {"p": ["a", "b", "c", "d"], "q": ["x"]}
    cohort = [{"orf_pair": ["a", "b"], "split_eligible": True, "reference_eligible": True,
               "reference_pillar": "p", "experimental_class": "High"}]
    return cohort, reference, owners


def test_separated_partial_coverage_foreign_unmapped_and_fragmentation():
    cohort, refs, owners = inputs()
    row = recalculate(cohort, {"one": ["a", "c", "x", "u"], "two": ["b"], "three": ["d"]}, refs, owners)[0]
    assert row["anchor_groups"] == ["one", "two"]
    assert row["anchor_group_sizes"] == [4, 1]
    assert row["separation_rate"] == 1 and row["supported_separation_rate"] == 0
    assert row["homolog_support_by_anchor"] == [1, 0]
    assert row["coverage_numerator"] == 1 and row["coverage_denominator"] == 2
    assert row["mean_non_scer_coverage"] == .5
    assert row["foreign_pillar_members"] == ["x"] and row["unmapped_members"] == ["u"]
    assert row["pillar_native_group_count"] == 3 and row["unassigned_pillar_members"] == []


@pytest.mark.parametrize("groups,state,coverage,supported", [
    ({"one": ["a", "b", "c", "d"]}, "merged", 1, 0),
    ({"one": ["a", "c"], "two": ["b", "d"]}, "separated", 1, 1),
    ({"one": ["a", "c"]}, "incomplete_assignment", .5, 0),
    ({"one": ["c", "d"]}, "incomplete_assignment", 0, 0),
])
def test_assignment_states_do_not_invent_missing_singletons(groups, state, coverage, supported):
    cohort, refs, owners = inputs()
    row = recalculate(cohort, groups, refs, owners)[0]
    assert row["assignment_state"] == state
    assert row["mean_non_scer_coverage"] == coverage
    assert row["supported_separation_rate"] == supported


def test_excluded_pair_has_no_fabricated_endpoints():
    cohort, refs, owners = inputs()
    cohort[0].update(split_eligible=False, reference_eligible=False)
    row = recalculate(cohort, {"one": ["a", "c"]}, refs, owners)[0]
    assert row["assignment_state"] == "input_excluded"
    assert row["separation_rate"] is None and row["mean_non_scer_coverage"] is None


def test_altered_report_row_is_rejected():
    cohort, refs, owners = inputs()
    rows = recalculate(cohort, {"one": ["a", "b", "c", "d"]}, refs, owners)
    changed = {**rows[0], "supported_separation_rate": 1}
    with pytest.raises(ValueError, match="rescore differs"):
        verify_entry({"rows": [changed]}, rows)


def test_duplicate_native_member_is_rejected():
    cohort, refs, owners = inputs()
    with pytest.raises(ValueError, match="Invalid native partition"):
        recalculate(cohort, {"one": ["a"], "two": ["a"]}, refs, owners)
