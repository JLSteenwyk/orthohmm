import gzip

import pytest

from benchmark_tools import audit_native_qfo_fragment_pairs as primary
from benchmark_tools import readback_native_qfo_fragment_pairs as reader


def annotation(flag=False, incomplete=False, later=False):
    return dict(fragment_flag=flag, incomplete_sequence_features=[{"type": "NON_TER"}] if incomplete else [],
                selection_class="later_sequence_version" if later else "baseline_release")


def panel():
    annotations = dict(a=annotation(flag=True), b=annotation(), c=annotation(later=True),
                       d=annotation(incomplete=True), e=annotation(flag=True, later=True), f=None)
    before = {("F", "a", "b"): "TP", ("F", "b", "c"): "FP", ("F", "c", "d"): "TP",
              ("F", "b", "e"): "FN", ("F", "b", "f"): "TN", ("F", "a", "f"): "FP"}
    after = {("F", "a", "b"): "FN", ("F", "b", "c"): "TN", ("F", "c", "d"): "TP",
             ("F", "b", "e"): "TP", ("F", "b", "f"): "FP", ("F", "a", "f"): "TN"}
    return before, after, annotations


def test_complete_pair_views_independent_sql_and_zero_cells():
    before, after, annotations = panel()
    actual, ledger = primary.project(before, after, annotations)
    expected, sql_ledger = reader.sql_projection(before, after, annotations)
    assert actual == expected and ledger == sql_ledger
    assert len(actual) == 6
    assert all(len(row["transitions"]) == 16 for row in actual)
    assert sum(row["relations"] for row in actual if row["view"] == "historical") == len(before)
    assert ledger[0][-2:] == ("annotation_positive", "annotation_positive")
    late = next(row for row in ledger if row[1:3] == ("b", "e"))
    assert late[-2:] == ("annotation_positive", "missing_without_positive")
    positive_overrides_missing = next(row for row in ledger if row[1:3] == ("a", "f"))
    assert positive_overrides_missing[-2:] == ("annotation_positive", "annotation_positive")
    base_unflagged = next(row for row in actual if row["view"] == "baseline_only" and row["bin"] == "all_matched_unflagged")
    assert base_unflagged["relations"] == 0
    assert base_unflagged["tp_removal_fraction"] is None and base_unflagged["fp_removal_fraction"] is None
    assert any(row["transitions"]["FN->TP"] for row in actual)
    assert any(row["transitions"]["TN->FP"] for row in actual)


@pytest.mark.parametrize("a,baseline,expected", [
    (None, False, None), (annotation(), False, False), (annotation(flag=True), False, True),
    (annotation(incomplete=True), False, True), (annotation(later=True), False, False),
    (annotation(flag=True, later=True), False, True), (annotation(later=True), True, None),
    (annotation(flag=True, later=True), True, None),
])
def test_annotation_state_and_baseline_missingness(a, baseline, expected):
    assert primary.state(a, baseline) is expected
    states = reader.annotation_states(a)
    assert states[int(baseline)] == expected


@pytest.mark.parametrize("change", ["missing_pair", "extra_pair", "truth", "label", "missing_annotation",
                                    "extra_annotation", "invalid_boolean", "invalid_features", "invalid_selection"])
def test_both_projections_refuse_changed_input(change):
    before, after, annotations = panel()
    if change == "missing_pair":
        after.pop(next(iter(after)))
    elif change == "extra_pair":
        after["F", "d", "e"] = "TP"
    elif change == "truth":
        after["F", "a", "b"] = "TN"
    elif change == "label":
        before["F", "a", "b"] = "unknown"
    elif change == "missing_annotation":
        annotations.pop("a")
    elif change == "extra_annotation":
        annotations["z"] = annotation()
    elif change == "invalid_boolean":
        annotations["a"]["fragment_flag"] = 1
    elif change == "invalid_features":
        annotations["a"]["incomplete_sequence_features"] = "NON_TER"
    else:
        annotations["a"]["selection_class"] = "current"
    for fn in (primary.project, reader.sql_projection):
        with pytest.raises(ValueError):
            fn(before, after, annotations)


def test_raw_parser_canonicalizes_and_preserves_unflagged_members(tmp_path):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(reader.HEADER + "\nF\tb\ta\tTP\n")
    assert reader.parse_labels(path, {"F": ["a", "b"]}) == {("F", "a", "b"): "TP"}


@pytest.mark.parametrize("body", ["F\ta\tb\tTP\nF\tb\ta\tFN\n", "F\ta\ta\tTP\n",
                                 "F\ta\tb\twrong\n", "G\ta\tb\tTP\n", "F\ta\tTP\n"])
def test_independent_raw_parser_refuses_invalid_rows(tmp_path, body):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(reader.HEADER + "\n" + body)
    with pytest.raises(ValueError):
        reader.parse_labels(path, {"F": ["a", "b"]})


def test_wrong_raw_header_refused(tmp_path):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write("wrong\nF\ta\tb\tTP\n")
    with pytest.raises(ValueError, match="Invalid raw header"):
        reader.parse_labels(path, {"F": ["a", "b"]})


def test_raw_member_mismatch_refused(tmp_path):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(reader.HEADER + "\nF\ta\tb\tTP\n")
    with pytest.raises(ValueError, match="Wrong represented"):
        reader.parse_labels(path, {"F": ["a", "b", "c"]})


def test_existing_output_refused_before_source_reads(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        primary.audit(tmp_path / "missing", tmp_path)
