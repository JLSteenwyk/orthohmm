import pytest

from benchmark_tools.audit_swiss_retained_mapping import compare, event, reconstruct


def leaf(name):
    return ("RETAINED_LEAF", name)


def node(annotation="[]"):
    return ("RETAINED_NODE", annotation)


@pytest.mark.parametrize("annotations, expected", [([], "S"), (["0.99"], "S"),
    ([":D=Y"], "D"), ([":GN=1:D=Y"], "D"), (["3[&&NHX:D=Y]"], "D"),
    ([":D=N"], "S"), ([":D=T"], "D"), ([":Ev=duplication"], "D"),
    ([":D=N:D=Y"], "S"), (["0.9", ":D=Y"], "D")])
def test_generator_event_semantics(annotations, expected):
    assert event(annotations) == expected


def test_unknown_annotation_rejected():
    with pytest.raises(ValueError):
        event(["unknown"])


def test_direct_mapping_precedes_split_alias():
    members, relations, details = reconstruct([node(), leaf("x_y"), leaf("z")],
        {"x_y": 1, "x": 7, "y": 8, "z": 2})
    assert members == {1, 2}
    assert relations == {(1, 2): "S"}
    assert details["mapped_labels"]["x_y"] == 1


def test_first_split_alias_and_absent_leaf():
    members, relations, details = reconstruct([node(), leaf("x_y"), leaf("z")], {"x": 7, "y": 8})
    assert members == {7} and relations == {}
    assert details["unmapped_leaf_occurrences"] == 1


def test_duplicate_left_minus_right_and_later_overwrite():
    rows = [node("[':D=Y']"), node(), leaf("x"), leaf("y"),
            node(), leaf("x"), leaf("z")]
    members, relations, details = reconstruct(rows, {"x": 1, "y": 2, "z": 3})
    assert members == {1, 2, 3}
    assert relations == {(1, 2): "D", (1, 3): "S", (2, 3): "D"}
    assert details["intersecting_child_mappings"] == [[1]]


def test_unconsumed_tree_rejected():
    with pytest.raises(ValueError):
        reconstruct([leaf("x"), leaf("y")], {"x": 1, "y": 2})


NATIVE = """RETAINED_CASE\tA\t2
RETAINED_NODE\tA\t[]
RETAINED_LEAF\tA\tx
RETAINED_LEAF\tA\ty
SWISS_MEMBER\tA\t1
SWISS_MEMBER\tA\t2
SWISS_RELATION\tA\t1\t2\t'S'
"""


def test_exact_relation_reconstruction():
    row = compare(NATIVE, {"x": 1, "y": 2}, {"A"})["A"]
    assert row["exact_match"] is True


def test_mapping_mismatch_remains_visible():
    row = compare(NATIVE, {"x": 1, "y": 3}, {"A"})["A"]
    assert row["exact_match"] is False
    assert row["missing_members"] == [2] and row["extra_members"] == [3]
    assert row["missing_pairs"] == row["extra_pairs"] == 1


def test_event_mismatch_remains_visible():
    row = compare(NATIVE.replace("[]", "[':D=Y']"), {"x": 1, "y": 2}, {"A"})["A"]
    assert row["exact_match"] is False and row["changed_events"] == 1


@pytest.mark.parametrize("text", [NATIVE.replace("SWISS_MEMBER\tA\t2\n", ""),
    NATIVE + "SWISS_MEMBER\tA\t2\n", NATIVE + "SWISS_RELATION\tA\t1\t2\t'S'\n",
    NATIVE.replace("\t1\t2\t'S'", "\t1\t3\t'S'"),
    NATIVE.replace("\t1\t2\t'S'", "\t2\t1\t'S'"),
    NATIVE.replace("\t'S'", "\t'unknown'")])
def test_invalid_reference_exports(text):
    with pytest.raises(ValueError):
        compare(text, {"x": 1, "y": 2}, {"A"})
