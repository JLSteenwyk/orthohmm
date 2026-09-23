import pytest

from benchmark_tools.inventory_swiss_retained_trees import parse


VALID = "RETAINED_CASE\tA\t2\nRETAINED_NODE\tA\t[':D=Y']\nRETAINED_LEAF\tA\tx\nRETAINED_LEAF\tA\ty\n"


def test_native_preorder():
    row = parse("Darwin banner\n" + VALID, {"A"})["A"]
    assert row["explicit_D_Y_nodes"] == 1
    assert row["leaf_count"] == 2
    assert row["duplicate_leaf_labels"] == {}


def test_duplicate_labels_retained_as_ambiguity():
    row = parse(VALID.replace("\ty\n", "\tx\n"), {"A"})["A"]
    assert row["duplicate_leaf_labels"] == {"x": 2}


def test_unannotated_is_not_explicit_event():
    row = parse(VALID.replace("[':D=Y']", "[]"), {"A"})["A"]
    assert row["explicit_D_Y_nodes"] == 0
    assert row["empty_annotation_nodes"] == 1


@pytest.mark.parametrize("annotation, count", [(":GN=1:D=Y", 1), (":D=YES", 0),
    (":GN=D=Y", 0), ("D=Y", 1)])
def test_compound_nhx_tokens(annotation, count):
    row = parse(VALID.replace(":D=Y", annotation), {"A"})["A"]
    assert row["explicit_D_Y_nodes"] == count


@pytest.mark.parametrize("text", [VALID.rsplit("RETAINED_LEAF", 1)[0],
    VALID + "RETAINED_LEAF\tA\tz\n", VALID + VALID,
    VALID.replace("RETAINED_NODE\tA", "RETAINED_NODE\tB"),
    VALID.replace("[':D=Y']", "[1]"), VALID.replace("RETAINED_NODE", "RETAINED_BAD"),
    VALID.replace("\t2\n", "\t-1\n"), ""])
def test_invalid_native_records(text):
    with pytest.raises(ValueError):
        parse(text, {"A"})
