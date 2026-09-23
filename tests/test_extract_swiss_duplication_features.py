import pytest

from benchmark_tools.extract_swiss_duplication_features import bins, features, node_class, validate_native


def leaf(name):
    return ("RETAINED_LEAF", name)


def node(annotation="[]"):
    return ("RETAINED_NODE", annotation)


@pytest.mark.parametrize("annotations, expected", [([], "default_speciation_nodes"),
    (["0.99"], "default_speciation_nodes"), ([":D=N"], "explicit_speciation_nodes"),
    (["3[&&NHX:D=Y]"], "explicit_duplication_nodes"),
    ([":D=N", ":D=Y"], "explicit_speciation_nodes")])
def test_annotation_categories(annotations, expected):
    assert node_class(annotations) == expected


def test_mapped_nodes_only():
    rows = [node("[':D=Y']"), leaf("absent"), node(), leaf("x"), leaf("y")]
    row = features(rows, {"x": 1, "y": 2})
    assert row["informative_nodes"] == 1
    assert row["default_speciation_nodes"] == 1
    assert row["explicit_duplication_nodes"] == 0
    assert row["duplication_fraction"] == "0"


def test_duplicate_child_does_not_create_event():
    row = features([node("[':D=Y']"), leaf("x"), leaf("alias")], {"x": 1, "alias": 1})
    assert row["informative_nodes"] == 0 and row["duplication_fraction"] is None
    assert row["child_overlap_nodes"] == 1


def test_overlap_counts_mapped_bifurcations_not_n_minus_one():
    rows = [node("[':D=Y']"), node(), leaf("x"), leaf("y"), node(), leaf("x"), leaf("z")]
    row = features(rows, {"x": 1, "y": 2, "z": 3})
    assert row["informative_nodes"] == 3 and row["mapped_members"] == 3
    assert row["duplication_fraction"] == "1/3"


def test_exact_mapping_precedes_alias():
    row = features([node(), leaf("x_y"), leaf("x")], {"x_y": 1, "x": 2})
    assert row["mapped_members"] == 2 and row["informative_nodes"] == 1


def test_ties_not_split_and_missing_preserved():
    median, groups = bins({k: {"duplication_fraction": v} for k, v in
                         {"a": "1/3", "b": "2/6", "c": "2/3", "d": None}.items()})
    assert median == "1/3"
    assert groups["lower_duplication_fraction"] == ["a", "b"]
    assert groups["upper_duplication_fraction"] == ["c"]
    assert groups["missing_duplication_fraction"] == ["d"]


def test_even_median_is_exact():
    median, _ = bins({k: {"duplication_fraction": v} for k, v in {"a": "1/3", "b": "1/2"}.items()})
    assert median == "5/12"


def test_all_missing():
    median, groups = bins({"a": {"duplication_fraction": None}})
    assert median is None and groups["missing_duplication_fraction"] == ["a"]


@pytest.mark.parametrize("rows", [[node(), leaf("x")], [leaf("x"), leaf("x")], [("bad", "x")]])
def test_malformed_tree(rows):
    with pytest.raises(ValueError):
        features(rows, {"x": 1})


def test_native_count_crosscheck():
    row = features([node(), leaf("x"), leaf("y")], {"x": 1, "y": 2})
    validate_native("DUPLICATION_FEATURE\tA\t2\t0\t0\t1\t0\n", {"A": row})
    for text in ("", "DUPLICATION_FEATURE\tA\t2\t1\t0\t0\t0\n",
                 "DUPLICATION_FEATURE\tB\t2\t0\t0\t1\t0\n"):
        with pytest.raises(ValueError):
            validate_native(text, {"A": row})
