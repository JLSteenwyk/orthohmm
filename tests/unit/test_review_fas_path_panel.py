import pytest

from benchmark_tools.review_fas_path_panel import validate


def row(protein="A", paths=10**15, rejected=False):
    return dict(protein=protein, paths=paths, linearized_instances=100, rejected_by_path_limit=rejected)


def test_complete_records_strict_boundary():
    value = validate([row(), row("B", 10**15 + 1, True)], {"A", "B"}, dict(paths_limit=10**15))
    assert value == dict(proteins=2, rejected=1, maximum_paths=10**15 + 1)


@pytest.mark.parametrize("rows", [[], [row("B")], [row(), row()]])
def test_wrong_identifier_inventory_rejected(rows):
    with pytest.raises(ValueError, match="protein records"):
        validate(rows, {"A"}, dict(paths_limit=10**15))


@pytest.mark.parametrize("field,value", [("paths", -1), ("paths", True), ("paths", 1.5),
    ("linearized_instances", -1), ("rejected_by_path_limit", 0), ("rejected_by_path_limit", True)])
def test_invalid_counts_and_flags(field, value):
    record = row()
    record[field] = value
    with pytest.raises(ValueError):
        validate([record], {"A"}, dict(paths_limit=10**15))


def test_wrong_limit_rejected():
    with pytest.raises(ValueError, match="path limit"):
        validate([row()], {"A"}, dict(paths_limit=15))
