import pytest

from benchmark_tools.freeze_wgd_application_cohort import parse_rows, selected_examples


def row(genes=("YKL072W", "YMR053C"), counts=(10, 1, 1), category="High"):
    total = sum(counts)
    fractions = [n / total for n in counts] if max(counts) >= 6 else ["NaN"] * 3
    return ["dm", "sm1", "sm2", *genes, "NAME1", "NAME2", *counts, *fractions,
            int(max(counts) >= 6), category]


def test_reproduces_fraction_and_canonical_gene_association():
    parsed = parse_rows([row(genes=("YMR053C", "YKL072W"))])[0]
    assert parsed["orf_pair"] == ["YKL072W", "YMR053C"]
    assert parsed["gene_names"] == ["NAME2", "NAME1"]
    assert parsed["trigenic_fraction"] == 10 / 12


def test_sparse_is_missing_not_zero():
    parsed = parse_rows([row(counts=(0, 4, 1), category=None)])[0]
    assert parsed["experimental_class"] == "Sparse"
    assert parsed["trigenic_fraction"] is None


@pytest.mark.parametrize("problem", ["fraction", "class", "degree", "negative", "self_pair"])
def test_inconsistent_source_rejected(problem):
    value = row()
    if problem == "fraction":
        value[10] = .1
    elif problem == "class":
        value[14] = "Low"
    elif problem == "degree":
        value[13] = 0
    elif problem == "negative":
        value[7] = -1
    else:
        value[4] = value[3]
    with pytest.raises(ValueError):
        parse_rows([value])


def test_gene_reuse_rejected():
    with pytest.raises(ValueError, match="reused"):
        parse_rows([row(), row()])


def test_example_selection_independent_of_input_order():
    rows = []
    for i, category in enumerate(("High", "Low", "Sparse")):
        for j in range(3):
            rows.append({"orf_pair": [str(i), str(j)], "experimental_class": category,
                         "example_rank_sha256": str(9 - j)})
    expected = selected_examples(rows)
    assert len(expected) == 6
    assert selected_examples(list(reversed(rows))) == expected
