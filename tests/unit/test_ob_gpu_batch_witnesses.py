import pytest

from benchmark_tools.audit_ob_gpu_batch_witnesses import witnesses


def test_boundary_and_missing_witness_are_not_conflated():
    rows = witnesses({("a", "b"): 1., ("b", "a"): 2.},
                     {"a": "A", "b": "B"}, {"a": 1998, "b": 1999})
    assert len(rows) == 4
    by_pair = {(r["query_species"], r["target_species"]): r for r in rows}
    assert by_pair["B", "A"]["witness"]["target_length"] == 1998
    assert by_pair["A", "B"]["retained_hits"] == 1
    assert by_pair["A", "B"]["eligible_target_hits"] == 0
    assert by_pair["A", "B"]["assessment"] == "unresolved_no_eligible_retained_target"
    assert by_pair["A", "A"]["witness"] is None


def test_witness_choice_independent_of_insertion_order():
    hits = {("a", "b"): 1., ("c", "a"): 2., ("b", "c"): 3.}
    owners = dict.fromkeys("abc", "A")
    lengths = {"a": 10, "b": 20, "c": 5}
    a = witnesses(hits, owners, lengths)
    assert a == witnesses(dict(reversed(list(hits.items()))), owners, lengths)
    assert a[0]["witness"]["target"] == "c"


@pytest.mark.parametrize("score", [0, -1, float("nan"), float("inf")])
def test_bad_score_rejected(score):
    with pytest.raises(ValueError):
        witnesses({("a", "a"): score}, {"a": "A"}, {"a": 10})


def test_unknown_hit_gene_rejected():
    with pytest.raises(ValueError):
        witnesses({("a", "b"): 1.}, {"a": "A"}, {"a": 10})
