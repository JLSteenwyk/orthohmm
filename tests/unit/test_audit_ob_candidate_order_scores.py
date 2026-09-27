import pytest

from benchmark_tools.audit_ob_candidate_order_scores import LABELS, compare_all, validate_summary


def test_all_ten_contrasts_ignore_group_labels():
    groups = {k:[frozenset("ab"),frozenset("c")] for k in LABELS}
    groups[LABELS[-1]].reverse()
    result = compare_all(groups)
    assert len(result) == 10
    assert all(r["label_invariant_equal"] for r in result.values())
    groups[LABELS[1]] = [frozenset("a"),frozenset("bc")]
    assert sum(not r["label_invariant_equal"] for r in compare_all(groups).values()) == 4


@pytest.mark.parametrize("keys", [LABELS[:-1], LABELS[::-1], (*LABELS,"extra")])
def test_wrong_arm_inventory(keys):
    with pytest.raises(ValueError):
        compare_all({k:[] for k in keys})


def summary():
    return dict(status="candidate_score_order_factorial_complete",plan={"p":1},source={"s":2},
                accuracy_evaluated=False,phylogeny_run=False,genes=251378,nonself_hits=18235373,
                removed_self_hits=251135,rows=[dict(label=k) for k in LABELS])


def test_summary_valid():
    validate_summary(summary(), {"p":1}, {"s":2})


@pytest.mark.parametrize("key,value", [("status","incomplete"),("plan",{}),("source",{}),
    ("accuracy_evaluated",True),("phylogeny_run",True),("genes",1),("nonself_hits",1),
    ("removed_self_hits",0),("rows",[])])
def test_wrong_summary_rejected(key,value):
    r = summary()
    r[key] = value
    with pytest.raises(ValueError):
        validate_summary(r, {"p":1}, {"s":2})
