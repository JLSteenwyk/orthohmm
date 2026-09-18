import pytest

from benchmark_tools.diagnose_vgnc_block_influence import incidents, metrics, deleted_scores


def test_cross_block_row_counted_once_per_deletion():
    rows = [("a","b","TP","A","A","s1","s2"), ("a","c","FP","A","C","s1","s2")]
    total, removed = incidents(rows, {"A":"A","C":"C","D":"D"})
    assert total == {"TP":1,"FP":1,"FN":0}
    assert removed["A"] == total
    assert removed["C"] == {"TP":0,"FP":1,"FN":0}
    assert removed["D"] == {"TP":0,"FP":0,"FN":0}
    assert deleted_scores(total, removed["C"])[1] == {"precision":1.,"recall":1.,"f1":1.}


def test_merged_labels_do_not_double_remove():
    _, removed = incidents([("a","b","TP","A","B","s1","s2")], {"A":"A","B":"A"})
    assert removed["A"]["TP"] == 1


def test_undefined_ratios_are_explicit():
    assert metrics(dict(TP=0,FP=0,FN=0)) == dict(precision=None,recall=None,f1=None)


@pytest.mark.parametrize("row", [("a","a","TP","A","A","s","s"), ("a","b","bad","A","A","s","s"), ("a","b","TP","X","A","s","s")])
def test_bad_rows(row):
    with pytest.raises(ValueError):
        incidents([row], {"A":"A"})


def test_conflicting_categories_rejected():
    with pytest.raises(ValueError):
        incidents([("a","b","TP","A","A","s1","s2"), ("b","a","FN","A","A","s2","s1")], {"A":"A"})


def test_negative_remaining_counts_rejected():
    with pytest.raises(ValueError):
        deleted_scores(dict(TP=1,FP=0,FN=0), dict(TP=2,FP=0,FN=0))
