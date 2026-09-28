import pytest

from benchmark_tools.compare_selectome_trees import signature


def test_child_order_does_not_change_tree_signature():
    assert signature("((a:1,b:2)7[&&NHX:D=Y],c:3);") == signature("(c:3,(b:2,a:1)7[&&NHX:D=Y]);")


def test_topology_and_events_are_separate_from_annotations():
    a = signature("((a,b)7[&&NHX:D=Y:B=90],c);")
    b = signature("((a,b)8[&&NHX:D=Y:B=95],c);")
    assert [k for k in a if a[k] != b[k]] == ["annotations", "internal_labels"]


def test_branch_length_change_is_not_topology_change():
    a, b = signature("(a:1,b:2);"), signature("(a:1,b:3);")
    assert [k for k in a if a[k] != b[k]] == ["lengths"]


def test_topology_event_and_identity_changes_detected():
    a = signature("((a,b)[&&NHX:D=Y],c);")
    assert a["topology"] != signature("((a,c)[&&NHX:D=Y],b);")["topology"]
    assert a["duplication_events"] != signature("((a,b)[&&NHX:D=N],c);")["duplication_events"]
    assert a["leaves"] != signature("((A,b)[&&NHX:D=Y],c);")["leaves"]


def test_leaf_case_is_preserved_but_exact_duplicates_rejected():
    assert signature("(A,a);")["leaves"] == ["A", "a"]
    with pytest.raises(ValueError, match="unique"):
        signature("(a,a);")
