import pytest

from benchmark_tools.compare_swisstree_nhx_events import candidates


def test_identical_lines_retained_without_selection():
    result = candidates("(X_a,Y_b)[&&NHX:D=Y]\n(X_a,Y_b)[&&NHX:D=Y]", {"a", "b"}, {("a", "b"): False})
    assert result["candidate_count"] == 2
    assert result["chosen_candidate"] is None
    assert result["identical_reference_observations"] is True
    assert result["candidates"][0]["appended_terminator_for_parse"] is True
    assert result["candidates"][0]["pair_counts"] == {"duplication_reference_nonortholog": 1}


def test_conflicting_candidates_preserved():
    result = candidates("(a,b)[&&NHX:D=Y];\n(a,b)[&&NHX:D=N];", {"a", "b"}, {("a", "b"): True})
    assert result["identical_reference_observations"] is False
    assert result["chosen_candidate"] is None
    assert len(result["candidates"][0]["disagreements"]) == 1
    assert not result["candidates"][1]["disagreements"]


@pytest.mark.parametrize("tag,label", [("", "unknown"), (":D=?", "unknown"),
    (":D=T", "duplication"), (":D=F", "speciation")])
def test_tag_semantics(tag, label):
    result = candidates(f"(a,b)[&&NHX{tag}];", {"a", "b"}, {("a", "b"): False})
    assert result["candidates"][0]["pair_counts"] == {label + "_reference_nonortholog": 1}


@pytest.mark.parametrize("text", ["", "(a,b)[&&NHX:D=INVALID];", "(a,b)[&&NHX:Ev=1>0>0];"])
def test_unsupported_semantics_rejected(text):
    with pytest.raises(ValueError):
        candidates(text, {"a", "b"}, {("a", "b"): False})
