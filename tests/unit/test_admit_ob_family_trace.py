import pytest

from benchmark_tools.admit_ob_family_trace import STAGES, verify_pairs


def evidence():
    refs = {"one": {"a", "b"}, "two": {"a", "c"}}
    indices = {stage: {"a": "x", "b": "x", "c": "y"} for stage in STAGES}
    owners = {"a": "s1", "b": "s1", "c": "s2"}
    rows = [{"refog": name, "left": "a", "right": gene, "same_species": str(gene == "b"),
             "forward_normalized_hit": "NA", "reverse_normalized_hit": "0.5",
             **{stage: str(gene == "b") for stage in STAGES}} for name, gene in (("one", "b"), ("two", "c"))]
    return rows, refs, indices, owners


def test_pair_admission_preserves_overlap_and_absent_search():
    rows, refs, indices, owners = evidence()
    result = verify_pairs(rows, refs, indices, owners)
    assert set(result) == {"one", "two"}
    assert result["one"][0]["forward_normalized_hit"] is None
    assert result["two"][0]["reverse_normalized_hit"] == .5


@pytest.mark.parametrize("change", ["missing", "duplicate", "reversed", "unknown", "columns", "boolean", "membership", "species", "score"])
def test_invalid_pair_evidence_rejected(change):
    rows, refs, indices, owners = evidence()
    if change == "missing":
        rows.pop()
    elif change == "duplicate":
        rows.append(rows[0])
    elif change == "reversed":
        rows[0]["left"], rows[0]["right"] = "b", "a"
    elif change == "unknown":
        rows[0]["refog"] = "other"
    elif change == "columns":
        rows[0]["extra"] = "x"
    elif change == "boolean":
        rows[0]["multipass"] = "1"
    elif change == "membership":
        rows[0]["multipass"] = "False"
    elif change == "species":
        rows[0]["same_species"] = "False"
    elif change == "score":
        rows[0]["reverse_normalized_hit"] = "nan"
    with pytest.raises(ValueError):
        verify_pairs(rows, refs, indices, owners)
