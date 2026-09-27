from benchmark_tools.explain_ob_input_differences import differences


def test_internal_terminal_and_other_changes_are_distinct():
    result = differences(dict(a="AC*D", b="ACD**", c="ACD", d="ACD"),
                         dict(a="ACD", b="ACD", c="ACE", d="ACD"), {"family": {"a", "c"}})
    a, b, c = result["changed"]
    assert a["deletion_of_all_asterisks_explains"] and not a["terminal_asterisks_only"]
    assert b["deletion_of_all_asterisks_explains"] and b["terminal_asterisks_only"]
    assert not c["deletion_of_all_asterisks_explains"]
    assert result["changed_count"] == 3 and result["direct_reference_genes"] == 2
    assert a["original_sequence_sha256"] != a["staged_sequence_sha256"]


def test_missing_added_and_unchanged_ids():
    result = differences(dict(a="ACD", b="DFG"), dict(a="ACD", c="DFG"), {})
    assert result["only_original"] == ["b"] and result["only_staged"] == ["c"]
    assert result["changed_count"] == result["direct_reference_genes"] == 0
