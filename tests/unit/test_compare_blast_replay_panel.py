from collections import Counter

import pytest

from benchmark_tools.compare_blast_replay_panel import hsp_rows, messages, compare_query

ROW = b"q\ts\t100.00\t3\t0\t0\t1\t3\t1\t3\t1e-9\t42\n"


def test_hsp_comparison_preserves_duplicate_rows():
    value = hsp_rows(ROW*2, ["q", "absent"])
    assert sum(value["q"].values()) == 2
    assert value["absent"] == Counter()
    assert value["q"] != hsp_rows(ROW, ["q"])["q"]


@pytest.mark.parametrize("content", [ROW[:-1], b"q\ts\n", ROW.replace(b"q\t", b"other\t"), ROW.replace(b"42", b"4\x00")])
def test_rejects_incomplete_malformed_and_unexpected_rows(content):
    with pytest.raises(ValueError):
        hsp_rows(content, ["q"])


def test_empty_no_hit_table_is_not_fabricated_hit():
    assert hsp_rows(b"", ["q"]) == {"q": Counter()}


def diagnostic(line=1, message="SetUpBlastSearch failed."):
    return {"q": {"query_failed": True, "messages": [dict(line=line, level="ERROR", category="setup_failure", message=message)]}}


def test_diagnostics_ignore_line_number_but_preserve_text_and_multiplicity():
    assert messages(diagnostic(1), "q") == messages(diagnostic(10), "q")
    assert messages(diagnostic(), "q") != messages(diagnostic(message="different"), "q")
    repeated = diagnostic()
    repeated["q"]["messages"] *= 2
    assert messages(repeated, "q") != messages(diagnostic(), "q")


@pytest.mark.parametrize("mode,compatible,relation", [
    ("equal", True, "equal"), ("reordered", True, "equal"),
    ("lost_duplicate", False, "equal"), ("original_mismatch", False, "mismatch"),
    ("partial", True, "partial_subset"), ("bad_partial", False, "mismatch"),
    ("original_absent", True, "unobserved"), ("log_change", False, "equal")])
def test_comparison_distinguishes_complete_and_partial_queries(mode, compatible, relation):
    combined = Counter([("a",), ("b",), ("a",)])
    single, original = combined.copy(), combined.copy()
    logs = diagnostic()
    incomplete = mode in {"partial", "bad_partial"}
    if mode == "reordered":
        single = Counter([("b",), ("a",), ("a",)])
    elif mode == "lost_duplicate":
        single[("a",)] -= 1
    elif mode in {"original_mismatch", "bad_partial"}:
        original[("c",)] += 1
    elif mode == "partial":
        original = Counter([("a",)])
    elif mode == "original_absent":
        original = None
    elif mode == "log_change":
        logs = diagnostic(message="different")
    result = compare_query("q", combined, single, original, incomplete, diagnostic(), logs, diagnostic())
    assert result["diagnostic_compatible"] is compatible
    assert result["original_hsp_relation"] == relation


def test_logged_failed_query_can_reproduce_with_no_hsp_rows():
    result = compare_query("q", Counter(), Counter(), None, False, diagnostic(), diagnostic(), diagnostic())
    assert result["diagnostic_compatible"]
    assert result["combined_hsp_rows"] == 0
    assert result["original_hsp_relation"] == "unobserved"
