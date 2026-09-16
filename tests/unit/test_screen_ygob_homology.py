from pathlib import Path

import pytest

from benchmark_tools.screen_ygob_homology import search_command, summarize_hits


def test_command_uses_frozen_thresholds():
    command = search_command(Path("diamond"), Path("q"), Path("db"), Path("out"), 32)
    for option, value in [("--id", "30"), ("--query-cover", "50"), ("--subject-cover", "50"), ("--evalue", "1e-5"), ("--max-target-seqs", "1")]:
        assert command[command.index(option) + 1] == value
    assert "--very-sensitive" in command


def test_counts_match_proteins_and_reference_pillars(tmp_path):
    hits = tmp_path / "hits.tsv"
    hits.write_text("a\tD0G1\t40\t70\t80\t1e-20\t50\n")
    result = summarize_hits(hits, {"a": "one", "b": "two"}, {"p": ["a", "b"]}, [2], ["development"])
    assert result["matching_query_fraction"] == 0.5
    assert result["reference_pillar_fraction_with_hit"] == 1
    assert result["by_species"]["two"]["with_hit"] == 0


@pytest.mark.parametrize("row", [
    "a\tD0G1\t29\t70\t80\t1e-20\t50\n",
    "a\tD0G3\t40\t70\t80\t1e-20\t50\n",
    "unknown\tD0G1\t40\t70\t80\t1e-20\t50\n",
    "a\tD0G1\tnan\t70\t80\t1e-20\t50\n",
    "a\tD0G1\t40\t70\t80\t1e-20\t50\n" * 2,
])
def test_invalid_hits_rejected(tmp_path, row):
    hits = tmp_path / "hits.tsv"
    hits.write_text(row)
    with pytest.raises(ValueError):
        summarize_hits(hits, {"a": "one"}, {"p": ["a"]}, [2], ["development"])
