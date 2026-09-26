import pytest

from benchmark_tools.prepare_matched_graph_inputs import convert


GENES = {"b": ("s2.fasta", 4), "a": ("s1.fasta", 9), "isolated": ("s2.fasta", 10)}


def test_hmm_not_normalized_twice_and_isolates_retained(tmp_path):
    path = tmp_path / "hmm.tsv"
    path.write_text("query_species\ttarget_species\tquery_id\ttarget_id\tscore\tevalue\n"
                    "s2.fasta\ts1.fasta\tb\ta\t2.5\t0\n"
                    "s1.fasta\ts1.fasta\ta\ta\t3.0\t0\n")
    result = convert([(path, None)], GENES, "hmm")
    assert result["gene_names"] == ["a", "b", "isolated"]
    assert result["gene_to_species"] == [0, 1, 1]
    assert result["hit_queries"] == [0, 1] and result["hit_targets"] == [0, 0]
    assert result["hit_scores"] == [3, 2.5]


def test_diamond_uses_raw_score_once_inclusive_cutoff(tmp_path):
    path = tmp_path / "diamond.tsv"
    path.write_text("a\tb\t9\t4\t60\t100\t1e-40\n"
                    "b\tb\t4\t4\t20\t100\t1e-39\n")
    result = convert([(path, "s2.fasta")], GENES, "diamond")
    assert result["hit_scores"] == [10]
    assert result["hit_queries"] == [0] and result["hit_targets"] == [1]


def test_duplicate_across_files(tmp_path):
    path = tmp_path / "diamond.tsv"
    path.write_text("a\tb\t9\t4\t60\t100\t0\n")
    with pytest.raises(ValueError, match="across files"):
        convert([(path, "s2.fasta")] * 2, GENES, "diamond")


def test_negative_score_and_unknown_arm(tmp_path):
    path = tmp_path / "diamond.tsv"
    path.write_text("a\tb\t9\t4\t-60\t100\t0\n")
    with pytest.raises(ValueError, match="graph score"):
        convert([(path, "s2.fasta")], GENES, "diamond")
    with pytest.raises(ValueError, match="Unknown"):
        convert([], GENES, "wrong")


def test_no_hits_retains_universe():
    result = convert([], GENES, "diamond")
    assert len(result["gene_names"]) == 3 and result["hit_scores"] == []
