import pytest

from benchmark_tools.summarize_three_kingdoms_parity import parse_score
from three_kingdoms.score_against_busco import read_ogs


def test_read_ogs_rejects_gene_in_multiple_groups(tmp_path):
    groups = tmp_path / "groups.txt"
    groups.write_text("a b\nb c\n")

    with pytest.raises(ValueError, match="occurs in multiple groups"):
        read_ogs(groups)


def test_read_ogs_rejects_duplicate_within_group(tmp_path):
    groups = tmp_path / "groups.txt"
    groups.write_text("a a\n")

    with pytest.raises(ValueError, match="Duplicate gene within"):
        read_ogs(groups)


def test_report_score_uses_exact_pair_counts(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        "reference OGs: 1\n"
        "reference genes: 3\n"
        "ref-genes in prediction: 3\n"
        "predicted OGs: 1\n"
        "TP gene pairs: 2\n"
        "FP gene pairs: 1\n"
        "FN gene pairs: 1\n"
        "precision: 0.6667\n"
        "recall: 0.6667\n"
        "F-score: 0.6667\n"
    )

    parsed = parse_score(score)

    assert parsed["precision"] == pytest.approx(2 / 3)
    assert parsed["recall"] == pytest.approx(2 / 3)
    assert parsed["f_score"] == pytest.approx(2 / 3)


def test_report_score_rejects_inconsistent_metric(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        "reference OGs: 1\n"
        "reference genes: 3\n"
        "ref-genes in prediction: 3\n"
        "predicted OGs: 1\n"
        "TP gene pairs: 2\n"
        "FP gene pairs: 1\n"
        "FN gene pairs: 1\n"
        "precision: 0.5000\n"
        "recall: 0.6667\n"
        "F-score: 0.6667\n"
    )

    with pytest.raises(ValueError, match="disagrees with TP/FP/FN"):
        parse_score(score)
