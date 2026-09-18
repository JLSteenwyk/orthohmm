from itertools import combinations
import json
import random

import pytest

from benchmark_tools.audit_three_kingdoms_pair_counts import audit_panel, count, compare, partition, record


def write(path, groups):
    path.write_text("".join(" ".join(group) + "\n" for group in groups))
    return path


def test_split_merge_missing_outside_and_singleton(tmp_path):
    ref = write(tmp_path / "ref", [["a", "b", "c"], ["d", "e"], ["f"]])
    pred = write(tmp_path / "pred", [["a", "b", "d", "outside"], ["e"], ["f"]])
    result = count(ref, pred)
    assert (result["true_positive_gene_pairs"], result["false_positive_gene_pairs"],
            result["false_negative_gene_pairs"]) == (1, 2, 3)
    assert result["reference_genes_in_prediction"] == 5
    assert result["f_score"] == pytest.approx(2 / 7)


@pytest.mark.parametrize("text", ["a a\n", "a b\nb c\n", "OG1:\n"])
def test_invalid_partition(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path)


def test_labels_blank_lines_and_no_predictions(tmp_path):
    ref = tmp_path / "ref"
    ref.write_text("\nOG0: a b\nOG1: c\n")
    pred = write(tmp_path / "pred", [])
    result = count(ref, pred)
    assert result["false_negative_gene_pairs"] == 1
    assert result["f_score"] == 0
    assert result["predicted_orthogroups"] == 0


def test_random_partitions_against_explicit_pairs(tmp_path):
    rng = random.Random(20260918)
    for _ in range(100):
        ref_groups, pred_groups = [[] for _ in range(5)], [[] for _ in range(7)]
        for i in range(30):
            ref_groups[rng.randrange(5)].append(str(i))
            if rng.random() < 0.8:
                pred_groups[rng.randrange(7)].append(str(i))
        pred_groups[0].append("outside")
        ref = write(tmp_path / "ref", ref_groups)
        pred = write(tmp_path / "pred", pred_groups)
        truth = {tuple(sorted(pair)) for g in ref_groups for pair in combinations(g, 2)}
        positive = {tuple(sorted(pair)) for g in pred_groups
                    for pair in combinations([x for x in g if x != "outside"], 2)}
        result = count(ref, pred)
        assert result["true_positive_gene_pairs"] == len(truth & positive)
        assert result["false_positive_gene_pairs"] == len(positive - truth)
        assert result["false_negative_gene_pairs"] == len(truth - positive)


@pytest.mark.parametrize("bad", [2, 1.0, True])
def test_integer_counts_must_match_exactly(bad):
    with pytest.raises(ValueError):
        compare({"count": 1}, {"count": bad})


def test_nonfinite_score_rejected():
    with pytest.raises(ValueError):
        compare({"f_score": 1.0}, {"f_score": float("nan")})


@pytest.fixture
def panel(tmp_path):
    ref = write(tmp_path / "ref", [["a", "b"]])
    groups = write(tmp_path / "orthogroups.txt", [["a", "b"]])
    score = tmp_path / "score.txt"
    score.write_text("retained historical score evidence\n")
    payload = {"dataset": {"reference_sha256": record(ref)["sha256"]}, "methods": [
        {"key": "test", "result_dir": ".", "score": count(ref, groups), "provenance": {
            "orthogroups_sha256": record(groups)["sha256"],
            "score_sha256": record(score)["sha256"]}}]}
    path = tmp_path / "panel.json"
    path.write_text(json.dumps(payload))
    return path, ref, payload


def test_panel_verified_but_not_accuracy_admitted(tmp_path, panel):
    path, ref, _ = panel
    report = audit_panel(tmp_path, path, ref)
    assert report["method_count"] == 1
    assert report["accuracy_admitted"] is False
    assert len(report["evidence"]) == 5


@pytest.mark.parametrize("filename", ["ref", "orthogroups.txt", "score.txt"])
def test_panel_source_mutation_rejected(tmp_path, panel, filename):
    path, ref, _ = panel
    (tmp_path / filename).write_text("changed\n")
    with pytest.raises(ValueError):
        audit_panel(tmp_path, path, ref)


def test_report_count_corruption_rejected(tmp_path, panel):
    path, ref, payload = panel
    payload["methods"][0]["score"]["true_positive_gene_pairs"] = 0
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Independent count mismatch"):
        audit_panel(tmp_path, path, ref)
