import json
import pickle

import pytest

from benchmark_tools.orthobench_stage_diagnostics import (
    cluster_reference_pairs,
    load_refogs,
    main,
    parse_labeled_path,
    read_cached_hits,
    read_clusters,
    relation_metrics,
)


def _write_refogs(tmp_path):
    refogs = tmp_path / "RefOGs"
    refogs.mkdir()
    (refogs / "RefOG001.txt").write_text("a\nb\nc\n")
    (refogs / "RefOG002.txt").write_text("d\ne\n")
    return refogs


def test_relation_metrics_distinguish_direct_transitive_and_cross_pairs():
    metrics = relation_metrics(
        [{"a", "b", "c"}, {"d", "e"}],
        [("a", "b"), ("b", "c"), ("a", "d"), ("a", "b")],
    )

    assert metrics["reference_pairs"] == 4
    assert metrics["direct_reference_pairs"] == 2
    assert metrics["direct_reference_pair_recall_micro"] == pytest.approx(0.5)
    assert metrics["transitive_reference_pairs"] == 3
    assert metrics["transitive_reference_pair_recall_micro"] == pytest.approx(0.75)
    assert metrics["fully_connected_refogs"] == 1
    assert metrics["cross_refog_reference_pairs"] == 1
    assert metrics["relation_records"] == 4


def test_cluster_pairs_leave_missing_reference_genes_isolated(tmp_path):
    clusters_path = tmp_path / "clusters.txt"
    clusters_path.write_text("OG0: a b\nOG1: d x\n")
    clusters = read_clusters(clusters_path)

    metrics = relation_metrics(
        [{"a", "b", "c"}, {"d", "e"}],
        cluster_reference_pairs(clusters, {"a", "b", "c", "d", "e"}),
    )

    assert metrics["transitive_reference_pairs"] == 1
    assert metrics["fully_connected_refogs"] == 0


def test_cached_hits_validate_expected_payload(tmp_path):
    valid = tmp_path / "hits.pkl"
    valid.write_bytes(pickle.dumps({"all_hits": {("a", "b"): 1.0}}))
    assert list(read_cached_hits(valid)) == [("a", "b")]

    invalid = tmp_path / "invalid.pkl"
    invalid.write_bytes(pickle.dumps({"other": {}}))
    with pytest.raises(ValueError, match="normalized-hit cache"):
        list(read_cached_hits(invalid))


def test_load_refogs_preserves_overlapping_reference_membership(tmp_path):
    refogs = _write_refogs(tmp_path)
    (refogs / "RefOG002.txt").write_text("a\ne\n")
    groups = load_refogs(refogs)

    assert groups == [{"a", "b", "c"}, {"a", "e"}]
    metrics = relation_metrics(groups, [("a", "b"), ("a", "e")])
    assert metrics["direct_reference_pairs"] == 2


def test_main_records_provenance_and_metrics(tmp_path, monkeypatch):
    refogs = _write_refogs(tmp_path)
    pairs = tmp_path / "edges.abc"
    pairs.write_text("a b 1.0\nb c 0.8\na d 0.1\n")
    clusters = tmp_path / "clusters.txt"
    clusters.write_text("a b c\nd e\n")
    output = tmp_path / "diagnostics.json"
    monkeypatch.setattr("sys.argv", ["orthobench_stage_diagnostics.py"])

    assert main([
        "--refogs", str(refogs),
        "--pairs", f"edges={pairs}",
        "--clusters", f"clusters={clusters}",
        "--json", str(output),
    ]) == 0

    payload = json.loads(output.read_text())
    assert payload["schema_version"] == 1
    assert payload["source"]["sha256"]
    assert payload["refogs"]["groups"] == 2
    assert payload["stages"][0]["input"]["sha256"]
    assert payload["stages"][0]["metrics"]["fully_connected_refogs"] == 1
    assert payload["stages"][1]["metrics"]["fully_connected_refogs"] == 2


def test_main_selects_predefined_partition(tmp_path, monkeypatch):
    refogs = _write_refogs(tmp_path)
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("a\tb\n")
    partition = tmp_path / "partition.json"
    partition.write_text(json.dumps({
        "development": ["RefOG001.txt"],
        "validation": ["RefOG002.txt"],
    }))
    output = tmp_path / "diagnostics.json"
    monkeypatch.setattr("sys.argv", ["orthobench_stage_diagnostics.py"])

    assert main([
        "--refogs", str(refogs),
        "--pairs", f"pairs={pairs}",
        "--partition-json", str(partition),
        "--partition", "development",
        "--json", str(output),
    ]) == 0

    payload = json.loads(output.read_text())
    assert payload["refogs"]["groups"] == 1
    assert payload["refogs"]["partition"] == "development"
    assert payload["refogs"]["partition_manifest"]["sha256"]
    assert payload["stages"][0]["metrics"]["reference_pairs"] == 3


def test_parse_labeled_path_requires_existing_file(tmp_path):
    with pytest.raises(Exception):
        parse_labeled_path("missing")
    with pytest.raises(Exception):
        parse_labeled_path(f"label={tmp_path / 'missing'}")
