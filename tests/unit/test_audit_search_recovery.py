from collections import Counter
import csv
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.audit_search_recovery import audit_direction, read_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.search_decision_trace import classify_search_result


def fixture(tmp_path):
    query_ids, target_ids = ["q"], ["x", "y", "z"]
    pairs = [("q", "x"), ("q", "y"), ("q", "z")]
    watched = {pair: {"historical_present": pair[1] == "x"} for pair in pairs}
    raw = tmp_path / "raw.npz"
    result = SimpleNamespace(query_indices=np.array([0, 0], dtype=np.int32),
        target_indices=np.array([0, 1], dtype=np.int32), scores=np.array([1., 2.]),
        evalues=np.array([1e-8, 1e-4]), candidate_count=2)
    np.savez_compressed(raw, query_ids=np.array(query_ids), target_ids=np.array(target_ids),
        **{key: getattr(result, key) for key in vars(result)})
    item = {"raw": record(raw), "candidate_count": 2, "query_count": 1,
            "target_count": 3, "watched_pairs": 3, "thresholds": []}
    for threshold in (1e-4, 1e-3):
        rows = classify_search_result(result, query_ids, target_ids, pairs, threshold)
        table = tmp_path / f"{threshold}.tsv"
        with table.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, delimiter="\t", fieldnames=[
                "query", "target", "decision", "score", "evalue", "historical_present"])
            writer.writeheader()
            for row in rows:
                writer.writerow({**row, **watched[row["query"], row["target"]]})
        item["thresholds"].append({"threshold": threshold,
            "decisions": dict(Counter(row["decision"] for row in rows)),
            "historically_absent_accepted": int(threshold == 1e-3),
            "all_significant_hits": 1 if threshold == 1e-4 else 2, "table": record(table)})
    return item, query_ids, target_ids, pairs, watched


def test_readback_preserves_strict_significance_boundary(tmp_path):
    item, queries, targets, pairs, watched = fixture(tmp_path)
    result = audit_direction(item, queries, targets, pairs, watched, [1e-4, 1e-3])
    assert result[1e-4]["decisions"] == {
        "accepted": 1, "scored_not_significant": 1, "not_selected_by_prefilter": 1}
    assert result[1e-3]["decisions"] == {"accepted": 2, "not_selected_by_prefilter": 1}


@pytest.mark.parametrize("field", ["candidate_count", "query_count", "target_count", "watched_pairs"])
def test_inconsistent_report_dimensions_rejected(tmp_path, field):
    item, queries, targets, pairs, watched = fixture(tmp_path)
    item[field] += 1
    with pytest.raises(ValueError, match="dimensions"):
        audit_direction(item, queries, targets, pairs, watched, [1e-4, 1e-3])


@pytest.mark.parametrize("field", ["decisions", "historically_absent_accepted", "all_significant_hits"])
def test_counter_tampering_rejected(tmp_path, field):
    item, queries, targets, pairs, watched = fixture(tmp_path)
    item["thresholds"][0][field] = {} if field == "decisions" else 999
    with pytest.raises(ValueError, match="counts"):
        audit_direction(item, queries, targets, pairs, watched, [1e-4, 1e-3])


def test_changed_table_rejected_even_if_hash_updated(tmp_path):
    item, queries, targets, pairs, watched = fixture(tmp_path)
    entry = item["thresholds"][0]
    path = tmp_path / "0.0001.tsv"
    path.write_text(path.read_text().replace("scored_not_significant", "accepted"))
    entry["table"] = record(path)
    with pytest.raises(ValueError, match="classification"):
        audit_direction(item, queries, targets, pairs, watched, [1e-4, 1e-3])


def test_fasta_order_and_threshold_set_preserved(tmp_path):
    item, queries, targets, pairs, watched = fixture(tmp_path)
    with pytest.raises(ValueError, match="universe/order"):
        audit_direction(item, queries, targets[::-1], pairs, watched, [1e-4, 1e-3])
    with pytest.raises(ValueError, match="thresholds"):
        audit_direction(item, queries, targets, pairs, watched, [1e-3, 1e-4])


@pytest.mark.parametrize("count", [-1, 2.0, np.array([2])])
def test_invalid_count_snapshot_rejected(tmp_path, count):
    path = tmp_path / "bad.npz"
    np.savez_compressed(path, query_ids=np.array(["q"]), target_ids=np.array(["t"]),
        query_indices=np.array([0]), target_indices=np.array([0]),
        scores=np.array([1.]), evalues=np.array([0.]), candidate_count=count)
    with pytest.raises(ValueError, match="candidate count"):
        read_raw(path)
