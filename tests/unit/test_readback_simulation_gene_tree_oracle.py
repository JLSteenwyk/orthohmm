"""Independent arithmetic and stage-boundary checks, without primary imports."""

from copy import deepcopy
import hashlib

import pytest

from benchmark_tools import readback_simulation_gene_tree_oracle as readback


def score(tp, fp, fn):
    return {"tp": tp, "fp": fp, "fn": fn, "predicted_pairs": tp + fp,
        "eligible_true_pairs": tp + fn, "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0,
        "precision": tp / (tp + fp) if tp + fp else 0, "recall": tp / (tp + fn) if tp + fn else 0}


def fixture():
    local = score(1, 0, 0)
    candidates = [{"family": "f1", "genes": 2, "ancestral_families": ["1"], "status": "unambiguous_bypass",
                   "arms": {a: dict(local) for a in readback.ARMS}},
                  {"family": "f2", "genes": 1, "ancestral_families": ["1"], "status": "unambiguous_bypass",
                   "arms": {a: score(0, 0, 0) for a in readback.ARMS}}]
    cell = {"candidates": candidates, "arms": {a: score(1, 0, 2) for a in readback.ARMS}}
    return cell, {"f1": {"a", "b"}, "f2": {"c"}}, {"a": "1", "b": "1", "c": "1"}, {("a", "b"), ("a", "c"), ("b", "c")}


def test_cross_candidate_false_negatives_are_separate_from_reconciliation():
    result = readback.decompose(*fixture())
    assert result["true_pairs_across_candidates"] == 2
    assert result["residual_by_arm"]["generating_root"]["unambiguous_bypass"]["fn"] == 0


@pytest.mark.parametrize("key,value", [("tp", -1), ("fp", True), ("fn", 0.0),
    ("predicted_pairs", 7), ("eligible_true_pairs", 8), ("f1", .8), ("recall", .3)])
def test_invalid_metrics_reject(key, value):
    value_score = score(2, 1, 1)
    value_score[key] = value
    with pytest.raises(ValueError):
        readback.validate_counts(value_score)


def test_undefined_zero_metrics_are_retained():
    readback.validate_counts(score(0, 0, 0))


def test_ineligible_change_rejects():
    cell, groups, ancestry, pairs = fixture()
    cell["candidates"][0]["arms"]["generating_root"] = score(0, 0, 1)
    with pytest.raises(ValueError, match="predictions changed"):
        readback.decompose(cell, groups, ancestry, pairs)


def test_mixed_ancestry_oracle_rejects():
    cell, groups, ancestry, pairs = fixture()
    cell["candidates"][0]["status"] = "oracle_eligible"
    cell["candidates"][0]["ancestral_families"] = ["1", "2"]
    ancestry["b"] = "2"
    with pytest.raises(ValueError, match="Mixed candidate"):
        readback.decompose(cell, groups, ancestry, pairs)


def test_global_count_mismatch_rejects():
    cell, groups, ancestry, pairs = fixture()
    cell["arms"]["generating_root"] = score(1, 1, 2)
    with pytest.raises(ValueError, match="decomposition mismatch"):
        readback.decompose(cell, groups, ancestry, pairs)


def test_local_reference_count_is_recomputed():
    cell, groups, ancestry, pairs = fixture()
    cell["candidates"][0]["arms"] = {a: score(1, 0, 1) for a in readback.ARMS}
    with pytest.raises(ValueError, match="Local truth"):
        readback.decompose(cell, groups, ancestry, pairs)


def test_repeated_partition_membership_rejects():
    cell, groups, ancestry, pairs = fixture()
    groups["f2"] = {"b"}
    with pytest.raises(ValueError, match="candidate partition"):
        readback.decompose(cell, groups, ancestry, pairs)


def test_record_and_pin_detect_same_size_change(tmp_path):
    path = tmp_path / "data.json"
    path.write_text("{}")
    data, ref = readback.load(path, hashlib.sha256(b"{}").hexdigest())
    assert data == {} and ref["bytes"] == 2
    path.write_text("[]")
    with pytest.raises(ValueError, match="Changed readback"):
        readback.load(path, ref["sha256"])


def test_readback_does_not_mutate_input_records():
    original = fixture()
    before = deepcopy(original)
    readback.decompose(*original)
    assert original == before
