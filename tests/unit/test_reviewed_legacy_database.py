import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import reviewed_legacy_database as module
from benchmark_tools.admit_blast_recovery_search import admit


def evidence():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_recovery_native_residue_deletions_20260925.json"
    review = json.loads(path.read_text())
    return dict(status="database_sequence_differences_require_review", content=copy.deepcopy(review["content"])), review


def test_exact_reviewed_scope():
    database, review = evidence()
    module.validate_content(database, review)


@pytest.mark.parametrize("key,value", [("input_sequences", 984136), ("database_sequences", 984136),
    ("exact_sequence_matches", 984129), ("input_residues", 440246927),
    ("database_residues", 440246934), ("exact_sequence_parity", True),
    ("header_order_verified", False), ("differences", [])])
def test_fresh_difference_rejected(key, value):
    database, review = evidence()
    database["content"][key] = value
    with pytest.raises(ValueError):
        module.validate_content(database, review)


def test_review_digest_is_bound():
    _, review = evidence()
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_recovery_native_residue_deletions_20260925.json"
    assert module.record(path)["sha256"] == module.REVIEW_SHA
    assert len(review["transformations"]) == 7


def test_original_mode_cannot_opt_in(tmp_path):
    with pytest.raises(ValueError, match="requires replacement"):
        admit(tmp_path, tmp_path / "new", reviewed_native_o=True)
