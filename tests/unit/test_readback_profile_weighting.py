import json

import numpy as np
import pytest

from benchmark_tools.diagnose_profile_weighting import evaluate
from benchmark_tools.hmm_generalization_holdout import generate_synthetic_dataset
from benchmark_tools.readback_profile_weighting import recount


def test_independent_reader_reproduces_tiny_real_scoring_results(tmp_path):
    result = evaluate(generate_synthetic_dataset(123, families=2), 4, 1, tmp_path)
    generated = json.loads((tmp_path / "seed_123_copies_4_generated.json").read_text())
    for row in result["rows"]:
        actual = recount(row["raw"]["path"], generated)
        assert all(row[key] == value for key, value in actual.items())


def test_independent_reader_rejects_wrong_training_minimum(tmp_path):
    raw = tmp_path / "bad.npz"
    np.savez(raw, scores=np.array([[10., 12.]]), thresholds=np.array([9.]), evalues=np.array([[0.]]))
    generated = {"clusters": [[0]], "queries": [{"family_id": 0, "category": "close"}], "names": ["a"]}
    with pytest.raises(ValueError, match="training minima"):
        recount(raw, generated)


def test_independent_reader_rejects_incomplete_ownership(tmp_path):
    with pytest.raises(ValueError, match="ownership"):
        recount(tmp_path / "not_read.npz", {"clusters": [[0], [0]], "queries": [], "names": ["a"]})
