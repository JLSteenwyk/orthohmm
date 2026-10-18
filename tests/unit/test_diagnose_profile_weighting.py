import json

import numpy as np
import pytest

from benchmark_tools.diagnose_profile_weighting import evaluate, main, summarize, training_database
from benchmark_tools.hmm_generalization_holdout import SyntheticQuery, generate_synthetic_dataset
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from orthohmm.search.sequence_weights import PROFILE_SEQUENCE_WEIGHTING


def test_expansion_keeps_species_and_cluster_ownership_explicit():
    dataset = generate_synthetic_dataset(123, families=2)
    names, sequences, species, clusters = training_database(dataset, 4)
    assert len(names) == len(sequences) == len(species) == 18
    assert len(set(names)) == 18
    assert sorted(gene for cluster in clusters for gene in cluster) == list(range(18))
    for cluster in clusters:
        assert [species[index] for index in cluster] == [0, 1, 2, 3, 4, 0, 0, 0, 0]
        assert all(sequences[index] == sequences[cluster[0]] for index in cluster[5:])


def test_assignment_summary_counts_wrong_negative_assignments_and_margins():
    queries = [SyntheticQuery("p0", "A", 0, "close"),
               SyntheticQuery("p1", "C", 1, "close"),
               SyntheticQuery("n", "D", None, "unrelated")]
    result = summarize(np.array([[12, 20, 11], [11, 25, 22]]), np.zeros((2, 3)),
                       np.array([10, 20]), queries)
    assert result["top1_positive_family_correct"] == 2
    assert result["strict_assignment"]["true_family_assignments"] == 2
    assert result["strict_assignment"]["predictions_made"] == 3
    assert result["strict_assignment"]["negative_rejection_rate"] == 0
    assert result["own_family_threshold_margin"] == {"minimum": 2, "median": 3.5, "maximum": 5}


def test_evalue_threshold_is_strict_and_no_eligible_profile_is_unassigned():
    query = [SyntheticQuery("p", "A", 0, "close")]
    result = summarize(np.array([[12]]), np.array([[1e-4]]), np.array([10]), query)
    assert result["strict_assignment"]["predictions_made"] == 0


def test_tiny_real_diagnostic_preserves_recountable_raw_arrays(tmp_path):
    dataset = generate_synthetic_dataset(123, families=2)
    result = evaluate(dataset, 4, 1, tmp_path)
    generated = json.loads((tmp_path / "seed_123_copies_4_generated.json").read_text())
    assert [row["mode"] for row in result["rows"]] == list(PROFILE_SEQUENCE_WEIGHTING)
    for pin in [result["generated"], *[row["raw"] for row in result["rows"]]]:
        check(pin)
    for row in result["rows"]:
        with np.load(row["raw"]["path"]) as raw:
            assert raw["scores"].shape == (2, 28)
            expected_thresholds = np.array([raw["scores"][index, cluster].min()
                                            for index, cluster in enumerate(generated["clusters"])])
            np.testing.assert_array_equal(raw["thresholds"], expected_thresholds)
            expected = summarize(raw["scores"][:, 18:], raw["evalues"], raw["thresholds"], dataset.queries)
            assert all(row[key] == value for key, value in expected.items())


def test_existing_output_rejected_before_work(tmp_path):
    with pytest.raises(ValueError, match="fresh output"):
        main(["--output", str(tmp_path)])


def test_cuda_visibility_must_be_explicitly_disabled(tmp_path, monkeypatch):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    with pytest.raises(ValueError, match="CUDA_VISIBLE_DEVICES"):
        main(["--output", str(tmp_path / "new")])
    assert not (tmp_path / "new").exists()
