import numpy as np

from benchmark_tools.hmm_generalization_holdout import (
    SyntheticQuery,
    classification_metrics,
    generate_synthetic_dataset,
    select_profile,
    species_sequences,
)


def test_synthetic_dataset_is_deterministic_and_balanced():
    left = generate_synthetic_dataset(17, families=4)
    right = generate_synthetic_dataset(17, families=4)

    assert left == right
    assert len(left.training_sequences) == 4
    assert all(len(members) == 5 for members in left.training_sequences.values())
    assert sum(query.family_id is not None for query in left.queries) == 16
    assert sum(query.family_id is None for query in left.queries) == 4


def test_species_sequences_preserves_ids_and_sequences():
    database = species_sequences(["a", "b"], ["ACDE", "FGHIK"], "test")

    assert database.ids == ["a", "b"]
    assert database.lengths.tolist() == [4, 5]
    assert database.get_sequence(0).shape == (4,)
    assert database.get_sequence(1).shape == (5,)


def test_margin_selection_calibrates_for_model_length():
    eligible = [
        (3, 900.0, 890.0, 500.0),
        (7, 240.0, 180.0, 200.0),
    ]

    assert select_profile(eligible, "raw_score") == 3
    assert select_profile(eligible, "margin_per_model_position") == 7
    assert select_profile([], "raw_score") is None


def test_classification_metrics_count_wrong_family_as_fp_and_fn():
    queries = [
        SyntheticQuery("a", "ACDE", 0, "standard"),
        SyntheticQuery("b", "FGHI", 1, "standard"),
        SyntheticQuery("c", "KLMN", None, "unrelated"),
    ]
    metrics = classification_metrics(queries, [0, 0, None])

    assert metrics["true_family_assignments"] == 1
    assert metrics["predictions_made"] == 2
    assert metrics["precision"] == 0.5
    assert metrics["recall"] == 0.5
    assert metrics["f_score"] == 0.5
    assert metrics["negative_rejection_rate"] == 1.0
    assert np.isclose(metrics["overall_accuracy"], 2 / 3)
