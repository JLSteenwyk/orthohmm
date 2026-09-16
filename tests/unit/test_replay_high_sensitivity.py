import numpy as np
import pytest

from benchmark_tools import replay_high_sensitivity
from benchmark_tools.replay_high_sensitivity import production_refinement_hits
from orthohmm.refinement import DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES


@pytest.mark.parametrize("count", [12, 49, 50, 51, 100])
def test_replay_matches_production_species_threshold(count):
    queries = np.array([0, 1], dtype=np.int32)
    targets = np.array([1, 0], dtype=np.int32)
    scores = np.array([5.0, 6.0])
    species = np.repeat(np.arange(count) * 3 + 100, 2)
    result = production_refinement_hits(queries, targets, scores, species)
    if count >= DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES:
        assert result == ([], [], [])
    else:
        assert all(actual is expected for actual, expected in zip(result, (queries, targets, scores)))
    np.testing.assert_array_equal(queries, [0, 1])
    np.testing.assert_array_equal(targets, [1, 0])
    np.testing.assert_array_equal(scores, [5.0, 6.0])


def test_gene_count_does_not_trigger_broad_species_branch():
    hits = np.arange(1000)
    result = production_refinement_hits(hits, hits, hits, np.zeros(1000, dtype=int))
    assert all(value is hits for value in result)

def test_replay_profile_iterations_default_to_one():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
    ])

    assert args.profile_iterations == 1
    assert args.jackknife_profile_thresholds is False
    assert args.jackknife_single_copy_profiles is False
    assert args.profile_min_species == 1


def test_replay_accepts_jackknife_profile_thresholds():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--jackknife-profile-thresholds",
    ])

    assert args.jackknife_profile_thresholds is True


def test_replay_accepts_single_copy_jackknife_thresholds():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--jackknife-single-copy-profiles",
    ])

    assert args.jackknife_single_copy_profiles is True


def test_replay_rejects_two_jackknife_modes():
    with pytest.raises(SystemExit, match="either global or single-copy"):
        replay_high_sensitivity.main([
            "--hits-pickle", "hits.pkl",
            "--output-directory", "output",
            "--json", "result.json",
            "--jackknife-profile-thresholds",
            "--jackknife-single-copy-profiles",
        ])


def test_replay_accepts_profile_min_species():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--profile-min-species", "3",
    ])

    assert args.profile_min_species == 3


def test_second_replay_profile_iteration_requires_fasta_directory():
    with pytest.raises(SystemExit, match="requires --fasta-directory"):
        replay_high_sensitivity.main([
            "--hits-pickle", "hits.pkl",
            "--output-directory", "output",
            "--json", "result.json",
            "--profile-iterations", "2",
        ])
