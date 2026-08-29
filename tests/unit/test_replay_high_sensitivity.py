import pytest

from benchmark_tools import replay_high_sensitivity


def test_replay_profile_iterations_default_to_one():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
    ])

    assert args.profile_iterations == 1
    assert args.jackknife_profile_thresholds is False
    assert args.profile_min_species == 1
    assert args.reciprocal_profile_merges is False
    assert args.reciprocal_profile_threshold_ratio == 0.7
    assert args.reciprocal_profile_min_support == 2


def test_replay_accepts_jackknife_profile_thresholds():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--jackknife-profile-thresholds",
    ])

    assert args.jackknife_profile_thresholds is True


def test_replay_accepts_profile_min_species():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--profile-min-species", "3",
    ])

    assert args.profile_min_species == 3


def test_replay_accepts_reciprocal_profile_merges():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--reciprocal-profile-merges",
        "--reciprocal-profile-threshold-ratio", "0.65",
        "--reciprocal-profile-min-support", "3",
    ])

    assert args.reciprocal_profile_merges is True
    assert args.reciprocal_profile_threshold_ratio == 0.65
    assert args.reciprocal_profile_min_support == 3


def test_second_replay_profile_iteration_requires_fasta_directory():
    with pytest.raises(SystemExit, match="requires --fasta-directory"):
        replay_high_sensitivity.main([
            "--hits-pickle", "hits.pkl",
            "--output-directory", "output",
            "--json", "result.json",
            "--profile-iterations", "2",
        ])
