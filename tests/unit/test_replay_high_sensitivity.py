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
