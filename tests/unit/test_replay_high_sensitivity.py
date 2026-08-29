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
    assert args.profile_score_per_target_residue is False
    assert args.profile_min_species == 1
    assert args.reciprocal_profile_merges is False
    assert args.reciprocal_profile_threshold_ratio == 0.7
    assert args.reciprocal_profile_min_support == 2
    assert args.profile_profile_merges is False
    assert args.profile_profile_similarity_threshold == 0.6
    assert args.profile_profile_max_combined_genes == 80
    assert args.component_split_high_duplication is False
    assert args.direct_profile_fallback is False


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


def test_replay_accepts_profile_score_per_target_residue():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--profile-score-per-target-residue",
    ])

    assert args.profile_score_per_target_residue is True


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


def test_replay_accepts_profile_profile_merges():
    args = replay_high_sensitivity.build_parser().parse_args([
        "--hits-pickle", "hits.pkl",
        "--output-directory", "output",
        "--json", "result.json",
        "--profile-profile-merges",
        "--profile-profile-similarity-threshold", "0.65",
        "--profile-profile-max-combined-genes", "60",
    ])

    assert args.profile_profile_merges is True
    assert args.profile_profile_similarity_threshold == 0.65
    assert args.profile_profile_max_combined_genes == 60


def test_second_replay_profile_iteration_requires_fasta_directory():
    with pytest.raises(SystemExit, match="requires --fasta-directory"):
        replay_high_sensitivity.main([
            "--hits-pickle", "hits.pkl",
            "--output-directory", "output",
            "--json", "result.json",
            "--profile-iterations", "2",
        ])
