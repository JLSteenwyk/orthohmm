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


def test_numeric_replay_keeps_gene_indices_self_hits_and_mmaps(tmp_path):
    from orthohmm.accuracy import write_accuracy_checkpoint
    from benchmark_tools.orthobench_stage_diagnostics import file_provenance
    checkpoint = write_accuracy_checkpoint(str(tmp_path), ["z", "a"], [9, 4], [0, 1], [0, 0], [3., 2.])
    checksum = file_provenance(checkpoint / "manifest.json")["sha256"]
    names, species, queries, targets, scores, evidence = replay_high_sensitivity.load_replay_input(
        checkpoint=checkpoint, checkpoint_sha256=checksum)
    assert names == ["z", "a"]
    for array in (species, queries, targets, scores):
        assert isinstance(array, np.memmap)
        assert not array.flags.writeable
    assert species.tolist() == [9, 4]
    assert queries.tolist() == [0, 1]
    assert targets.tolist() == [0, 0]
    assert scores.tolist() == [3., 2.]
    assert evidence["summary"]["self_hits"] == 1


def test_numeric_replay_rejects_wrong_manifest_hash(tmp_path):
    from orthohmm.accuracy import write_accuracy_checkpoint
    checkpoint = write_accuracy_checkpoint(str(tmp_path), ["a"], [0], [0], [0], [1.])
    with pytest.raises(ValueError):
        replay_high_sensitivity.load_replay_input(checkpoint=checkpoint, checkpoint_sha256="0" * 64)


@pytest.mark.parametrize("flags", [[], ["--hits-pickle", "x", "--accuracy-checkpoint", "y"]])
def test_replay_requires_exactly_one_input(flags):
    with pytest.raises(SystemExit):
        replay_high_sensitivity.build_parser().parse_args([*flags, "--output-directory", "out", "--json", "out.json"])


def test_numeric_manifest_required_before_outputs(tmp_path):
    output = tmp_path / "out"
    with pytest.raises(SystemExit, match="required together"):
        replay_high_sensitivity.main(["--accuracy-checkpoint", "missing", "--output-directory", str(output), "--json", str(tmp_path / "out.json")])
    assert not output.exists()
