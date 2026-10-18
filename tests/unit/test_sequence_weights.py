import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser
from orthohmm.search import profile_expansion
from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_profile import build_msa_profile, compute_pssm, encode_alignment
from orthohmm.search.sequence_weights import (
    PROFILE_SEQUENCE_WEIGHTING, effective_cluster_count, position_based_weights,
    profile_sequence_weights, validate_profile_sequence_weighting,
)


MATRIX = get_matrix("BLOSUM62")
BACKGROUND = get_background_freqs("BLOSUM62")


def test_hand_calculated_position_weights_and_effective_counts():
    msa = encode_alignment(["AA", "AA", "CC"])
    np.testing.assert_allclose(position_based_weights(msa), [.25, .25, .5])
    np.testing.assert_allclose(profile_sequence_weights(msa, "henikoff"), [.75, .75, 1.5])
    np.testing.assert_allclose(profile_sequence_weights(msa, "henikoff_cluster"), [.5, .5, 1])
    assert effective_cluster_count(msa) == 2


def test_species_totals_are_equal_without_gene_name_heuristics():
    msa = encode_alignment(["AA", "AA", "CC"])
    labels = ["unrelated_name", "same_species", "same_species"]
    weights = profile_sequence_weights(msa, "species_henikoff_cluster", labels)
    np.testing.assert_allclose(weights, [1, 1 / 3, 2 / 3])
    assert weights[0] == pytest.approx(weights[1:].sum())
    assert weights.sum() == pytest.approx(2)


def test_single_linkage_is_transitive_and_counts_only_usable_rows():
    msa = encode_alignment(["AAAA", "AAAC", "AACC", "----", "XXXX"])
    assert effective_cluster_count(msa, .75) == 1
    assert effective_cluster_count(msa, 1) == 3


@pytest.mark.parametrize("identity", [0, -1, 1.1, np.nan, np.inf, True, None, "0.62", []])
def test_invalid_identity_thresholds_rejected(identity):
    with pytest.raises(ValueError, match="Identity"):
        effective_cluster_count(encode_alignment(["AA", "CC"]), identity)


def test_fragment_length_normalization_and_unknowns():
    msa = encode_alignment(["AA--", "AAAA", "XXXX"])
    np.testing.assert_allclose(position_based_weights(msa), [.5, .75, 0])
    assert effective_cluster_count(msa) == 1
    weights = profile_sequence_weights(msa, "henikoff_cluster")
    np.testing.assert_allclose(weights, [.4, .6, 0])


@pytest.mark.parametrize("mode", PROFILE_SEQUENCE_WEIGHTING)
def test_weighting_row_permutation_and_column_permutation(mode):
    msa = encode_alignment(["AA-C", "CAAC", "CCAA", "AA-C"])
    labels = np.array([1, 2, 3, 1])
    order = np.array([3, 1, 0, 2])
    expected = profile_sequence_weights(msa, mode, labels)
    np.testing.assert_allclose(
        profile_sequence_weights(msa[order][:, ::-1], mode, labels[order]), expected[order]
    )


@pytest.mark.parametrize("mode", PROFILE_SEQUENCE_WEIGHTING[1:])
def test_empty_and_all_unknown_alignments_return_finite_zero_weights(mode):
    for msa in (np.empty((0, 3), dtype=np.int8), np.full((3, 4), 20, dtype=np.int8)):
        weights = profile_sequence_weights(msa, mode, np.arange(len(msa)))
        assert weights.shape == (len(msa),)
        assert np.isfinite(weights).all() and not weights.any()
        assert effective_cluster_count(msa) == 0


@pytest.mark.parametrize("value", [None, True, 1, [], {}, "other"])
def test_invalid_modes_rejected_before_empty_profile(value):
    with pytest.raises(ValueError, match="profile sequence weighting"):
        validate_profile_sequence_weighting(value)
    with pytest.raises(ValueError, match="profile sequence weighting"):
        build_msa_profile([], [], MATRIX, BACKGROUND, sequence_weighting=value)


@pytest.mark.parametrize("msa", [[], np.zeros((2, 2, 2), dtype=int), [[-1]], [[21]], [[1.0]]])
def test_invalid_encoded_alignments_rejected(msa):
    with pytest.raises(ValueError, match="2D integer alignment"):
        position_based_weights(msa)


@pytest.mark.parametrize("labels", [None, [1], [[1], [2]], [1, np.nan], [None, None]])
def test_invalid_species_labels_rejected(labels):
    with pytest.raises(ValueError, match="species label"):
        profile_sequence_weights(encode_alignment(["AA", "CC"]), "species_henikoff_cluster", labels)


@pytest.mark.parametrize("weights", [[-1, 2], [0, 0], [np.nan, 1], [np.inf, 1], [1], [[1, 1]]])
def test_invalid_explicit_counts_rejected(weights):
    with pytest.raises(ValueError, match="Sequence weights"):
        compute_pssm(encode_alignment(["AA", "CC"]), MATRIX, BACKGROUND, sequence_weights=weights)


def test_uniform_counts_are_exact_default_and_weighting_preserves_column_mask():
    msa = encode_alignment(["AA-C", "CC-C", "CCAA"])
    legacy = compute_pssm(msa, MATRIX, BACKGROUND)
    uniform = compute_pssm(msa, MATRIX, BACKGROUND, sequence_weights=np.ones(3))
    weighted = compute_pssm(msa, MATRIX, BACKGROUND, sequence_weights=[.1, .1, 1.8])
    for left, right in zip(legacy, uniform):
        assert left.dtype == right.dtype and left.tobytes() == right.tobytes()
    np.testing.assert_array_equal(legacy[1], weighted[1])
    assert not np.array_equal(legacy[0], weighted[0])


def test_effective_count_changes_pseudocount_strength_not_only_relative_weights():
    msa = encode_alignment(["AAAA", "AAAA", "CCCC"])
    full = compute_pssm(msa, MATRIX, BACKGROUND, sequence_weights=profile_sequence_weights(msa, "henikoff"))
    reduced = compute_pssm(msa, MATRIX, BACKGROUND, sequence_weights=profile_sequence_weights(msa, "henikoff_cluster"))
    assert not np.array_equal(full[0], reduced[0])
    np.testing.assert_array_equal(full[1], reduced[1])


def profile_fixture(tmp_path):
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 3
    (tmp_path / "a.fa").write_text(
        f">a\n{sequence}\n>b\n{sequence}\n>c\n{sequence.replace('A', 'C')}\n"
    )
    names = ["a", "b", "c"]
    database = profile_expansion.load_global_sequence_database(str(tmp_path), ["a.fa"], names)
    return names, database


@pytest.mark.parametrize("mode", PROFILE_SEQUENCE_WEIGHTING)
def test_real_profile_build_serial_spawn_and_self_scoring(tmp_path, mode):
    names, database = profile_fixture(tmp_path)
    # Two tasks ensure the spawn path, with identical deterministic inputs.
    arguments = ([[0, 1, 2], [0, 1, 2]], names, database, "BLOSUM62")
    serial = profile_expansion.build_cluster_profiles(
        *arguments, cpu=1, sequence_weighting=mode, gene_to_species=[0, 0, 1]
    )
    parallel = profile_expansion.build_cluster_profiles(
        *arguments, cpu=2, sequence_weighting=mode, gene_to_species=[0, 0, 1]
    )
    assert set(serial) == set(parallel) == {0, 1}
    for key in serial:
        left, right = serial[key], parallel[key]
        assert left.sequence_weighting == right.sequence_weighting == mode
        assert left.effective_sequence_count == right.effective_sequence_count
        for field in ("match_emissions", "insert_emissions", "transitions", "consensus"):
            np.testing.assert_array_equal(getattr(left, field), getattr(right, field))
    thresholds = profile_expansion.compute_profile_self_thresholds(serial, names, database, cpu=1)
    assert all(np.isfinite(value) and value > 0 for value in thresholds.values())


@pytest.mark.parametrize("mode", PROFILE_SEQUENCE_WEIGHTING[1:])
def test_jackknife_rebuild_preserves_policy_and_species(tmp_path, monkeypatch, mode):
    names, database = profile_fixture(tmp_path)
    profiles = profile_expansion.build_cluster_profiles(
        [[0, 1, 2]], names, database, "BLOSUM62", cpu=1,
        sequence_weighting=mode, gene_to_species=[0, 0, 1],
    )
    original = profile_expansion._build_profile_worker
    observed = []

    def capture(task):
        result = original(task)
        observed.append((task, result[1]))
        return result

    monkeypatch.setattr(profile_expansion, "_build_profile_worker", capture)
    strict = profile_expansion.compute_profile_self_thresholds(profiles, names, database, cpu=1)
    calibrated = profile_expansion.compute_profile_self_thresholds(
        profiles, names, database, cpu=1, matrix_name="BLOSUM62", calibrate_weakest_member=True
    )
    assert len(observed) == 1 and observed[0][1].sequence_weighting == mode
    assert calibrated[0] <= strict[0]
    if mode == "species_henikoff_cluster":
        task, rebuilt = observed[0]
        mapping = dict(zip(names, [0, 0, 1]))
        assert rebuilt.member_species == [mapping[name] for name in task[1]]


def test_cli_default_unchanged_and_weighted_option_passed(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    explicit = process_args(parser.parse_args([str(tmp_path), "--profile-sequence-weighting", "none"]))
    assert default == explicit and "profile_sequence_weighting" not in default
    for mode in PROFILE_SEQUENCE_WEIGHTING[1:]:
        result = process_args(parser.parse_args([
            str(tmp_path), "--accuracy_profile", "high_sensitivity", "--profile-sequence-weighting", mode,
        ]))
        assert result["profile_sequence_weighting"] == mode


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_mode", "phmmer"]])
def test_cli_rejects_inactive_scope(tmp_path, monkeypatch, extra):
    processing = importlib.import_module("orthohmm.args_processing")
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit) as caught:
        process_args(create_parser().parse_args([str(tmp_path), "--profile-sequence-weighting", "henikoff", *extra]))
    assert caught.value.code == 1


@pytest.mark.parametrize("accuracy,search,start", [
    ("standard", "builtin", None), ("high_sensitivity", "phmmer", None),
    ("high_sensitivity", "builtin", StartStep.search_res),
])
def test_api_rejects_inactive_scope_before_search(tmp_path, accuracy, search, start):
    module = importlib.import_module("orthohmm.orthohmm")
    with pytest.raises(ValueError, match="Profile sequence weighting requires"):
        module.execute(str(tmp_path), str(tmp_path / "out"), "phmmer", 1, .5,
            "mcl", 1.5, start, None, SubstitutionMatrix.blosum62, 1e-4,
            accuracy_profile=accuracy, search_mode=search, profile_sequence_weighting="henikoff")


def test_replay_rejects_weighting_without_sequences_before_creating_outputs(tmp_path):
    from benchmark_tools.replay_high_sensitivity import main
    with pytest.raises(SystemExit, match="requires --fasta-directory"):
        main(["--hits-pickle", str(tmp_path / "missing"), "--output-directory", str(tmp_path / "out"),
              "--json", str(tmp_path / "result.json"), "--profile-sequence-weighting", "henikoff"])
    assert not (tmp_path / "out").exists()


def test_benchmark_harness_rejects_inactive_weighting():
    from benchmark_tools.benchmark_production import parse_args
    with pytest.raises(SystemExit):
        parse_args(["in", "out", "metrics.json", "--profile-sequence-weighting", "henikoff"])


def test_expand_profiles_reports_effective_counts_for_real_scoring(tmp_path):
    names, database = profile_fixture(tmp_path)
    result = profile_expansion.expand_profiles(
        [[0, 1, 2]], names, str(tmp_path), ["a.fa"], "BLOSUM62", 1, 1e-4,
        np.array([0, 1, 2]), np.array([1, 2, 0]), np.array([1., 1., 1.]),
        sequence_weighting="species_henikoff_cluster", gene_to_species=[0, 0, 1],
    )
    assert result.profiles_built == 1 and result.profile_candidates > 0
    assert result.profile_effective_counts == ((0, 1.),)
    json.dumps(result.profile_effective_counts)


@pytest.mark.parametrize("mode", PROFILE_SEQUENCE_WEIGHTING[1:])
def test_complete_native_phylogeny_fixture_records_active_weighting(tmp_path, mode):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from benchmark_tools.run_full_search_hillclimb import native_command, validate_native, PRODUCTS
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    source = Path(__file__).resolve().parents[2]
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools):
        pytest.skip("Requires MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    folder = tmp_path / "native"
    folder.mkdir()
    command = native_command(source, fasta, folder, tools, 100, 2)
    command.extend(["--search-result-order", "canonical", "--profile-sequence-weighting", mode])
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "search_max_candidates_per_query": 100,
            "inputs": [record(path) for path in sorted(fasta.iterdir())],
            "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    actual = validate_native(case, 2)
    assert actual["metrics"]["metadata"]["profile_sequence_weighting"] == mode
    counts = actual["metrics"]["counts"]
    assert counts["high_sensitivity_weighted_profiles"] == counts["high_sensitivity_profiles"] > 0
    diagnostics = json.loads((folder / "output/orthohmm_working_res/profile_sequence_weighting.json").read_text())
    assert diagnostics["mode"] == mode and len(diagnostics["profiles"]) == counts["high_sensitivity_profiles"]
    assert all(row["effective_sequence_count"] > 0 for row in diagnostics["profiles"])
