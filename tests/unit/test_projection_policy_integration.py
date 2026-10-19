import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from benchmark_tools.hmm_generalization_holdout import species_sequences
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser
from orthohmm.search import profile_expansion as expansion
from orthohmm.search.matrices import get_matrix
from orthohmm.search.msa_center_star import center_star_msa, validate_profile_projection_policy


BASE = "ACDEFGHIKLMNPQRSTVWY" * 3
POLICY = "residue_preserving"


def family():
    names = ["a", "b", "c"]
    return names, species_sequences(names, [BASE, BASE, "GGG" + BASE], "fixture")


@pytest.mark.parametrize("value", [None, True, [], "unknown"])
def test_invalid_policy_stops_before_build_or_input_loading(value):
    with pytest.raises(ValueError, match="projection policy"):
        validate_profile_projection_policy(value)
    with pytest.raises(ValueError, match="projection policy"):
        expansion.build_cluster_profiles([], [], None, "BLOSUM62", 1, projection_policy=value)
    with pytest.raises(ValueError, match="projection policy"):
        expansion.expand_profiles([], [], "missing", [], "BLOSUM62", 1, 1e-4,
                                  [], [], [], projection_policy=value)


@pytest.mark.parametrize("sequences,expected", [([], []), (["ABX-"], ["AXXX"])])
def test_degenerate_msa_normalizes_only_in_opt_in_policy(sequences, expected):
    assert center_star_msa(sequences, get_matrix("BLOSUM62"), projection_policy=POLICY) == expected
    assert center_star_msa(sequences, get_matrix("BLOSUM62")) == sequences


def test_real_worker_builds_preserve_policy_across_spawn_and_jackknife(monkeypatch):
    names, database = family()
    arguments = ([[0, 1, 2], [0, 1, 2]], names, database, "BLOSUM62")
    serial = expansion.build_cluster_profiles(*arguments, cpu=1, projection_policy=POLICY)
    parallel = expansion.build_cluster_profiles(*arguments, cpu=2, projection_policy=POLICY)
    assert set(serial) == set(parallel) == {0, 1}
    for key in serial:
        assert serial[key].projection_policy == parallel[key].projection_policy == POLICY
        for field in ("match_emissions", "insert_emissions", "transitions", "consensus"):
            np.testing.assert_array_equal(getattr(serial[key], field), getattr(parallel[key], field))
    original = expansion._build_profile_worker
    observed = []

    def capture(task):
        result = original(task)
        observed.append(result[1])
        return result

    monkeypatch.setattr(expansion, "_build_profile_worker", capture)
    strict = expansion.compute_profile_self_thresholds(serial, names, database, cpu=1)
    calibrated = expansion.compute_profile_self_thresholds(
        serial, names, database, cpu=1, matrix_name="BLOSUM62", calibrate_weakest_member=True)
    assert len(observed) == 2 and all(p.projection_policy == POLICY for p in observed)
    assert all(np.isfinite(calibrated[k]) and calibrated[k] <= strict[k] for k in strict)


def test_expansion_passes_policy_to_real_profiles_and_scoring(tmp_path, monkeypatch):
    names, _ = family()
    (tmp_path / "a.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in
                                        zip(names, [BASE, BASE, "GGG" + BASE])))
    observed = []
    original = expansion.build_cluster_profiles

    def capture(*args, **kwargs):
        profiles = original(*args, **kwargs)
        observed.extend(profiles.values())
        return profiles

    monkeypatch.setattr(expansion, "build_cluster_profiles", capture)
    result = expansion.expand_profiles([[0, 1, 2]], names, str(tmp_path), ["a.fa"], "BLOSUM62", 1,
                                      1e-4, [0, 1, 2], [1, 2, 0], [1., 1., 1.], projection_policy=POLICY)
    assert result.profiles_built == 1 and result.profile_candidates > 0
    assert observed[0].projection_policy == POLICY


def test_cli_default_call_shape_is_unchanged_and_opt_in_is_passed(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    explicit = process_args(parser.parse_args([str(tmp_path), "--profile-projection-policy", "legacy"]))
    assert default == explicit and "profile_projection_policy" not in default
    active = process_args(parser.parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity",
                                            "--profile-projection-policy", POLICY]))
    assert active["profile_projection_policy"] == POLICY


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_mode", "phmmer"]])
def test_cli_rejects_inactive_scopes(tmp_path, monkeypatch, extra):
    processing = importlib.import_module("orthohmm.args_processing")
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit):
        process_args(create_parser().parse_args([str(tmp_path), "--profile-projection-policy", POLICY, *extra]))


@pytest.mark.parametrize("accuracy,search,start", [("standard", "builtin", None),
    ("high_sensitivity", "phmmer", None), ("high_sensitivity", "builtin", StartStep.search_res)])
def test_api_rejects_inactive_scopes_before_search(tmp_path, accuracy, search, start):
    module = importlib.import_module("orthohmm.orthohmm")
    with pytest.raises(ValueError, match="Profile projection repair requires"):
        module.execute(str(tmp_path), str(tmp_path / "out"), "phmmer", 1, .5, "mcl", 1.5,
                       start, None, SubstitutionMatrix.blosum62, 1e-4, accuracy_profile=accuracy,
                       search_mode=search, profile_projection_policy=POLICY)


@pytest.mark.parametrize("option,value,api_key", [("--profile-sequence-weighting", "henikoff", "sequence_weighting"),
    ("--profile-alignment-policy", "strict_span_guard", "alignment_policy")])
def test_unregistered_combinations_rejected_at_all_entry_points(tmp_path, capsys, option, value, api_key):
    from benchmark_tools.benchmark_production import parse_args
    from benchmark_tools.replay_high_sensitivity import main
    with pytest.raises(SystemExit):
        process_args(create_parser().parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity",
            "--profile-projection-policy", POLICY, option, value]))
    with pytest.raises(SystemExit):
        parse_args(["in", "out", "metrics.json", "--cpu", "1", "--accuracy-profile", "high_sensitivity",
                    "--profile-projection-policy", POLICY, option, value])
    assert "separate experiment" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="separate experiment"):
        main(["--hits-pickle", "missing", "--output-directory", str(tmp_path / "out"),
              "--json", str(tmp_path / "out.json"), "--fasta-directory", str(tmp_path),
              "--profile-projection-policy", POLICY, option, value])
    with pytest.raises(ValueError, match="separate experiment"):
        expansion.expand_profiles([], [], "missing", [], "BLOSUM62", 1, 1e-4, [], [], [],
                                  projection_policy=POLICY, **{api_key: value})
    assert not (tmp_path / "out").exists()
    module = importlib.import_module("orthohmm.orthohmm")
    with pytest.raises(ValueError, match="separate experiment"):
        module.execute(str(tmp_path), str(tmp_path / "out"), "phmmer", 1, .5, "mcl", 1.5,
                       None, None, SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity",
                       profile_projection_policy=POLICY, **{"profile_" + api_key: value})
    # The production API creates output/working directories before validation,
    # but must not run inference or write biological products for this rejection.
    assert not [p for p in (tmp_path / "out").rglob("*") if p.is_file()]


def test_replay_and_harness_reject_inactive_policy_before_outputs(tmp_path, capsys):
    from benchmark_tools.benchmark_production import parse_args
    from benchmark_tools.replay_high_sensitivity import main
    with pytest.raises(SystemExit):
        parse_args(["in", "out", "metrics.json", "--cpu", "1", "--profile-projection-policy", POLICY])
    assert "requires --accuracy-profile high_sensitivity" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="requires --fasta-directory"):
        main(["--hits-pickle", "missing", "--output-directory", str(tmp_path / "out"),
              "--json", str(tmp_path / "out.json"), "--profile-projection-policy", POLICY])
    assert not (tmp_path / "out").exists()


def test_numeric_hit_replay_builds_repaired_profiles_and_records_policy(tmp_path, monkeypatch):
    from orthohmm.accuracy import write_accuracy_checkpoint
    from benchmark_tools.orthobench_stage_diagnostics import file_provenance
    from benchmark_tools.replay_high_sensitivity import main
    names = [f"gene{i}" for i in range(6)]
    fasta = tmp_path / "input"
    fasta.mkdir()
    sequences = [BASE, BASE, "GGG" + BASE, BASE[::-1], BASE[::-1], "GGG" + BASE[::-1]]
    for i in range(3):
        (fasta / f"taxon{i}.fa").write_text("".join(f">{names[j]}\n{sequences[j]}\n" for j in (i, i + 3)))
    pairs = [(q, t) for q in range(6) for t in range(6) if q // 3 == t // 3]
    checkpoint = write_accuracy_checkpoint(str(tmp_path / "checkpoint"), names, [0, 1, 2, 0, 1, 2],
                                          [p[0] for p in pairs], [p[1] for p in pairs], [100.] * len(pairs))
    checksum = file_provenance(checkpoint / "manifest.json")["sha256"]
    observed = []
    original = expansion.build_cluster_profiles

    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        observed.extend(result.values())
        return result

    monkeypatch.setattr(expansion, "build_cluster_profiles", capture)
    output = tmp_path / "graph"
    result_json = tmp_path / "result.json"
    assert main(["--accuracy-checkpoint", str(checkpoint), "--checkpoint-sha256", checksum,
                 "--output-directory", str(output), "--json", str(result_json), "--fasta-directory", str(fasta),
                 "--profile-projection-policy", POLICY, "--cpu", "1"]) == 0
    assert observed and all(p.projection_policy == POLICY for p in observed)
    assert json.loads(result_json.read_text())["parameters"]["profile_projection_policy"] == POLICY


def test_complete_native_pipeline_propagates_opt_in_through_harness_and_cli(tmp_path):
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
    command.extend(["--search-result-order", "canonical", "--profile-projection-policy", POLICY])
    environment = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(command, cwd=source, env=environment, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "search_max_candidates_per_query": 100,
            "inputs": [record(path) for path in sorted(fasta.iterdir())],
            "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    actual = validate_native(case, 2)
    assert actual["metrics"]["metadata"]["profile_projection_policy"] == POLICY
    counts = actual["metrics"]["counts"]
    assert counts["high_sensitivity_residue_preserving_profiles"] == counts["high_sensitivity_profiles"] > 0
