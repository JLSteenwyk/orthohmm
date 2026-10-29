import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from tests.support import species_sequences
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


@pytest.mark.parametrize("weighting", ["none", "henikoff", "henikoff_cluster"])
def test_real_worker_builds_preserve_policy_across_spawn_and_jackknife(monkeypatch, weighting):
    names, database = family()
    arguments = ([[0, 1, 2], [0, 1, 2]], names, database, "BLOSUM62")
    serial = expansion.build_cluster_profiles(*arguments, cpu=1, projection_policy=POLICY, sequence_weighting=weighting)
    parallel = expansion.build_cluster_profiles(*arguments, cpu=2, projection_policy=POLICY, sequence_weighting=weighting)
    assert set(serial) == set(parallel) == {0, 1}
    for key in serial:
        assert serial[key].projection_policy == parallel[key].projection_policy == POLICY
        assert serial[key].sequence_weighting == parallel[key].sequence_weighting == weighting
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
    assert all(p.sequence_weighting == weighting for p in observed)
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


@pytest.mark.parametrize("option,value,api_key", [
    ("--profile-sequence-weighting", "species_henikoff_cluster", "sequence_weighting"),
    ("--profile-alignment-policy", "strict_span_guard", "alignment_policy"),
])
def test_unregistered_combinations_rejected_at_application_entry_points(tmp_path, option, value, api_key):
    with pytest.raises(SystemExit):
        process_args(create_parser().parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity",
            "--profile-projection-policy", POLICY, option, value]))
    with pytest.raises(ValueError, match="separate experiment"):
        expansion.expand_profiles([], [], "missing", [], "BLOSUM62", 1, 1e-4, [], [], [],
                                  projection_policy=POLICY, **{api_key: value})
    assert not (tmp_path / "out").exists()
    module = importlib.import_module("orthohmm.orthohmm")
    with pytest.raises(ValueError, match="separate experiment"):
        module.execute(str(tmp_path), str(tmp_path / "out"), "phmmer", 1, .5, "mcl", 1.5,
                       None, None, SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity",
                       profile_projection_policy=POLICY, **{"profile_" + api_key: value})
    assert not [path for path in (tmp_path / "out").rglob("*") if path.is_file()]
