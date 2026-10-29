import importlib
import hashlib
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
from orthohmm.search.msa_profile import MSAProfile
from orthohmm.search.profile_spans import validate_profile_alignment_policy


def model(length, members, match=5):
    emissions = np.full((length, 20), -20, dtype=np.int8)
    emissions[:, 0] = match
    return MSAProfile(length, emissions, np.full(20, -1, dtype=np.int8),
                      np.array([0, -12, -12, -1, -3, -1, -3], dtype=np.int32),
                      np.zeros(length, dtype=np.uint8), list(members))


def fixture():
    names = ["a0", "a1", "a2", "full", "fragment", "partial"]
    database = species_sequences(names, ["A" * 20] * 4 + ["A" * 10, "C" * 10 + "A" * 5 + "C" * 10], "fixture")
    profiles = {0: model(20, names[:3])}
    hits = expansion.ProfileHits(np.array([0, 0, 0], dtype=np.int32), np.array([3, 4, 5], dtype=np.int32),
                                  np.array([100., 50., 25.]), np.zeros(3), 3)
    arguments = (names, [[0, 1, 2], [3], [4], [5]], hits, {0: 25.},
                 np.array([3, 4, 5]), np.array([0, 0, 0]), np.array([.9, .8, .7]))
    return arguments, profiles, database


def records(edges):
    return list(zip(edges.sources.tolist(), edges.targets.tolist(), edges.weights.tolist()))


def randomized_default_edge_digest(builder):
    """Seeded selection/tie/duplicate fixtures bound to pre-guard76fea958."""
    rng = np.random.default_rng(18102026)
    digest = hashlib.sha256()
    names = [f"gene_{index}" for index in range(12)]
    clusters = [list(range(index, index + 3)) for index in range(0, 12, 3)]
    for _ in range(256):
        hits = expansion.ProfileHits(rng.integers(0, 4, 64, dtype=np.int32),
            rng.integers(0, 12, 64, dtype=np.int32), rng.integers(1, 51, 64).astype(np.float64),
            np.zeros(64), 64)
        thresholds = {index: float(rng.integers(1, 51)) for index in range(4)}
        edges = builder(names, clusters, hits, thresholds,
            rng.integers(0, 12, 96, dtype=np.int32), rng.integers(0, 12, 96, dtype=np.int32),
            rng.integers(1, 11, 96).astype(np.float64))
        digest.update(len(edges.sources).to_bytes(8, "little"))
        for field, dtype in (("sources", "<i4"), ("targets", "<i4"), ("weights", "<f8")):
            digest.update(getattr(edges, field).astype(dtype).tobytes())
    return digest.hexdigest()


def test_256_randomized_default_edge_fixtures_match_committed_pre_guard_digest():
    assert randomized_default_edge_digest(expansion.build_strict_profile_edges) == (
        "dd778e42fb41d137899000ffbf36f67ca5495a228581120b639c89ad87646fe0"
    )


def test_default_and_explicit_off_are_byte_identical():
    arguments, _, _ = fixture()
    original = expansion.build_strict_profile_edges(*arguments)
    explicit = expansion.build_strict_profile_edges(*arguments, alignment_policy="off")
    assert records(original) == [(0, 3, .9), (0, 4, .8), (0, 5, .7)]
    for field in ("sources", "targets", "weights"):
        assert getattr(original, field).tobytes() == getattr(explicit, field).tobytes()


@pytest.mark.parametrize("cpu", [1, 4])
def test_guard_preserves_full_and_fragment_links_but_removes_partial_shared_region(cpu):
    arguments, profiles, database = fixture()
    trace = []
    guarded = expansion.build_strict_profile_edges(*arguments, alignment_policy="strict_span_guard",
        profiles=profiles, sequence_database=database, alignment_trace=trace, cpu=cpu)
    assert records(guarded) == [(0, 3, .9), (0, 4, .8)]
    assert [row["gene_id"] for row in trace] == [3, 4, 5]
    assert [row["category"] for row in trace] == ["both_sides_covered", "asymmetric_coverage", "partial_shared_region"]
    assert [row["score"] for row in trace] == [100, 50, 25]
    assert [row["passes_guard"] for row in trace] == [True, True, False]


def test_serial_and_parallel_evidence_and_edge_bytes_are_identical():
    arguments, profiles, database = fixture()
    traces, edges = [], []
    for cpu in (1, 4):
        trace = []
        edges.append(expansion.build_strict_profile_edges(*arguments, alignment_policy="strict_span_guard",
            profiles=profiles, sequence_database=database, alignment_trace=trace, cpu=cpu))
        traces.append(trace)
    assert traces[0] == traces[1]
    for field in ("sources", "targets", "weights"):
        assert getattr(edges[0], field).tobytes() == getattr(edges[1], field).tobytes()


def test_guard_does_not_switch_to_a_lower_scoring_profile_after_rejection():
    names = ["a0", "a1", "a2", "b0", "b1", "b2", "query"]
    profiles = {0: model(20, names[:3], match=6), 1: model(5, names[3:6])}
    database = species_sequences(names, ["A" * 20] * 3 + ["A" * 5] * 3 + ["C" * 10 + "A" * 5 + "C" * 10], "fixture")
    hits = expansion.ProfileHits(np.array([0, 1]), np.array([6, 6]), np.array([30., 25.]), np.zeros(2), 2)
    trace = []
    edges = expansion.build_strict_profile_edges(names, [[0, 1, 2], [3, 4, 5], [6]], hits,
        {0: 20., 1: 20.}, np.array([6, 6]), np.array([0, 3]), np.array([1., .9]),
        alignment_policy="strict_span_guard", profiles=profiles, sequence_database=database, alignment_trace=trace)
    assert not records(edges)
    assert len(trace) == 1 and trace[0]["profile_cluster_id"] == 0 and not trace[0]["passes_guard"]


def test_score_mismatch_stops_the_intervention_without_partial_trace():
    arguments, profiles, database = fixture()
    arguments[2].scores[0] += 1
    trace = []
    with pytest.raises(ValueError, match="score mismatch"):
        expansion.build_strict_profile_edges(*arguments, alignment_policy="strict_span_guard",
            profiles=profiles, sequence_database=database, alignment_trace=trace)
    assert trace == []


@pytest.mark.parametrize("policy", [None, True, [], "unknown"])
def test_invalid_policy_rejected_before_any_inference(policy):
    with pytest.raises(ValueError, match="profile alignment policy"):
        validate_profile_alignment_policy(policy)
    with pytest.raises(ValueError, match="profile alignment policy"):
        expansion.build_strict_profile_edges([], [], None, {}, [], [], [], alignment_policy=policy)


def test_active_guard_requires_models_sequences_and_explicit_empty_trace():
    arguments, _, _ = fixture()
    with pytest.raises(ValueError, match="requires profiles"):
        expansion.build_strict_profile_edges(*arguments, alignment_policy="strict_span_guard")


@pytest.mark.parametrize("cpu", [0, -1, True, 1.5])
def test_trace_requires_a_real_positive_cpu_budget(cpu):
    arguments, profiles, database = fixture()
    with pytest.raises(ValueError, match="CPU budget"):
        expansion.build_strict_profile_edges(*arguments, alignment_policy="strict_span_guard",
            profiles=profiles, sequence_database=database, alignment_trace=[], cpu=cpu)


def test_real_profile_expansion_traces_and_rejects_nonempty_candidate_set(tmp_path):
    names = ["a0", "a1", "a2", "full", "fragment", "partial"]
    sequences = ["A" * 100, "A" * 100, "A" * 20, "A" * 100, "A" * 20, "C" * 50 + "A" * 20 + "C" * 50]
    (tmp_path / "input.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in zip(names, sequences)))
    arguments = ([[0, 1, 2], [3], [4], [5]], names, str(tmp_path), ["input.fa"], "BLOSUM62", 2,
                 1e-4, np.array([3, 4, 5]), np.array([0, 0, 0]), np.array([.9, .8, .7]))
    control = expansion.expand_profiles(*arguments)
    guarded = expansion.expand_profiles(*arguments, alignment_policy="strict_span_guard")
    assert records(control.edges) == [(0, 3, .9), (0, 4, .8), (0, 5, .7)]
    assert records(guarded.edges) == [(0, 3, .9), (0, 4, .8)]
    assert len(guarded.alignment_evidence) == 3
    assert sum(not row["passes_guard"] for row in guarded.alignment_evidence) == 1
    assert guarded.alignment_profiles == ((0, tuple(names[:3])),)


def test_cli_default_call_shape_unchanged_and_active_policy_is_explicit(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    explicit = process_args(parser.parse_args([str(tmp_path), "--profile-alignment-policy", "off"]))
    assert default == explicit and "profile_alignment_policy" not in default
    guarded = process_args(parser.parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity",
        "--profile-alignment-policy", "strict_span_guard"]))
    assert guarded["profile_alignment_policy"] == "strict_span_guard"


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_mode", "phmmer"],
    ["--accuracy_profile", "high_sensitivity", "--profile-sequence-weighting", "henikoff"]])
def test_cli_rejects_inactive_scope_and_unvalidated_combination(tmp_path, monkeypatch, extra):
    processing = importlib.import_module("orthohmm.args_processing")
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit) as error:
        process_args(create_parser().parse_args([str(tmp_path), "--profile-alignment-policy", "strict_span_guard", *extra]))
    assert error.value.code == 1


@pytest.mark.parametrize("accuracy,search,start", [("standard", "builtin", None),
    ("high_sensitivity", "phmmer", None), ("high_sensitivity", "builtin", StartStep.search_res)])
def test_api_rejects_inactive_scope_before_search(tmp_path, accuracy, search, start):
    module = importlib.import_module("orthohmm.orthohmm")
    with pytest.raises(ValueError, match="Profile alignment guard requires"):
        module.execute(str(tmp_path), str(tmp_path / "out"), "phmmer", 1, .5, "mcl", 1.5,
            start, None, SubstitutionMatrix.blosum62, 1e-4, accuracy_profile=accuracy,
            search_mode=search, profile_alignment_policy="strict_span_guard")




def test_unvalidated_weighting_guard_combination_fails_before_loading_sequences():
    with pytest.raises(ValueError, match="separate experiment"):
        expansion.expand_profiles([], [], "/missing", [], "BLOSUM62", 1, 1e-4, [], [], [],
            sequence_weighting="henikoff", alignment_policy="strict_span_guard")


@pytest.mark.parametrize("policy", ["off", "strict_span_guard"])
def test_production_passes_policy_only_when_requested(tmp_path, monkeypatch, policy):
    from tests.support import fixture as native_fixture
    engine = importlib.import_module("orthohmm.search.engine")
    pipeline = importlib.import_module("orthohmm.orthohmm")
    fasta = tmp_path / "input"
    fasta.mkdir()
    species_ids = {}
    for taxon, genes in native_fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
        species_ids[f"{taxon}.faa"] = list(genes)
    monkeypatch.setattr(engine, "execute_builtin_search", lambda *args, **kwargs: engine.IndexedSearchResults({}, species_ids))
    captured = []

    class ExpansionReached(Exception):
        pass

    def stop(*args, **kwargs):
        captured.append(kwargs)
        raise ExpansionReached()

    monkeypatch.setattr(expansion, "expand_profiles", stop)
    with pytest.raises(ExpansionReached):
        pipeline.execute(str(fasta), str(tmp_path / "output"), None, 2, .5, "mcl", 1.5, None, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity", profile_alignment_policy=policy)
    assert captured == ([{}] if policy == "off" else [{"alignment_policy": policy}])
