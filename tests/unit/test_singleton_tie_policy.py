import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from orthohmm.accuracy import build_singleton_assignment_edges, validate_singleton_tie_policy
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser


def data():
    names = ["s", "a", "b", "c", "d", "unique"]
    clusters = [[0], [1, 2], [3, 4], [5]]
    q = np.array([0, 0, 0, 0, 5, 5], dtype=np.int32)
    t = np.array([1, 3, 2, 4, 3, 4], dtype=np.int32)
    s = np.array([.8, .8, .7, .6, .9, .5], dtype=np.float64)
    return names, clusters, q, t, s


def records(edges):
    return list(zip(edges.sources.tolist(), edges.targets.tolist(), edges.weights.tolist()))


def test_default_and_explicit_first_hit_are_byte_identical():
    args = data()
    a = build_singleton_assignment_edges(*args)
    b = build_singleton_assignment_edges(*args, tie_policy="first_hit")
    assert a.gene_names == b.gene_names
    for field in ("sources", "targets", "weights"):
        assert getattr(a, field).tobytes() == getattr(b, field).tobytes()


def test_defer_ambiguous_removes_only_ambiguous_query_edges_and_preserves_genes():
    args = data()
    result = build_singleton_assignment_edges(*args, tie_policy="defer_ambiguous")
    assert records(result) == [(3, 5, .9), (4, 5, .5)]
    assert result.gene_names == args[0]
    assert records(result) == [r for r in records(build_singleton_assignment_edges(*args)) if 0 not in r[:2]]


def test_fixed_partition_assignment_is_permutation_invariant():
    names, clusters, q, t, s = data()
    expected = records(build_singleton_assignment_edges(names, clusters, q, t, s, tie_policy="defer_ambiguous"))
    rng = np.random.default_rng(123)
    for _ in range(20):
        order = rng.permutation(len(q))
        assert records(build_singleton_assignment_edges(names, clusters, q[order], t[order], s[order],
            tie_policy="defer_ambiguous")) == expected


def test_tied_genes_within_one_cluster_and_cluster_zero_remain_assigned():
    result = build_singleton_assignment_edges(["a", "b", "s"], [[0, 1], [2]],
        np.array([2, 2, 2]), np.array([0, 1, 1]), np.array([.8, .8, .7]), tie_policy="defer_ambiguous")
    assert records(result) == [(0, 2, .8), (1, 2, .8)]


def test_near_ties_are_not_treated_as_exact_ties():
    names, clusters, q, t, s = data()
    s[1] = np.nextafter(s[0], 0.)
    assert records(build_singleton_assignment_edges(names, clusters, q, t, s, tie_policy="defer_ambiguous")) == records(
        build_singleton_assignment_edges(names, clusters, q, t, s))


def test_empty_inputs_return_typed_empty_graph():
    result = build_singleton_assignment_edges(["s"], [[0]], [], [], [], tie_policy="defer_ambiguous")
    assert result.gene_names == ["s"] and len(result.sources) == 0
    assert result.sources.dtype == result.targets.dtype == np.int32 and result.weights.dtype == np.float64


@pytest.mark.parametrize("policy", [None, True, [], {}, "unknown", "legacy"])
def test_invalid_policies_are_rejected_even_for_empty_input(policy):
    with pytest.raises(ValueError, match="singleton tie policy"):
        validate_singleton_tie_policy(policy)
    with pytest.raises(ValueError, match="singleton tie policy"):
        build_singleton_assignment_edges([], [], [], [], [], tie_policy=policy)


def test_cli_default_call_shape_is_unchanged_and_variant_is_explicit(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    explicit = process_args(parser.parse_args([str(tmp_path), "--singleton_tie_policy", "first_hit"]))
    assert default == explicit and "singleton_tie_policy" not in default
    variant = process_args(parser.parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity",
        "--singleton_tie_policy", "defer_ambiguous"]))
    assert variant["singleton_tie_policy"] == "defer_ambiguous"


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_mode", "phmmer"]])
def test_cli_rejects_inactive_or_resumed_scope(tmp_path, monkeypatch, extra):
    import orthohmm.args_processing as processing
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit) as error:
        process_args(create_parser().parse_args([str(tmp_path), "--singleton_tie_policy", "defer_ambiguous", *extra]))
    assert error.value.code == 1


@pytest.mark.parametrize("profile,search,start", [("standard", "builtin", None),
    ("high_sensitivity", "phmmer", None), ("high_sensitivity", "builtin", StartStep.search_res)])
def test_api_rejects_inactive_scope_before_search(tmp_path, profile, search, start):
    from orthohmm.orthohmm import execute
    with pytest.raises(ValueError, match="fresh built-in high_sensitivity"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, .5, "mcl", 1.5, start, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile=profile, search_mode=search,
            singleton_tie_policy="defer_ambiguous")


@pytest.mark.parametrize("policy", ["first_hit", "defer_ambiguous"])
def test_policy_metadata_does_not_change_hmm_search_configuration(tmp_path, monkeypatch, policy):
    from orthohmm.orthohmm import execute
    import orthohmm.search.engine as engine
    class SearchReached(Exception):
        pass
    captured = {}
    def search(*args, **kwargs):
        captured.update(kwargs)
        raise SearchReached()
    monkeypatch.setattr(engine, "execute_builtin_search", search)
    (tmp_path / "s.faa").write_text(">gene\nACDEFGHIK\n")
    metrics = tmp_path / "metrics.json"
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, .5, "mcl", 1.5, None, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity",
            singleton_tie_policy=policy, metrics_json=str(metrics))
    assert "singleton_tie_policy" not in captured and "tie_policy" not in captured
    assert captured["kmer_k"] == 4 and captured["max_candidates_per_query"] == 100
    assert json.loads(metrics.read_text())["metadata"]["singleton_tie_policy"] == policy


def test_cached_replay_parser_and_production_wrapper_share_policy_names():
    from benchmark_tools.replay_high_sensitivity import build_parser
    from benchmark_tools.benchmark_production import parse_args
    replay = build_parser().parse_args(["--accuracy-checkpoint", "input", "--checkpoint-sha256", "hash",
        "--output-directory", "out", "--json", "result.json", "--singleton-tie-policy", "defer_ambiguous"])
    assert replay.singleton_tie_policy == "defer_ambiguous"
    with pytest.raises(SystemExit):
        parse_args(["inputs", "output", "metrics.json", "--singleton-tie-policy", "defer_ambiguous"])


@pytest.mark.parametrize("policy", ["first_hit", "defer_ambiguous"])
def test_pipeline_passes_policy_to_both_graph_attachment_stages(tmp_path, monkeypatch, policy):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    import orthohmm.orthohmm as pipeline
    import orthohmm.search.engine as engine
    import orthohmm.search.profile_expansion as expansion
    from orthohmm.helpers import IndexedEdges
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    fasta = tmp_path / "input"
    fasta.mkdir()
    species_ids = {}
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
        species_ids[f"{taxon}.faa"] = list(genes)
    monkeypatch.setattr(engine, "execute_builtin_search", lambda *a, **k: engine.IndexedSearchResults({}, species_ids))
    def profiles(clusters, names, *args, **kwargs):
        edges = IndexedEdges(names, np.array([], dtype=np.int32), np.array([], dtype=np.int32), np.array([], dtype=np.float64))
        return expansion.ProfileExpansionResult(edges, 0, 0, 0)
    monkeypatch.setattr(expansion, "expand_profiles", profiles)
    captured = []
    original = pipeline.build_singleton_assignment_edges
    class BothStagesReached(Exception):
        pass
    def attach(*args, **kwargs):
        captured.append(dict(kwargs))
        result = original(*args, **kwargs)
        if len(captured) == 2:
            raise BothStagesReached()
        return result
    monkeypatch.setattr(pipeline, "build_singleton_assignment_edges", attach)
    with pytest.raises(BothStagesReached):
        pipeline.execute(str(fasta), str(tmp_path / "output"), None, 2, .5, "mcl", 1.5, None, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity", singleton_tie_policy=policy)
    expected = {} if policy == "first_hit" else {"tie_policy": policy}
    assert captured == [expected, expected]


def test_real_native_phylogeny_pipeline_records_opt_in_policy(tmp_path):
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
    command.extend(["--search-result-order", "canonical", "--singleton-tie-policy", "defer_ambiguous"])
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr + (folder / "output/benchmark.log").read_text()
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "search_max_candidates_per_query": 100,
            "inputs": [record(p) for p in sorted(fasta.iterdir())],
            "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    actual = validate_native(case, 2)
    assert actual["metrics"]["metadata"]["singleton_tie_policy"] == "defer_ambiguous"
    assert (actual["genes"], actual["species"]) == (16, 4)
