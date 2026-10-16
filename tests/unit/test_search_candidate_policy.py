import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.accuracy import (
    load_accuracy_checkpoint, resolve_search_candidate_policy,
    validate_search_candidate_policy, write_accuracy_checkpoint,
)
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser
from orthohmm.search import engine
from orthohmm.search.prefilter import union_candidate_sets


@pytest.mark.parametrize("value", [None, True, 1, "", "union", [], {}])
def test_unknown_policies_refused(value):
    with pytest.raises(ValueError, match="Unsupported search candidate policy"):
        validate_search_candidate_policy(value)


@pytest.mark.parametrize("cap", [0, -1, None, True, np.bool_(True), 1., "100", 2**31])
def test_union_requires_fixed_positive_integer_cap(cap):
    with pytest.raises(ValueError):
        resolve_search_candidate_policy("protected_union", max_candidates_per_query=cap)


@pytest.mark.parametrize("minimum", [1, 3, 5, True, 4., "4"])
def test_union_baseline_is_exactly_min4(minimum):
    with pytest.raises(ValueError):
        resolve_search_candidate_policy("protected_union", min_total_hits=minimum, max_candidates_per_query=100)


def test_single_retains_adaptive_settings_and_union_allows_generic_fixed_caps():
    assert resolve_search_candidate_policy() == "single"
    for cap in (1, 100, np.int64(400)):
        assert resolve_search_candidate_policy("protected_union", max_candidates_per_query=cap) == "protected_union"


def test_csr_union_is_sorted_deduplicated_and_does_not_retruncate():
    ids, offsets = union_candidate_sets(np.array([2, 0, 1]), np.array([0, 2, 3]),
        np.array([1, 2, 0]), np.array([0, 2, 3]), 3, 2)
    assert ids.tolist() == [0, 1, 2, 0, 1] and offsets.tolist() == [0, 3, 5]
    assert ids.dtype == np.int32 and offsets.dtype == np.int64


def test_csr_union_matches_python_sets_on_random_bounded_inputs():
    rng = np.random.default_rng(20261016)
    for targets in range(21):
        groups = [[rng.choice(targets, size=rng.integers(min(targets, 5) + 1), replace=False).tolist()
                   for _ in range(7)] for _ in range(2)]
        raw = [(np.array([i for g in arm for i in g], dtype=np.int32),
                np.array([0, *np.cumsum([len(g) for g in arm])], dtype=np.int64)) for arm in groups]
        ids, offsets = union_candidate_sets(*raw[0], *raw[1], targets, 5)
        for q, (left, right) in enumerate(zip(*groups)):
            assert ids[offsets[q]:offsets[q + 1]].tolist() == sorted(set(left) | set(right))


@pytest.mark.parametrize("left,offsets,right,right_offsets,targets,cap", [
    ([0], [0, 1], [], [0], 2, 1), ([0], [1, 1], [0], [0, 1], 2, 1),
    ([0], [0, 2], [0], [0, 1], 2, 1), ([0], [0, 1, 0], [0], [0, 1, 1], 2, 1),
    ([-1], [0, 1], [0], [0, 1], 2, 1), ([2], [0, 1], [0], [0, 1], 2, 1),
    ([0, 0], [0, 2], [0], [0, 1], 2, 2), ([0, 1], [0, 2], [0], [0, 1], 2, 1),
    ([0], [0, 1], [0], [0, 1], 2, 0), ([0], [0, 1], [0], [0, 1], True, 1),
])
def test_malformed_or_unbounded_csr_refused(left, offsets, right, right_offsets, targets, cap):
    with pytest.raises(ValueError):
        union_candidate_sets(np.array(left, dtype=np.int32), np.array(offsets),
            np.array(right, dtype=np.int32), np.array(right_offsets), targets, cap)


def test_empty_csr_query_and_target_universes():
    empty = np.array([], dtype=np.int32)
    ids, offsets = union_candidate_sets(empty, [0], empty, [0], 0, 1)
    assert not len(ids) and offsets.tolist() == [0]
    ids, offsets = union_candidate_sets(empty, [0, 0, 0], empty, [0, 0, 0], 0, 1)
    assert not len(ids) and offsets.tolist() == [0, 0, 0]


def test_unsigned_scalar_dimensions_are_normalized_before_packing():
    ids, offsets = union_candidate_sets(np.array([2], dtype=np.int32), [0, 1],
        np.array([1], dtype=np.int32), [0, 1], np.uint64(3), np.uint64(1))
    assert ids.tolist() == [1, 2] and offsets.tolist() == [0, 2]


def test_engine_rejects_union_with_adaptive_cap_even_for_empty_inputs():
    empty = SimpleNamespace(num_sequences=0)
    with pytest.raises(ValueError, match="fixed positive candidate cap"):
        engine.search_species_pair_indexed(empty, empty, "BLOSUM62", candidate_policy="protected_union")


@pytest.mark.parametrize("policy", ["single", "protected_union"])
def test_engine_reuses_one_index_builds_profiles_once_and_scores_unique_candidates_once(monkeypatch, policy):
    query = SimpleNamespace(species_file="q", num_sequences=2, lengths=np.array([3, 4]),
                            flat_sequences=np.array([], dtype=np.int32), offsets=np.array([0, 3]))
    target = SimpleNamespace(species_file="t", num_sequences=3, lengths=np.array([3, 4, 5]),
                             flat_sequences=np.array([], dtype=np.int32), offsets=np.array([0, 3, 7]))
    calls, indexed, profiled, scored, databases = [], [], [], [], []
    index = object()
    def prepare(*args, **kwargs):
        indexed.append(True)
        return index
    def prefilter(*args, **kwargs):
        assert kwargs["prepared_index"] is index and kwargs["max_candidates_per_query"] == 2
        calls.append(kwargs)
        return (np.array([2, 0, 1], dtype=np.int32) if kwargs["min_total_hits"] == 4
                else np.array([1, 2, 0], dtype=np.int32)), np.array([0, 2, 3], dtype=np.int64)
    def profiles(*args):
        profiled.append(True)
        return (None,) * 5
    def viterbi(*args, **kwargs):
        pairs = args[8]
        scored.append(pairs.copy())
        return 10 * pairs[:, 0] + pairs[:, 1] + 1
    def evalues(scores, lengths, db, *args):
        databases.append(db)
        return np.full(len(scores), 1e-5)
    monkeypatch.setattr(engine, "prepare_kmer_index", prepare)
    monkeypatch.setattr(engine, "prefilter_candidates", prefilter)
    monkeypatch.setattr(engine, "build_profiles_batch", profiles)
    monkeypatch.setattr(engine, "batch_viterbi_multipair_c", viterbi)
    monkeypatch.setattr(engine, "batch_estimate_evalues", evalues)
    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    monkeypatch.setattr(engine, "is_c_available", lambda: True)
    result = engine.search_species_pair_indexed(query, target, "BLOSUM62", max_candidates_per_query=2,
        candidate_policy=policy)
    assert indexed == profiled == [True] and len(scored) == 1 and databases == [12]
    assert [c["min_total_hits"] for c in calls] == ([4] if policy == "single" else [4, 1])
    expected = [(0, 2), (0, 0), (1, 1)] if policy == "single" else [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    assert list(zip(result.query_indices, result.target_indices)) == expected
    assert result.candidate_count == len(expected)
    for q, t, score in zip(result.query_indices, result.target_indices, result.scores):
        assert score == (10 * q + t + 1) / np.sqrt(query.lengths[q] * target.lengths[t])


def test_cli_default_dictionary_is_unchanged_and_both_policy_spellings_work(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    assert default == process_args(parser.parse_args([str(tmp_path), "--search-candidate-policy", "single"]))
    assert "search_candidate_policy" not in default
    for option in ("--search-candidate-policy", "--search_candidate_policy"):
        args = [str(tmp_path), "--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical",
                option, "protected_union"]
        assert process_args(parser.parse_args(args))["search_candidate_policy"] == "protected_union"


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity"],
    ["--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical", "--search_max_candidates_per_query", "0"],
    ["--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical", "--search_min_total_hits", "1"],
    ["--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical", "--search_mode", "phmmer"]])
def test_cli_refuses_unsupported_combinations(tmp_path, monkeypatch, extra):
    import orthohmm.args_processing as processing
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit):
        process_args(create_parser().parse_args([str(tmp_path), "--search-candidate-policy", "protected_union", *extra]))


@pytest.mark.parametrize("change", [{"accuracy_profile": "standard"}, {"search_mode": "phmmer"},
    {"start": StartStep.search_res}, {"search_result_order": "arrival"},
    {"search_max_candidates_per_query": 0}, {"search_min_total_hits": 1}])
def test_api_refuses_unsupported_combinations(tmp_path, change):
    from orthohmm.orthohmm import execute
    options = {"accuracy_profile": "high_sensitivity", "search_mode": "builtin", "start": None,
        "search_result_order": "canonical", "search_candidate_policy": "protected_union", **change}
    start = options.pop("start")
    with pytest.raises(ValueError):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, .5, "mcl", 1.5, start, None,
            SubstitutionMatrix.blosum62, 1e-4, **options)


@pytest.mark.parametrize("policy", ["single", "protected_union"])
def test_pipeline_forwards_only_opt_in_policy_and_records_actual_value(tmp_path, monkeypatch, policy):
    from orthohmm.orthohmm import execute
    class SearchReached(Exception): pass
    captured = {}
    def search(*args, **kwargs):
        captured.update(kwargs)
        raise SearchReached()
    monkeypatch.setattr(engine, "execute_builtin_search", search)
    (tmp_path / "s.faa").write_text(">g\nACDEFGHIK\n")
    metrics = tmp_path / "metrics.json"
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, .5, "mcl", 1.5, None, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile="high_sensitivity", search_result_order="canonical",
            search_candidate_policy=policy, metrics_json=str(metrics))
    assert captured.get("candidate_policy", "single") == policy
    assert ("candidate_policy" in captured) == (policy != "single")
    assert json.loads(metrics.read_text())["metadata"]["search_candidate_policy"] == policy


def test_worker_handoff_uses_explicit_eleven_fields_only_for_union(tmp_path, monkeypatch):
    captured = []
    result = engine.IndexedSpeciesPairResults("A", "A", np.array([], dtype=np.int32), np.array([], dtype=np.int32),
        np.array([]), np.array([]), 0)
    class Pool:
        def __init__(self, **kwargs): pass
        def imap_unordered(self, function, args, chunksize):
            captured.extend(args)
            return iter([(result, ["gene"])])
        def close(self): pass
        def join(self): pass
    monkeypatch.setattr(engine.multiprocessing, "Pool", Pool)
    engine.execute_builtin_search(["A"], str(tmp_path), str(tmp_path), 2, SubstitutionMatrix.blosum62,
        max_candidates_per_query=100, candidate_policy="protected_union")
    assert len(captured[0]) == 11 and captured[0][-2:] == (4, "protected_union")
    options = {}
    monkeypatch.setattr(engine.SpeciesSequences, "from_fasta", lambda *_: SimpleNamespace(ids=["g"]))
    def search(*args, **kwargs):
        options.update(kwargs)
        return SimpleNamespace(evalues=np.array([0.]))
    monkeypatch.setattr(engine, "search_species_pair_indexed", search)
    engine._search_pair_worker(captured[0])
    assert options["candidate_policy"] == "protected_union" and "min_total_hits" not in options


def test_checkpoint_policy_isolated_and_default_manifest_unchanged(tmp_path):
    checkpoints = []
    for label, options in (("implicit", {}), ("single", {"search_candidate_policy": "single"}),
                            ("union", {"search_candidate_policy": "protected_union"})):
        checkpoints.append(write_accuracy_checkpoint(str(tmp_path / label), ["a", "b"], [0, 1], [0], [1], [1.], **options))
    assert (checkpoints[0] / "manifest.json").read_bytes() == (checkpoints[1] / "manifest.json").read_bytes()
    assert "search_candidate_policy" not in json.loads((checkpoints[0] / "manifest.json").read_text())
    with pytest.raises(ValueError, match="candidate policy mismatch"):
        load_accuracy_checkpoint(checkpoints[2])
    with pytest.raises(ValueError, match="candidate policy mismatch"):
        load_accuracy_checkpoint(checkpoints[0], expected_candidate_policy="protected_union")
    assert load_accuracy_checkpoint(checkpoints[2], expected_candidate_policy="protected_union")[0] == ["a", "b"]
    from benchmark_tools.audit_accuracy_checkpoint import audit
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    pin = record(checkpoints[2] / "manifest.json")
    with pytest.raises(ValueError, match="candidate policy mismatch"):
        audit(checkpoints[2], pin["sha256"])
    assert audit(checkpoints[2], pin["sha256"], expected_candidate_policy="protected_union")["summary"]["genes"] == 2


def test_real_cpu_union_scores_match_independent_single_searches(tmp_path, monkeypatch):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from orthohmm.search.prefilter import is_c_prefilter_available
    if not engine.is_c_available() or not is_c_prefilter_available():
        pytest.skip("Requires compiled CPU HMM and prefilter")
    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    path = tmp_path / "proteins.faa"
    path.write_text("".join(f">{gene}\n{seq}\n" for genes in fixture().values() for gene, seq in genes.items()))
    species = engine.SpeciesSequences.from_fasta(str(path), path.name)
    def search(**options):
        return engine.search_species_pair_indexed(species, species, "BLOSUM62", max_candidates_per_query=3, **options)
    def values(result):
        return {(int(q), int(t)): (np.float64(s).tobytes(), np.float64(e).tobytes())
                for q, t, s, e in zip(result.query_indices, result.target_indices, result.scores, result.evalues)}
    left, right = search(n_threads=1), search(n_threads=1, min_total_hits=1)
    a, b = values(left), values(right)
    assert all(a[k] == b[k] for k in a.keys() & b.keys())
    for threads in (1, 2):
        union = search(n_threads=threads, candidate_policy="protected_union")
        assert values(union) == {**a, **b}
        assert list(zip(union.query_indices, union.target_indices)) == sorted(a.keys() | b.keys())
        assert union.candidate_count == len(a.keys() | b.keys())


def test_real_native_pipeline_records_and_validates_union_without_losing_gene_universe(tmp_path):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from benchmark_tools.run_full_search_hillclimb import native_command, validate_native, PRODUCTS
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    source = Path(__file__).resolve().parents[2]
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools): pytest.skip("Requires MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    folder = tmp_path / "native"
    folder.mkdir()
    command = native_command(source, fasta, folder, tools, 3, 2)
    command.extend(["--search-result-order", "canonical", "--search-candidate-policy", "protected_union"])
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr + (folder / "output/benchmark.log").read_text()
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "search_max_candidates_per_query": 3,
        "search_candidate_policy": "protected_union", "inputs": [record(p) for p in sorted(fasta.iterdir())],
        "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    native = validate_native(case, 2)
    assert (native["genes"], native["species"]) == (16, 4)
    assert native["metrics"]["metadata"]["search_candidate_policy"] == "protected_union"
    case.pop("search_candidate_policy")
    with pytest.raises(ValueError, match="search candidate policy"):
        validate_native(case, 2)
