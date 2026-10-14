import json
from pathlib import Path
import os
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.accuracy import resolve_search_min_total_hits
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.parser import create_parser
from orthohmm.search import engine


@pytest.mark.parametrize("value", [1, 4, np.int64(1), 2**31 - 1])
def test_positive_int32_thresholds(value):
    assert resolve_search_min_total_hits(value) == int(value)


@pytest.mark.parametrize("value", [0, -1, 2**31, True, np.bool_(False), None, 1., "1", float("nan")])
def test_invalid_thresholds(value):
    with pytest.raises(ValueError, match="positive int32"):
        resolve_search_min_total_hits(value)


def test_cli_keeps_default_dictionary_shape(tmp_path):
    parser = create_parser()
    default = process_args(parser.parse_args([str(tmp_path)]))
    explicit = process_args(parser.parse_args([str(tmp_path), "--search_min_total_hits", "4"]))
    assert "search_min_total_hits" not in default and default == explicit
    configured = process_args(parser.parse_args([str(tmp_path), "--accuracy_profile", "high_sensitivity", "--search_min_total_hits", "1"]))
    assert configured["search_min_total_hits"] == 1


@pytest.mark.parametrize("extra", [[], ["--accuracy_profile", "high_sensitivity", "--search_max_candidates_per_query", "0"],
    ["--accuracy_profile", "high_sensitivity", "--start", "search_res"],
    ["--accuracy_profile", "high_sensitivity", "--search_mode", "phmmer"]])
def test_cli_rejects_inactive_or_adaptive_interventions(tmp_path, monkeypatch, extra):
    import orthohmm.args_processing as processing
    monkeypatch.setattr(processing, "find_executable", lambda name: name)
    with pytest.raises(SystemExit) as error:
        process_args(create_parser().parse_args([str(tmp_path), "--search_min_total_hits", "1", *extra]))
    assert error.value.code == 1


@pytest.mark.parametrize("profile,search,start,cap", [("standard", "builtin", None, None),
    ("high_sensitivity", "phmmer", None, None), ("high_sensitivity", "builtin", StartStep.search_res, None),
    ("high_sensitivity", "builtin", None, 0)])
def test_api_rejects_inactive_interventions_before_search(tmp_path, profile, search, start, cap):
    from orthohmm.orthohmm import execute
    with pytest.raises(ValueError, match="fixed positive cap"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 2, .5, "mcl", 1.5, start, None,
            SubstitutionMatrix.blosum62, 1e-4, accuracy_profile=profile, search_mode=search,
            search_min_total_hits=1, **({"search_max_candidates_per_query": cap} if cap is not None else {}))


@pytest.mark.parametrize("value", [4, 1])
def test_pipeline_records_actual_value_and_forwards_only_nondefault(tmp_path, monkeypatch, value):
    from orthohmm.orthohmm import execute
    class SearchReached(Exception): pass
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
            search_min_total_hits=value, metrics_json=str(metrics))
    assert captured.get("min_total_hits", 4) == value
    assert ("min_total_hits" in captured) == (value != 4)
    assert json.loads(metrics.read_text())["metadata"]["search_min_total_hits"] == value


@pytest.mark.parametrize("value", [4, 1])
def test_pool_worker_tuple_is_unchanged_by_default(tmp_path, monkeypatch, value):
    captured = []
    result = engine.IndexedSpeciesPairResults("A", "A", np.array([], dtype=np.int32), np.array([], dtype=np.int32), np.array([]), np.array([]), 0)
    class Pool:
        def __init__(self, **kwargs): pass
        def imap_unordered(self, function, args, chunksize):
            captured.extend(args)
            return iter([(result, ["gene"])])
        def close(self): pass
        def join(self): pass
    monkeypatch.setattr(engine.multiprocessing, "Pool", Pool)
    engine.execute_builtin_search(["A"], str(tmp_path), str(tmp_path), 2,
        SubstitutionMatrix.blosum62, max_candidates_per_query=100, min_total_hits=value)
    assert len(captured[0]) == (9 if value == 4 else 10)
    if value == 1:
        assert captured[0][-1] == 1


@pytest.mark.parametrize("extra", [(), (4,), (1,)])
def test_worker_forwards_nondefault_and_preserves_default_pair_call(monkeypatch, extra):
    captured = {}
    monkeypatch.setattr(engine.SpeciesSequences, "from_fasta", lambda *_: SimpleNamespace(ids=["g"]))
    def search(*args, **kwargs):
        captured.update(kwargs)
        return SimpleNamespace(evalues=np.array([0.]))
    monkeypatch.setattr(engine, "search_species_pair_indexed", search)
    engine._search_pair_worker(("A", "A", "/input", "BLOSUM62", 64, 4, 100, 1e-4, 1) + extra)
    assert captured.get("min_total_hits", 4) == (extra[0] if extra else 4)
    assert ("min_total_hits" in captured) == (extra == (1,))


def test_real_fresh_native_pipeline_records_the_nondefault_threshold(tmp_path):
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
    command = native_command(source, fasta, folder, tools, 100, 2)
    command.extend(["--search-result-order", "canonical", "--search-min-total-hits", "1"])
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr + (folder / "output/benchmark.log").read_text()
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "search_max_candidates_per_query": 100,
            "search_min_total_hits": 1, "inputs": [record(p) for p in sorted(fasta.iterdir())],
            "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    actual = validate_native(case, 2)
    assert actual["metrics"]["metadata"]["search_min_total_hits"] == 1
    assert (actual["genes"], actual["species"]) == (16, 4)
    case["search_min_total_hits"] = 4
    with pytest.raises(ValueError, match="shared-hit threshold"):
        validate_native(case, 2)
