import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_boundary_clear_pipeline import route_totals, compare_products, product_records
from benchmark_tools.run_full_search_hillclimb import native_command, validate_native, PRODUCTS
from benchmark_tools.strict_native_cli import install_guard


def fake_modules(tmp_path, scores=None, cuda=False):
    library_path = tmp_path / "private.so"
    library_path.write_bytes(b"test-only fake library identity")
    lib = SimpleNamespace(_name=str(library_path), batch_hmm_viterbi_multipair_avx2_c=True)
    vector = np.array([17, 23], dtype=np.int32) if scores is None else scores
    def score(*args, **kwargs):
        return vector
    def module():
        return SimpleNamespace(batch_viterbi_multipair_c=score, is_cuda_available=lambda: cuda,
            batch_viterbi_c=score, batch_viterbi_score=score, batch_viterbi_cuda=score)
    return module(), module(), SimpleNamespace(_load_c_hmm=lambda: (lib, True)), vector


def test_guard_records_actual_routes_without_changing_scores(tmp_path):
    engine, profiles, viterbi, scores = fake_modules(tmp_path)
    routes = tmp_path / "routes"
    pin = install_guard(engine, profiles, viterbi, routes)
    args = [None] * 8 + [np.zeros((2, 2), dtype=np.int32), 64]
    assert engine.batch_viterbi_multipair_c(*args) is scores
    assert profiles.batch_viterbi_multipair_c(*args) is scores
    totals = route_totals(routes, pin, {"search_candidates": 2, "high_sensitivity_profile_candidates": 2})
    assert totals["totals"] == {name: {"calls": 1, "scored_pairs": 2} for name in ("initial_search", "profile_scoring")}
    assert len(totals["files"]) == 1
    for name in ("batch_viterbi_c", "batch_viterbi_score", "batch_viterbi_cuda"):
        with pytest.raises(AssertionError, match="fallback"):
            getattr(engine, name)(*args)
        with pytest.raises(AssertionError, match="fallback"):
            getattr(profiles, name)(*args)


@pytest.mark.parametrize("scores", [np.array([1., 2.]), np.array([1], dtype=np.int32), np.array([[1, 2]], dtype=np.int32)])
def test_guard_rejects_invalid_native_return_arrays(tmp_path, scores):
    engine, profiles, viterbi, _ = fake_modules(tmp_path, scores)
    routes = tmp_path / "routes"
    install_guard(engine, profiles, viterbi, routes)
    with pytest.raises(ValueError, match="score array"):
        engine.batch_viterbi_multipair_c(*([None] * 8 + [np.zeros((2, 2)), 64]))
    assert not list(routes.iterdir())


def test_guard_rejects_cuda_before_creating_evidence(tmp_path):
    engine, profiles, viterbi, _ = fake_modules(tmp_path, cuda=True)
    with pytest.raises(ValueError, match="hidden CUDA"):
        install_guard(engine, profiles, viterbi, tmp_path / "routes")
    assert not (tmp_path / "routes").exists()


@pytest.mark.parametrize("counts", [
    {"search_candidates": 3, "high_sensitivity_profile_candidates": 2},
    {"search_candidates": 2, "high_sensitivity_profile_candidates": 3},
])
def test_route_counts_must_cover_actual_pipeline_candidates(tmp_path, counts):
    engine, profiles, viterbi, _ = fake_modules(tmp_path)
    routes = tmp_path / "routes"
    pin = install_guard(engine, profiles, viterbi, routes)
    args = [None] * 8 + [np.zeros((2, 2)), 64]
    engine.batch_viterbi_multipair_c(*args)
    profiles.batch_viterbi_multipair_c(*args)
    with pytest.raises(ValueError, match="cover"):
        route_totals(routes, pin, counts)


def test_native_equivalence_requires_inventory_and_content(tmp_path):
    path = tmp_path / "a"
    path.write_text("one")
    pin = record(path)
    left = {"orthogroups": pin}
    right = {"orthogroups": {**pin, "path": "/other"}}
    assert compare_products(left, right) == {"files": 1, "all_hashes_equal": True}
    with pytest.raises(ValueError, match="inventories"):
        compare_products(left, {})
    with pytest.raises(ValueError, match="equivalence"):
        compare_products(left, {"orthogroups": {**pin, "sha256": "changed"}})


def test_real_guarded_cli_and_full_native_reader_on_small_fixture(tmp_path):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    source = Path(__file__).resolve().parents[2]
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools):
        pytest.skip("Requires MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{sequence}\n" for gene, sequence in genes.items()))
    folder = tmp_path / "native"
    folder.mkdir()
    trace = folder / "routes"
    command = native_command(source, fasta, folder, tools, 100, 2)
    command.extend(["--search-result-order", "canonical", "--native-kernel-trace", str(trace)])
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "NUMBA_NUM_THREADS": "2", "PYTHONPATH": str(source)}
    actual = subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert actual.returncode == 0, actual.stdout + actual.stderr + (folder / "output/benchmark.log").read_text()
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"output": str(folder), "fasta": str(fasta), "search_max_candidates_per_query": 100,
            "inputs": [record(p) for p in sorted(fasta.iterdir())],
            "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    result = validate_native(case, 2)
    assert (result["genes"], result["species"]) == (16, 4)
    totals = route_totals(trace, record(source / "orthohmm/search/csrc/hmm_viterbi.so"), result["metrics"]["counts"])
    assert totals["totals"]["initial_search"]["calls"] >= 16
    assert totals["totals"]["profile_scoring"]["calls"] > 0
    assert result["metrics"]["metadata"]["search_result_order"] == "canonical"
    assert result["metrics"]["harness"]["command"][1].endswith("strict_native_cli.py")
    products = product_records(folder / "output")
    assert compare_products(products, products)["all_hashes_equal"]


def test_failed_full_clear_control_does_not_launch_optimized_arm(tmp_path, monkeypatch):
    import benchmark_tools.run_boundary_clear_pipeline as runner
    from benchmark_tools.run_membership_hillclimb import write_json
    source = Path(runner.__file__).resolve().parent.parent
    output = tmp_path / "output"
    output.mkdir()
    plan = output / "execution_plan.json"
    write_json(plan, {"repo": str(tmp_path), "source_root": str(source), "output": str(output),
                     "source_revision": "mocked", "inputs": [], "sources": [], "runtime": {"files": []},
                     "scope": "mocked execution", "tool_path": "/unused",
                     "cases": [{"label": name, "source_root": str(source), "output": str(output / name),
                                "command": [name]} for name in ("full_clear", "boundary_only")]})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: "mocked host")
    monkeypatch.setattr(runner, "verify_input_inventory", lambda *_: None)
    calls = []
    monkeypatch.setattr(runner, "timed_stage", lambda command, *args: calls.append(command) or {"returncode": 1})
    result = runner.run(plan)
    assert calls == [["full_clear"]]
    assert result["status"] == "finished_with_explicit_gaps"
    assert [row["status"] for row in result["arms"]] == ["failed", "not_run_failed_control"]
    assert not (output / "boundary_only").exists()
