import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import run_full_search_hillclimb as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_input_inventory_matches_production_reader_not_staging_metadata(tmp_path):
    (tmp_path / "species.fasta").write_text(">gene\nACDE\n")
    (tmp_path / "staging_manifest.json").write_text("{}")
    inputs = runner.input_records(tmp_path, 1)
    assert [Path(item["path"]).name for item in inputs] == ["species.fasta"]
    case = {"fasta": str(tmp_path), "inputs": inputs}
    runner.verify_input_inventory(case)
    (tmp_path / "extra.faa").write_text(">foreign\nACDE\n")
    with pytest.raises(ValueError, match="inventory changed"):
        runner.verify_input_inventory(case)
    with pytest.raises(ValueError, match="Unexpected input proteome inventory"):
        runner.input_records(tmp_path, 1)


def test_command_changes_only_initial_cap(tmp_path):
    tools = ["mafft", "FastTree"]
    control = runner.native_command(tmp_path, tmp_path, tmp_path, tools, 100, 32)
    variant = runner.native_command(tmp_path, tmp_path, tmp_path, tools, 400, 32)
    index = control.index("--search-max-candidates-per-query") + 1
    assert variant[:index] == control[:index]
    assert variant[index + 1:] == control[index + 1:]
    assert (control[index], variant[index]) == ("100", "400")
    assert "--full-output" in control
    assert "--checkpoint-source" not in control
    assert control[control.index("--species-tree-mode") + 1] == "infer"


@pytest.mark.parametrize("body", [
    "a\tA\tb\tB\thigh\n", "a\tA\tb\tB\tmedium\n", "a\tA\tb\tB\tlow\n",
])
def test_matching_confidence_rows(tmp_path, body):
    pairs = tmp_path / "pairs.tsv"
    confidence = tmp_path / "confidence.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tB\n")
    confidence.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n" + body)
    assert runner.validate_confidence(confidence, pairs) == 1


@pytest.mark.parametrize("body", [
    "", "a\tA\tc\tB\thigh\n", "a\tA\tb\tB\tunknown\n",
    "a\tA\tb\tB\thigh\nextra\n", "a\tA\tb\tB\n",
])
def test_confidence_differences_fail(tmp_path, body):
    pairs = tmp_path / "pairs.tsv"
    confidence = tmp_path / "confidence.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tB\n")
    confidence.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n" + body)
    with pytest.raises(ValueError, match="outputs differ"):
        runner.validate_confidence(confidence, pairs)


def fake_plan(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    plan = {"repo": str(tmp_path), "source_root": str(source), "source_revision": "fixture",
        "output": str(tmp_path / "results"), "cpu": 32, "tool_path": os.environ["PATH"],
        "scope": "test fixture", "runtime": {"files": []}, "cases": []}
    Path(plan["output"]).mkdir()
    for label in ("control100", "cap400"):
        folder = Path(plan["output"]) / "OrthoBench" / label
        plan["cases"].append({"dataset": "OrthoBench", "label": label,
            "output": str(folder), "command": [label], "fasta": str(tmp_path), "inputs": []})
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setattr(runner, "checked_records", lambda *args: [])
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: "fixture activity")
    monkeypatch.setattr(runner, "require_profile_runtime", lambda *args: {"status": "passed"})
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *args, **kwargs: "")
    monkeypatch.setattr(runner, "validate_native", lambda *args: {"status": "fixture"})
    monkeypatch.setattr(runner, "score_orthobench", lambda *args: {"scores": {"OrthoBench F1": .75}})
    return path


def test_runner_preserves_control_failure_and_skips_variant(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    calls = []
    def child(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=1)
    monkeypatch.setattr(runner.subprocess, "run", child)
    result = runner.run(path, "OrthoBench")
    assert calls == [["control100"]]
    assert [r["status"] for r in result["arms"]] == ["failed", "not_run_failed_control"]
    assert result["status"] == "finished_with_explicit_gaps"
    with pytest.raises(FileExistsError):
        runner.run(path, "OrthoBench")


def test_runner_requires_scored_control_before_variant(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", lambda command, **kwargs:
                        calls.append(command) or SimpleNamespace(returncode=0))
    def fail(*args):
        raise ValueError("scoring failed")
    monkeypatch.setattr(runner, "score_orthobench", fail)
    result = runner.run(path, "OrthoBench")
    assert calls == [["control100"]]
    assert result["arms"][0]["error"] == "scoring failed"
    assert result["arms"][1]["status"] == "not_run_failed_control"


def test_runner_completed_cells_and_retry_gate(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    result = runner.run(path, "OrthoBench")
    assert result["status"] == "finished_all_arms_scored"
    assert [r["scores"] for r in result["arms"]] == [{"OrthoBench F1": .75}] * 2
    assert all(r["inputs_sources_runtime_reverified"] for r in result["arms"])
    assert result["defaults_promoted"] is False
    with pytest.raises(FileExistsError):
        runner.run(path, "OrthoBench")


def test_runner_runtime_mutation_invalidates_score_and_gates_variant(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    def score(*args):
        monkeypatch.setattr(runner, "runtime", lambda: {"files": [], "changed": True})
        return {"scores": {"OrthoBench F1": .75}}
    monkeypatch.setattr(runner, "score_orthobench", score)
    result = runner.run(path, "OrthoBench")
    assert all(r["status"] == "failed" for r in result["arms"])
    assert "Runtime changed during run" in result["arms"][0]["error"]
    assert result["status"] == "finished_with_explicit_gaps"


def test_real_small_full_pipeline_and_independent_native_validation(tmp_path):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    source = Path(runner.__file__).resolve().parent.parent
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools):
        pytest.skip("Requires installed MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{g}\n{s}\n" for g, s in genes.items()))
    folder = tmp_path / "control"
    folder.mkdir()
    command = runner.native_command(source, fasta, folder, tools, 100, 2)
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "NUMBA_NUM_THREADS": "2"}
    child = runner.subprocess.run(command, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"output": str(folder), "fasta": str(fasta), "search_max_candidates_per_query": 100,
        "inputs": [record(p) for p in sorted(fasta.iterdir())],
        "golden": {name: record(phylo / name) for name in (*runner.PRODUCTS, "species_tree.rooted.nwk")}}
    result = runner.validate_native(case, 2)
    assert (result["genes"], result["species"]) == (16, 4)
    assert result["pairs"] > 0
    assert all(result["historical_product_matches"].values())
    assert result["summary"]["checkpoint_hits"] == 0
    from benchmark_tools.shared_qfo_export import normalize, verify_projection
    pairs, confidence = phylo / runner.PRODUCTS[1], phylo / runner.PRODUCTS[2]
    pair_bytes, confidence_bytes = pairs.read_bytes(), confidence.read_bytes()
    pair_rows, confidence_rows = pair_bytes.splitlines(keepends=True), confidence_bytes.splitlines(keepends=True)
    pairs.write_bytes(pair_rows[0] + pair_rows[1] + b"".join(pair_rows[1:]))
    confidence.write_bytes(confidence_rows[0] + confidence_rows[1] + b"".join(confidence_rows[1:]))
    invalid_pins = record(pairs), record(confidence)
    with pytest.raises(ValueError, match="canonical, unique and sorted"):
        runner.validate_native(case, 2)
    proof = normalize(pairs, confidence, tmp_path / "canonical_pairs", tmp_path / "canonical_confidence",
        tuple(pair_rows[1].decode().rstrip("\n").split("\t")), len(pair_rows), len(pair_rows) - 1)
    assert verify_projection(proof)["distinct_rows"] == result["pairs"]
    normalized = runner.validate_native(case, 2, pair_export=proof)
    assert normalized["status"] == "derived_native_content_verified"
    assert normalized["pairs"] == result["pairs"] and normalized["genes"] == result["genes"]
    assert normalized["products"][runner.PRODUCTS[1]]["sha256"] == case["golden"][runner.PRODUCTS[1]]["sha256"]
    assert (record(pairs), record(confidence)) == invalid_pins
    with pytest.raises(ValueError, match="canonical, unique and sorted"):
        runner.validate_native(case, 2)
    pairs.write_bytes(pair_bytes)
    confidence.write_bytes(confidence_bytes)
    metrics_path = folder / "metrics.json"
    data = json.loads(metrics_path.read_text())
    data["metadata"]["search_max_candidates_per_query"] = 400
    metrics_path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="parameters differ"):
        runner.validate_native(case, 2)
    data["metadata"]["search_max_candidates_per_query"] = 100
    metrics_path.write_text(json.dumps(data))
    manifest_path = phylo / "provenance_manifest.json"
    native = json.loads(manifest_path.read_text())
    native["input_proteomes"][0]["sha256"] = "foreign"
    manifest_path.write_text(json.dumps(native))
    with pytest.raises(ValueError, match="native species inventory"):
        runner.validate_native(case, 2)
