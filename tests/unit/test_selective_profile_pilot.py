from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools import run_profile_weighting_pilot as pilot
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_selective_profile_evidence import validate
from orthohmm.search.selective_profiles import selective_profile_edges
from tests.unit.test_profile_weighting_pilot import fake_execution
from tests.unit.test_selective_profiles import real_stage


ARMS = ("off", "builtin_reselect", "hmmer_score_only", "hmmer_reselect")


@pytest.mark.parametrize("mode", ARMS)
def test_only_selective_graph_flag_changes_and_native_policies_stay_fixed(mode):
    source, folder = Path("/frozen"), Path("/out")
    plan = {"cpu": 16, "experiment": "selective_hmmer", "case": {
        "checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"}, "fasta_directory": "/proteomes",
        "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/old/native"}}
    graph, native = pilot.arm_commands(source, plan, folder, mode)
    off_graph, off_native = pilot.arm_commands(source, plan, folder, "off")
    assert native == off_native
    assert graph == (off_graph if mode == "off" else off_graph + ["--selective-profile-mode", mode])
    assert pilot.experiment_settings("selective_hmmer")[0] == ARMS


def fake_selective_execution(tmp_path, monkeypatch, *, failed_stage=None, changed_digest=False, changed_graph=False):
    path, executions = fake_execution(tmp_path, monkeypatch, failed_stage=failed_stage)
    plan = json.loads(path.read_text())
    plan.update(experiment="selective_hmmer", arms=list(ARMS), hmmer_tools={})
    plan["case"]["baseline_output"] = "/baseline"
    path.write_text(json.dumps(plan))
    def activation(graph, mode, experiment, **kwargs):
        assert experiment == "selective_hmmer"
        return {"candidate_digest": "changed" if changed_digest and mode == "hmmer_score_only" else "same",
                "completely_measured_queries": 2, "applicability": "completely_measured"}
    monkeypatch.setattr(pilot, "validate_activation", activation)
    monkeypatch.setattr(pilot, "native_pair_turnover", lambda *a: {"added_pairs": 0, "lost_pairs": 0})
    if changed_graph:
        original = pilot.timed_stage
        def timed(command, folder, name, source, env):
            result = original(command, folder, name, source, env)
            if command[0] == "hmmer_score_only" and name == "graph":
                graph = json.loads((folder / "graph_metrics.json").read_text())
                graph["stages"] = [{"label": "changed", "output": {"bytes": 1, "sha256": "changed"}}]
                (folder / "graph_metrics.json").write_text(json.dumps(graph))
            return result
        monkeypatch.setattr(pilot, "timed_stage", timed)
    return path, executions


def test_four_arm_runner_requires_exact_measurement_control_and_scores_both_reselects(tmp_path, monkeypatch):
    path, executions = fake_selective_execution(tmp_path, monkeypatch)
    result = pilot.run(path)
    assert result["status"] == "finished_all_arms_verified"
    assert executions == [(mode, stage) for mode in ARMS for stage in ("graph", "native")]
    assert result["arms"][2]["measurement_control_matches"]
    assert not result["arms"][1]["advancement"]["eligible"]
    assert not result["arms"][3]["advancement"]["eligible"]
    assert "advancement" not in result["arms"][2]


@pytest.mark.parametrize("kind", ["graph_failure", "native_failure", "candidate_digest", "graph_products"])
def test_measurement_control_failure_blocks_active_hmmer_without_restarting_or_fallback(tmp_path, monkeypatch, kind):
    failed = ("hmmer_score_only", "graph" if kind == "graph_failure" else "native") if kind.endswith("failure") else None
    path, executions = fake_selective_execution(tmp_path, monkeypatch, failed_stage=failed,
        changed_digest=kind == "candidate_digest", changed_graph=kind == "graph_products")
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert result["arms"][1]["status"] == "independent_content_readback_verified"
    assert result["arms"][2]["status"] == "failed"
    assert result["arms"][3]["status"] == "not_run_failed_measurement_control_gate"
    assert not any(mode == "hmmer_reselect" for mode, _ in executions)


@pytest.mark.parametrize("mode", ["builtin_reselect", "hmmer_score_only", "hmmer_reselect"])
def test_raw_measurement_readback_reproduces_actual_nonempty_stage(real_stage, tmp_path, mode):
    args, database = real_stage
    graphdir = tmp_path / "graph"
    graphdir.mkdir()
    (graphdir / "orthogroups_multipass.txt").write_text("".join(" ".join(args[0][g] for g in group) + "\n" for group in args[1]))
    edges, evidence = selective_profile_edges(*args, database, mode=mode, output_directory=graphdir / "selective_profiles", cpu=2)
    graph = {"parameters": {"selective_profile_mode": mode}, "counts": {"strict_profile_edges": len(edges)},
             "selective_profile_evidence": [evidence]}
    actual = validate(graph, mode, tmp_path, args[0])
    assert actual["selected_queries"] == actual["completely_measured_queries"] == 1
    assert actual["selected_profiles"] == 2 and actual["checked"]


@pytest.mark.parametrize("kind", ["receipt", "digest", "source", "measurements", "decision", "changes", "database", "own"])
def test_raw_readback_rejects_tampering_even_when_graph_claims_success(real_stage, tmp_path, kind):
    args, database = real_stage
    graphdir = tmp_path / "graph"
    graphdir.mkdir()
    partition = graphdir / "orthogroups_multipass.txt"
    partition.write_text("".join(" ".join(args[0][g] for g in group) + "\n" for group in args[1]))
    edges, evidence = selective_profile_edges(*args, database, mode="builtin_reselect",
        output_directory=graphdir / "selective_profiles", cpu=1)
    evidence = deepcopy(evidence)
    graph = {"parameters": {"selective_profile_mode": "builtin_reselect"}, "counts": {"strict_profile_edges": len(edges)},
             "selective_profile_evidence": [evidence]}
    receipt = graphdir / "selective_profiles/selection.json"
    if kind == "receipt":
        receipt.write_text("{}")
    elif kind == "digest":
        evidence["selection"]["candidate_digest"] = "wrong"
    elif kind == "source":
        lines = partition.read_text().splitlines()
        partition.write_text("\n".join([lines[1], lines[0], *lines[2:]]) + "\n")
    elif kind in ("measurements", "database"):
        path = graphdir / "selective_profiles/builtin_measurements.json"
        data = json.loads(path.read_text())
        if kind == "measurements":
            data["0"]["measurements"]["query"]["full_sequence_score"] += 100
        else:
            data["0"]["parameters"]["database_genes"] = 4
        path.write_text(json.dumps(data))
        if kind == "database":
            evidence["measurements"] = {p: record(path) for p in evidence["measurements"]}
    elif kind == "decision":
        evidence["decisions"][0]["action"] = "reject"
    elif kind == "changes":
        evidence["applied_changed_query_anchors"] = []
    else:
        evidence["selection"]["selected"][0]["own_profile_id"] = 1
    if kind != "receipt":
        receipt.write_text(json.dumps(evidence))
    with pytest.raises(ValueError):
        validate(graph, "builtin_reselect", tmp_path, args[0])


def test_off_readback_needs_no_selective_files_and_rejects_hidden_intervention(tmp_path):
    graph = {"parameters": {}}
    assert validate(graph, "off", tmp_path, []) == {"selective_mode": "off"}
    graph["selective_profile_evidence"] = [{}]
    with pytest.raises(ValueError, match="Off control"):
        validate(graph, "off", tmp_path, [])


@pytest.mark.parametrize("experiment", ["selective_hmmer", "projection_recovery"])
def test_actual_four_arm_graph_candidates_and_inferred_native_control(tmp_path, monkeypatch, experiment):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from benchmark_tools.run_full_search_hillclimb import native_command
    tools = [shutil.which(tool) for tool in ("mafft", "FastTree", "hmmbuild", "hmmsearch")]
    if not all(tools):
        pytest.skip("Requires actual MAFFT/FastTree/HMMER3.4")
    source = Path(pilot.__file__).resolve().parents[1]
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{g}\n{s}\n" for g, s in genes.items()))
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    command = native_command(source, fasta, baseline, tools[:2], 100, 2) + ["--search-result-order", "canonical"]
    environment = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    done = subprocess.run(command, cwd=source, env=environment, capture_output=True, text=True, timeout=180)
    assert done.returncode == 0, done.stdout + done.stderr
    native = baseline / "output"
    phylo = native / "orthohmm_phylogeny"
    work = native / "orthohmm_working_res"
    checkpoint = work / "high_sensitivity_checkpoint"
    output = tmp_path / "pilot"
    output.mkdir()
    case = {"dataset": "OrthoBench", "fasta_directory": str(fasta), "baseline_output": str(native),
        "checkpoint": str(checkpoint), "checkpoint_manifest": record(checkpoint / "manifest.json"),
        "fastas": [record(p) for p in sorted(fasta.iterdir())],
        "candidate": record(work / "phylogeny_candidate_superfamilies.txt"),
        "constraints": record(work / "phylogeny_candidate_merges.json"), "aligner": tools[0], "tree_builder": tools[1],
        "golden": {name: record(phylo / name) for name in (*pilot.PRODUCTS, "species_tree.rooted.nwk")},
        "expected": {"f1": .7, "precision": .8, "recall": .6}}
    arms = pilot.experiment_settings(experiment)[0]
    plan = {"source_root": str(source), "repo": str(source), "output": str(output), "cpu": 16,
        "experiment": experiment, "arms": list(arms), "sources": [], "inputs": [],
        "runtime": {"files": []}, "source_revision": "unit_fixture_not_frozen_biological_run", "tool_path": os.environ["PATH"],
        "scope": "Actual complete16-gene/four-species integration, scorer mocked; no biological accuracy claim",
        "hmmer_tools": {name: record(path) for name, path in zip(("hmmbuild", "hmmsearch"), tools[2:])}, "case": case}
    plan_path = output / "execution_plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(pilot, "UNIVERSE", (16, 4))
    monkeypatch.setattr(pilot, "runtime", lambda: {"files": []})
    # This synthetic fixture tests real inference/writers, not benchmark truth.
    monkeypatch.setattr(pilot, "score_orthobench", lambda *a: {"scores": {"OrthoBench F1": .7}})
    monkeypatch.setattr(pilot, "readback_scores", lambda *a: (
        {"OrthoBench F1": .7}, {"precision": .8, "recall": .6}, []))
    result = pilot.run(plan_path)
    assert result["status"] == "finished_all_arms_verified", result["arms"]
    for row in result["arms"]:
        assert row["native_summary"]["candidate_families"] > 0
        assert row["coverage"]["all_candidate_genes_preserved"]
    measurement = result["arms"][2 if experiment == "selective_hmmer" else 3]
    assert measurement["measurement_control_matches"]
    assert all(measurement["golden_matches"].values())
    # Separate real profile-stage fixtures require nonempty HMMER activation;
    # this complete native fixture can be inactive and remains so labeled.
    assert result["arms"][-1]["activation"]["applicability"] in ("completely_measured", "no_complete_measurements")
