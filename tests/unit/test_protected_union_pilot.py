from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools import run_protected_union_pilot as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, write_json

ACCOUNTING = ("JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n"
    "24240|FAILED|1:0|16|bizon\n24360|COMPLETED|0:0|8|bizon\n24362|COMPLETED|0:0|8|bizon\n")


@pytest.fixture
def control():
    return {"scores": {"OrthoBench F1": .74}, "details": {"precision": .82, "recall": .68}}


@pytest.mark.parametrize("f1,precision,recall,passed", [(.75, .82, .68, True),
    (.74, .82, .68, False), (.73, .82, .68, False),
    (.75, .82 - .005, .68, True), (.75, .82 - .005001, .68, False),
    (.75, .82, .68 - .005, True), (.75, .82, .68 - .005001, False)])
def test_prespecified_gate_is_strict_and_not_retuned(control, f1, precision, recall, passed):
    actual = runner.pilot_gate({"OrthoBench F1": f1}, {"precision": precision, "recall": recall}, control)
    assert actual["passes_development_gate"] is passed
    assert actual["deltas"]["F1"] == f1 - .74


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), -1., 82.])
def test_gate_refuses_unormalized_or_nonfinite_axes(control, value):
    with pytest.raises(ValueError):
        runner.pilot_gate({"OrthoBench F1": .75}, {"precision": value, "recall": .68}, control)


def test_command_uses_fresh_full_inference_and_opt_in_policy(tmp_path):
    argv = runner.command(tmp_path / "executor", tmp_path / "input", tmp_path / "native",
                          ["/mafft", "/FastTree"], tmp_path / "routes")
    for flag, value in (("--cpu", "16"), ("--search-max-candidates-per-query", "100"),
        ("--search-candidate-policy", "protected_union"), ("--search-result-order", "canonical"),
        ("--phylogeny", "reconcile"), ("--species-tree-mode", "infer"), ("--refinement-profile", "default")):
        assert argument(argv, flag) == value
    assert "--full-output" in argv
    assert all(flag not in argv for flag in ("--start", "--checkpoint-source", "--accuracy-checkpoint"))


@pytest.fixture
def planned(tmp_path):
    paths = [tmp_path / "mafft", tmp_path / "FastTree"]
    for path in paths:
        path.write_text("trusted executable fixture\n")
    tools = [record(p) for p in paths]
    source, output, fasta = tmp_path / "executor", tmp_path / "output", tmp_path / "input"
    old = {"dataset": "OrthoBench", "fasta_directory": str(fasta), "fastas": [],
           "aligner": str(paths[0]), "tree_builder": str(paths[1])}
    case = {"dataset": "OrthoBench", "label": "protected_union", "output": str(output / "native"),
        "fasta": str(fasta), "inputs": [], "search_candidate_policy": "protected_union",
        "search_min_total_hits": 4, "search_max_candidates_per_query": 100,
        "command": runner.command(source, fasta, output / "native", [str(p) for p in paths], output / "native_routes")}
    plan = {"case": case, "output": str(output), "source_root": str(source), "tools": tools, "cpu": 16,
            "native_routes": str(output / "native_routes")}
    return plan, {"cases": [old], "inputs": tools}


def test_worker_plan_remains_bound_to_original_input_tools_and_namespace(planned):
    runner.validate_plan(*planned)


@pytest.mark.parametrize("field,value", [("dataset", "QfO"), ("label", "other"),
    ("fasta", "/different/proteomes"), ("inputs", [{"other": "proteome"}]),
    ("output", "/different/output"), ("command", ["other-command"]),
    ("search_candidate_policy", "single"), ("search_min_total_hits", 1),
    ("search_max_candidates_per_query", 400)])
def test_worker_rejects_changed_intervention_or_input(planned, field, value):
    plan, parent = planned
    plan["case"][field] = value
    with pytest.raises(ValueError, match="Pilot plan differs"):
        runner.validate_plan(plan, parent)


@pytest.mark.parametrize("field,value", [("cpu", 32), ("tools", []), ("native_routes", "/other/routes")])
def test_worker_rejects_changed_resource_tool_or_trace_plan(planned, field, value):
    plan, parent = planned
    plan[field] = value
    with pytest.raises(ValueError, match="Pilot plan differs"):
        runner.validate_plan(plan, parent)


def test_terminal_dependencies_preserve_failed_control_job_and_successful_admissions(monkeypatch):
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: ACCOUNTING)
    rows = runner.terminal_dependencies()
    assert rows[0]["State"] == "FAILED" and rows[-1]["JobIDRaw"] == "24362"


@pytest.mark.parametrize("text", [ACCOUNTING.replace("24362|COMPLETED", "24362|RUNNING"),
    ACCOUNTING.replace("24362|COMPLETED|0:0", "24362|FAILED|1:0"),
    ACCOUNTING.replace("24360|COMPLETED|0:0|8", "24360|COMPLETED|0:0|16"),
    ACCOUNTING.replace("bizon", "other"), ACCOUNTING + "24362|COMPLETED|0:0|8|bizon\n",
    ACCOUNTING.replace("24240|FAILED|1:0|16|bizon\n", "")])
def test_wrong_or_live_dependency_refused(monkeypatch, text):
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: text)
    with pytest.raises(ValueError, match="terminal original control"):
        runner.terminal_dependencies()


@pytest.fixture
def scientific(tmp_path):
    roots = tmp_path / "old", tmp_path / "new"
    names = ("orthohmm/__init__.py", "orthohmm/search/engine.py", "orthohmm/search/csrc/hmm_viterbi.c",
             "orthohmm/search/csrc/hmm_viterbi.so")
    pins = []
    for name in names:
        for root in roots:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((name + " fixture").encode())
        pins.append(record(roots[0] / name))
    return roots[1], {"source_root": str(roots[0]), "sources": pins}


def test_scientific_source_and_cpu_library_equivalence_required(scientific):
    source, fresh = scientific
    assert len(runner.scientific_source(source, fresh)) == 4


@pytest.mark.parametrize("change", ["python", "library", "extra", "empty"])
def test_changed_scientific_source_refused(scientific, change):
    source, fresh = scientific
    if change == "empty":
        fresh["sources"] = []
    else:
        name = {"python": "search/engine.py", "library": "search/csrc/hmm_viterbi.so", "extra": "new.py"}[change]
        (source / "orthohmm" / name).write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed after the fresh gate"):
        runner.scientific_source(source, fresh)


def test_prerequisite_receipt_must_match_frozen_git_content(tmp_path, monkeypatch):
    path = tmp_path / "receipt.json"
    write_json(path, {"status": "verified"})
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: path.read_bytes())
    value, pin = runner.committed(tmp_path, tmp_path / "executor", "revision", "receipt.json")
    assert value == {"status": "verified"} and pin == record(path)
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: b"wrong committed receipt")
    with pytest.raises(ValueError, match="committed admission"):
        runner.committed(tmp_path, tmp_path / "executor", "revision", "receipt.json")


@pytest.fixture
def prerequisites(tmp_path, monkeypatch):
    artifact = tmp_path / "artifact.json"
    write_json(artifact, {})
    pin = record(artifact)
    fresh_path = tmp_path / "fresh_plan.json"
    fresh = {"source_revision": "fresh", "inputs": [pin], "sources": [pin], "runtime": {"files": []}}
    write_json(fresh_path, fresh)
    fresh_pin = record(fresh_path)
    report_path = tmp_path / "fresh_report.json"
    report = {"status": runner.FRESH_SUCCESS, "postflight_pins_verified": True,
        "source_revision": "fresh", "execution_plan": fresh_pin, "fresh_initial_search_executed": True,
        "scored_candidates": 4560971, "actual_kernel_calls": 288,
        "directions": [{"control": {"raw": pin}, "protected_union": {"raw": pin}} for _ in range(144)]}
    write_json(report_path, report)
    gate = {"status": "fresh_protected_union_all_144_directions_independently_verified",
        "advance_to_orthobench_pilot": True, "directions": 144, "raw_outputs_independently_recounted": 288,
        "observed_cpu_scored_candidates": 4560971, "defaults_promoted": False, "source_revision": "fresh",
        "producer": record(report_path), "execution_plan": fresh_pin, "kernel_routes": pin, "time": pin}
    admission = {"status": "shared_qfo_all_six_scores_admitted_and_independently_crosschecked",
        "complete_eight_endpoint_contrast": True, "native_inference_rerun": False, "defaults_promoted": False,
        **{k: pin for k in ("summary", "score_table", "producer", "normalization", "reporter_tests")}}
    parent_path = tmp_path / "parent.json"
    parent = {"repo": str(tmp_path), "source_root": "old_source", "source_revision": "old",
              "runtime": {"files": []}, "sources": [pin]}
    write_json(parent_path, parent)
    launch = {"execution_plan": record(parent_path), "frozen_source": "old_source", "source_revision": "old"}
    values = {runner.FRESH_GATE: gate, runner.QFO_ADMISSION: admission, runner.CONTROL_LAUNCH: launch}
    monkeypatch.setattr(runner, "committed", lambda repo, source, revision, name: (deepcopy(values[name]), pin))
    monkeypatch.setattr(runner, "terminal_dependencies", lambda: [{"state": "terminal"}])
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "control_readback", lambda *a: {"coverage": {"genes": 251378, "species": 12}, "checked": []})
    return tmp_path, values, report_path


def test_fresh_equivalence_and_closed_qfo_required(prerequisites):
    repo, _, _ = prerequisites
    parent, fresh, control, pins, accounting = runner.prerequisites(repo, repo / "source", "revision")
    assert parent["source_revision"] == "old" and fresh["source_revision"] == "fresh"
    assert control["coverage"]["genes"] == 251378 and accounting and pins


@pytest.mark.parametrize("name,field,value", [(runner.FRESH_GATE, "advance_to_orthobench_pilot", False),
    (runner.FRESH_GATE, "directions", 143), (runner.FRESH_GATE, "raw_outputs_independently_recounted", 287),
    (runner.FRESH_GATE, "observed_cpu_scored_candidates", 1),
    (runner.QFO_ADMISSION, "complete_eight_endpoint_contrast", False),
    (runner.QFO_ADMISSION, "native_inference_rerun", True), (runner.QFO_ADMISSION, "defaults_promoted", True)])
def test_failed_or_incomplete_prerequisites_refused(prerequisites, name, field, value):
    repo, values, _ = prerequisites
    values[name][field] = value
    with pytest.raises(ValueError, match="not independently admitted"):
        runner.prerequisites(repo, repo / "source", "revision")


def test_changed_fresh_producer_refused(prerequisites):
    repo, _, path = prerequisites
    path.write_text("{}\n")
    with pytest.raises(ValueError, match="Frozen input/source identity changed"):
        runner.prerequisites(repo, repo / "source", "revision")


@pytest.fixture
def admitted_control(tmp_path, monkeypatch):
    folder = tmp_path / "control/OrthoBench/canonical"
    phylo = folder / "native/orthohmm_phylogeny"
    phylo.mkdir(parents=True)
    for name in (*PRODUCTS, "species_tree.rooted.nwk"):
        (phylo / name).write_text(name + "\n")
    write_json(folder / "preflight.json", {"native": ["native-command"]})
    products = {name: record(phylo / name) for name in PRODUCTS}
    native = {"genes": 3, "species": 2, "groups": 2, "pairs": 1, "products": products, "checked": [], "summary": {}}
    cell = {"dataset": "OrthoBench", "status": "independent_content_readback_verified", "products": products,
        "checked": [record(folder / "preflight.json"), *products.values()],
        "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "scores": {"OrthoBench: OrthoBench F1": .74}, "details": {"precision": .82, "recall": .68},
        "source_revision": "original", "resources": {"cached": True}}
    table = {"methods": [{"key": "canonical_collection:canonical", "cells": [cell]}]}
    monkeypatch.setattr(runner, "replay_native", lambda row, case: deepcopy(native))
    monkeypatch.setattr(runner, "readback_scores", lambda *a: ({"OrthoBench F1": .74}, deepcopy(cell["details"]), []))
    return table, {"output": str(tmp_path / "control"), "source_revision": "original"}, cell, native


def test_control_revalidated_without_new_inference(tmp_path, admitted_control):
    table, parent, cell, _ = admitted_control
    actual = runner.control_readback(tmp_path, table, parent)
    assert actual["scores"] == {"OrthoBench F1": .74}
    assert actual["coverage"] == cell["coverage"] and actual["resources"] == {"cached": True}


@pytest.mark.parametrize("change", ["coverage", "source", "duplicate", "preflight"])
def test_changed_or_ambiguous_control_refused(tmp_path, admitted_control, change):
    table, parent, cell, native = admitted_control
    if change == "coverage":
        native["pairs"] = 2
    elif change == "source":
        parent["source_revision"] = "other"
    elif change == "duplicate":
        table["methods"][0]["cells"].append(deepcopy(cell))
    else:
        cell["checked"] = [p for p in cell["checked"] if not p["path"].endswith("preflight.json")]
    with pytest.raises(ValueError):
        runner.control_readback(tmp_path, table, parent)


def test_real_small_guarded_full_pipeline_validates_actual_policy_and_routes(tmp_path):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from benchmark_tools.run_full_search_hillclimb import validate_native
    source = Path(__file__).resolve().parents[2]
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools): pytest.skip("Requires MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    folder, routes = tmp_path / "native", tmp_path / "routes"
    folder.mkdir()
    argv = runner.command(source, fasta, folder, tools, routes, cpu=2)
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "PYTHONPATH": str(source)}
    child = subprocess.run(argv, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr + (folder / "output/benchmark.log").read_text()
    phylo = folder / "output/orthohmm_phylogeny"
    case = {"fasta": str(fasta), "output": str(folder), "inputs": [record(p) for p in sorted(fasta.iterdir())],
        "search_candidate_policy": "protected_union", "search_min_total_hits": 4, "search_max_candidates_per_query": 100,
        "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}}
    harness = json.loads((folder / "metrics.json").read_text())["harness"]
    plan = {"case": case, "cpu": 2, "source_root": str(source), "source_revision": harness["git_commit"],
        "control": {"coverage": {"genes": 16, "species": 4}}, "native_routes": str(routes)}
    native, route, pins = runner.validate_variant(plan)
    assert (native["genes"], native["species"]) == (16, 4)
    assert route["totals"]["initial_search"]["scored_pairs"] == native["metrics"]["counts"]["search_candidates"]
    assert pins and validate_native(case, 2)["status"] == "fresh_native_content_verified"
    plan["source_revision"] = "different"
    with pytest.raises(ValueError, match="scientific source"):
        runner.validate_variant(plan)
    plan["source_revision"] = harness["git_commit"]
    plan["control"]["coverage"]["genes"] = 17
    with pytest.raises(ValueError, match="Actual intervention"):
        runner.validate_variant(plan)
