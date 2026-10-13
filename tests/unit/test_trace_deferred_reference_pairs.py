from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import trace_deferred_reference_pairs as wrapper
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return record(path)


def row(a="A", b="B", category="TP", stage="predicted", cross="True"):
    return {"left": a, "right": b, "reference_label": "family", "native_gene_left": "sp|" + a,
        "native_gene_right": "sp|" + b, "reference_is_cross_species": cross, "category": category,
        "stage": stage, "initial_hit_directions": "2", "forward_normalized_score": "4.2",
        "reverse_normalized_score": "4.3", "same_root_hog": "True"}


def helper_sources(tmp_path, monkeypatch):
    repo, source = tmp_path / "repo", tmp_path / "source"
    relative = Path("benchmark_tools/run_membership_hillclimb.py")
    historical = ("import os\nPRODUCTS = ('pairs',)\n"
        "def argument(argv, flag):\n    return argv[argv.index(flag) + 1]\n"
        "def score_qfo(repo, folder):\n    return folder\n")
    for root in (repo, source):
        path = root / relative
        path.parent.mkdir(parents=True)
        path.write_text(historical)
    old = record(repo / relative)
    current = source / relative
    current.write_text(historical.replace("return folder", "return repo / folder"))
    commands = []
    def git_show(command):
        commands.append(command)
        return historical.encode()
    monkeypatch.setattr(wrapper.subprocess, "check_output", git_show)
    return repo, source, old, current, historical, commands


def test_source_equivalence_records_historical_blob_and_only_unused_body_change(tmp_path, monkeypatch):
    repo, source, old, current, _, commands = helper_sources(tmp_path, monkeypatch)
    executed, proof = wrapper.tracer_source_identity(repo, source, old, "old-revision")
    assert executed == record(current) and executed["sha256"] != old["sha256"]
    assert proof["historical_source"] == old and proof["executed_source"] == executed
    assert proof["tracer_import"] == "argument"
    assert commands == [["git", "-C", str(source), "show", "old-revision:benchmark_tools/run_membership_hillclimb.py"]]


def test_only_verified_historical_source_leaves_live_checks(tmp_path, monkeypatch):
    repo, source, old, _, _, _ = helper_sources(tmp_path, monkeypatch)
    executed, proof = wrapper.tracer_source_identity(repo, source, old, "old-revision")
    data = write(tmp_path / "input.json", {"original": True})
    assert wrapper.live_checks([data, old], [proof]) == [data]
    assert wrapper.live_checks([data], []) == [data]
    with pytest.raises(ValueError, match="differs"):
        wrapper.live_checks([data, {**old, "sha256": executed["sha256"]}], [proof])
    with pytest.raises(ValueError, match="absent"):
        wrapper.live_checks([data], [proof])


@pytest.mark.parametrize("mutation", ["argument", "imports", "constant", "decorator", "default", "signature", "duplicate", "missing", "old_hash", "old_size", "revision", "other_module"])
def test_source_equivalence_rejects_changed_dependency_or_unproven_history(tmp_path, monkeypatch, mutation):
    repo, source, old, current, _, _ = helper_sources(tmp_path, monkeypatch)
    code = current.read_text()
    revision = "old-revision"
    if mutation == "argument":
        code = code.replace("argv.index(flag) + 1", "argv.index(flag)")
    elif mutation == "imports":
        code = code.replace("import os", "import sys")
    elif mutation == "constant":
        code = code.replace("('pairs',)", "('different',)")
    elif mutation == "decorator":
        code = code.replace("def score_qfo", "@decorator\ndef score_qfo")
    elif mutation == "default":
        code = code.replace("repo, folder):", "repo, folder=side_effect()):")
    elif mutation == "signature":
        code = code.replace("repo, folder):", "repo, folder, extra):")
    elif mutation == "duplicate":
        code += "\ndef score_qfo(repo, folder):\n    return folder\n"
    elif mutation == "missing":
        code = code.replace("def score_qfo", "def renamed")
    elif mutation == "old_hash":
        old["sha256"] = "wrong"
    elif mutation == "old_size":
        old["bytes"] += 1
    elif mutation == "revision":
        revision = None
    else:
        path = source / "benchmark_tools/other.py"
        path.write_text("changed")
        old["path"] = str(repo / "benchmark_tools/other.py")
    current.write_text(code)
    with pytest.raises(ValueError):
        wrapper.tracer_source_identity(repo, source, old, revision)


def test_paired_transitions_keep_rescues_losses_and_within_species_rows():
    left = [row("A", "B", "FN", "different_candidate_families"), row("C", "D"),
            row("E", "F", "FN", "reference_within_species", "False")]
    right = [row("A", "B"), row("C", "D", "FN", "within_candidate_not_predicted"), deepcopy(left[2])]
    result = wrapper.paired_transitions(left, right)
    assert result["rescued_asserted_pairs"] == result["lost_asserted_pairs"] == 1
    assert result["net_recalled_pairs"] == 0 and result["asserted_pairs"] == 3
    assert result["stages"]["reference_within_species->reference_within_species"] == 1


@pytest.mark.parametrize("mutation", ["length", "duplicate", "order", "label", "gene", "hits", "score", "eligibility", "category", "within_species"])
def test_transition_diagnostic_rejects_changed_cohort_or_search_data(mutation):
    left, right = [row()], [row()]
    if mutation == "length":
        right = []
    elif mutation == "duplicate":
        left *= 2
        right *= 2
    elif mutation == "order":
        left[0]["left"], left[0]["right"] = "B", "A"
        right = deepcopy(left)
    elif mutation == "within_species":
        left[0]["reference_is_cross_species"] = "False"
        right = deepcopy(left)
    elif mutation == "category":
        right[0]["category"] = "FP"
    else:
        key = {"label": "reference_label", "gene": "native_gene_left", "hits": "initial_hit_directions",
               "score": "forward_normalized_score", "eligibility": "reference_is_cross_species"}[mutation]
        right[0][key] = "changed"
    with pytest.raises(ValueError):
        wrapper.paired_transitions(left, right)


def setup_prepare(tmp_path, monkeypatch):
    source = Path(wrapper.__file__).resolve().parent.parent
    old_candidate, old_pairs, old_roots = (tmp_path / name for name in ("old_candidates", "old_pairs", "old_roots"))
    for path in (old_candidate, old_pairs, old_roots):
        path.write_text(path.name)
    previous = {"inputs": {"candidates": record(old_candidate), "pairs": record(old_pairs), "roots": record(old_roots)},
        "sources": [record(source / "benchmark_tools/trace_qfo_reference_stages.py")], "expected_universe": {"genes": 984137}}
    previous_path = tmp_path / "previous_plan.json"
    previous_identity = write(previous_path, previous)
    table = tmp_path / "old_table.tsv"
    table.write_text("fixture\n")
    trace_path = tmp_path / "previous_trace.json"
    write(trace_path, {"status": "native_reference_pair_stage_trace_verified", "reference_pairs": 23934,
        "within_species_reference_pairs_preserved": 2, "checked": [previous_identity], "table": record(table)})
    cases, cells = [], []
    for label in wrapper.LABELS:
        folder = tmp_path / label
        products = {}
        for name, original in ((PRODUCTS[0], old_roots), (PRODUCTS[1], old_pairs), (PRODUCTS[2], old_pairs)):
            path = folder / "output/orthohmm_phylogeny" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(original.read_bytes() if label == "control" else b"variant")
            products[name] = record(path)
        candidate = folder / "candidate.txt"
        candidate.write_bytes(old_candidate.read_bytes() if label == "control" else b"variant")
        raw = folder / "qfo_scoring/results/VGNC/fixture_raw.txt.gz"
        raw.parent.mkdir(parents=True)
        raw.write_text("raw")
        db = folder / "qfo_scoring/other/fixture.db"
        db.parent.mkdir(parents=True)
        db.write_text("database")
        cases.append({"label": label, "output": str(folder), "command": ["python", "replay.py", "--candidate-clusters", str(candidate),
            "--species-tree-mode", "infer", "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy"]})
        cells.append({"method_key": "deferred_inferred:" + label, "status": "independent_content_readback_verified",
            "dataset": "QfO", "coverage": {"genes": 984137, "species": 78, "pairs": 5}, "products": products, "checked": [],
            "details": {"VGNC F": {"precision": .99, "recall": .82}}, "scores": {"QfO: VGNC F": .9}})
    parent_path = tmp_path / "execution_plan.json"
    parent_identity = write(parent_path, {"repo": str(source), "cases": cases})
    full = tmp_path / "full.json"
    full_identity = write(full, {"cells": cells, "captures": [{"source_path": str(parent_path), "snapshot": parent_identity}]})
    receipt = tmp_path / "receipt.json"
    write(receipt, {"full_report": full_identity})
    monkeypatch.setattr(wrapper.subprocess, "check_output", lambda *args, **kwargs: "fixture\n")
    return source, receipt, previous_path, trace_path, full, cases, cells


def test_preparation_binds_ordered_admitted_arms_and_keeps_old_source_inputs_read_only(tmp_path, monkeypatch):
    source, receipt, previous, trace, full, cases, cells = setup_prepare(tmp_path, monkeypatch)
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    plan = wrapper.prepare(source, source, receipt, previous, trace, tmp_path / "new")
    assert tuple(c["label"] for c in plan["cases"]) == wrapper.LABELS
    assert plan["cpus"] == 1 and plan["memory_gib"] == 8
    assert plan["cases"][1]["trace_plan"]["expected_f1"] == .9
    assert all(p.read_bytes() == data for p, data in before.items())
    with pytest.raises(ValueError, match="new diagnostic"):
        wrapper.prepare(source, source, receipt, previous, trace, tmp_path / "new")


@pytest.mark.parametrize("mutation", ["status", "coverage", "species", "order", "products", "control_bytes", "raw", "source", "plan_order", "tree_mode"])
def test_unadmitted_changed_or_ambiguous_inputs_fail_before_preparing_output(tmp_path, monkeypatch, mutation):
    source, receipt, previous, trace, full, cases, cells = setup_prepare(tmp_path, monkeypatch)
    content = json.loads(full.read_text())
    if mutation == "status":
        content["cells"][0]["status"] = "native_running"
    elif mutation in ("coverage", "species"):
        content["cells"][1]["coverage"]["genes" if mutation == "coverage" else "species"] -= 1
    elif mutation == "order":
        content["cells"].reverse()
    elif mutation == "products":
        content["cells"][1]["products"][PRODUCTS[0]] = cells[0]["products"][PRODUCTS[0]]
    elif mutation == "control_bytes":
        path = Path(cases[0]["command"][3])
        path.write_text("changed")
    elif mutation == "raw":
        raw = tmp_path / "control/qfo_scoring/results/VGNC/second_raw.txt.gz"
        raw.write_text("second")
    elif mutation == "source":
        payload = json.loads(previous.read_text())
        payload["sources"][0]["sha256"] = "wrong"
        write(previous, payload)
        previous_trace = json.loads(trace.read_text())
        previous_trace["checked"][0] = record(previous)
        write(trace, previous_trace)
    else:
        parent_path = Path(content["captures"][0]["snapshot"]["path"])
        parent = json.loads(parent_path.read_text())
        if mutation == "plan_order":
            parent["cases"].reverse()
        else:
            command = parent["cases"][0]["command"]
            command[command.index("--species-tree-mode") + 1] = "supplied"
        content["captures"][0]["snapshot"] = write(parent_path, parent)
    write(receipt, {"full_report": write(full, content)})
    with pytest.raises(ValueError):
        wrapper.prepare(source, source, receipt, previous, trace, tmp_path / "new")
    assert not (tmp_path / "new").exists()


def setup_worker(tmp_path, monkeypatch):
    source = Path(wrapper.__file__).resolve().parent.parent
    initial = tmp_path / "control_table.tsv"
    rows = [row(str(index).zfill(6), "Z") for index in range(23932)] + [
        row("Y1", "Z", "FN", "reference_within_species", "False"), row("Y2", "Z", "FN", "reference_within_species", "False")]
    def table(path, rows):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return record(path)
    initial_identity = table(initial, rows)
    plan = {"source_root": str(source), "source_revision": "fixture", "output": str(tmp_path), "inputs": [], "sources": [],
        "original_control_table": initial_identity, "cases": [{"label": label, "trace_plan": {}} for label in wrapper.LABELS], "scope": "test"}
    path = tmp_path / "execution_plan.json"
    write(path, plan)
    calls = []
    def trace_run(child, output):
        calls.append(output.name)
        output.mkdir()
        result = {"reference_pairs": 23934, "within_species_reference_pairs_preserved": 2,
            "table": table(output / "reference_pair_stages.tsv", rows), "summary": {}, "raw_categories": {}, "resources": {}}
        write(output / "report.json", result)
        return result
    monkeypatch.setattr(wrapper.tracing, "run", trace_run)
    return path, calls, trace_run


def test_worker_runs_control_then_variant_and_verifies_physical_paired_cohort(tmp_path, monkeypatch):
    path, calls, _ = setup_worker(tmp_path, monkeypatch)
    result = wrapper.run(path)
    assert calls == ["control", "deferred500"]
    assert result["status"] == "complete_paired_reference_localization_verified"
    assert result["paired_transitions"]["asserted_pairs"] == 23934
    assert result["paired_transitions"]["categories"] == {"TP->TP": 23932, "FN->FN": 2}
    with pytest.raises(ValueError, match="untouched"):
        wrapper.run(path)


@pytest.mark.parametrize("mutation", ["control_table", "cohort", "failure"])
def test_control_failure_is_preserved_and_prevents_variant_analysis(tmp_path, monkeypatch, mutation):
    path, calls, original = setup_worker(tmp_path, monkeypatch)
    def fail(child, output):
        if mutation == "failure":
            calls.append(output.name)
            raise ValueError("broken native audit")
        result = original(child, output)
        if mutation == "cohort":
            result["reference_pairs"] -= 1
        else:
            result["table"]["sha256"] = "wrong"
        return result
    monkeypatch.setattr(wrapper.tracing, "run", fail)
    result = wrapper.run(path)
    assert calls == ["control"] and result["status"] == "failed_with_explicit_gap"
    assert result["arms"][0]["status"] == "failed" and result["error"]
