import json

import pytest

from benchmark_tools import report_deferred_cells as deferred
from benchmark_tools import report_hillclimb as reporting
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return record(path)


def fixture_reports(tmp_path, dataset="QfO", kind="deferred_fixed"):
    root = tmp_path / kind
    root.mkdir()
    row = {"label": "control", "status": "independent_content_readback_verified", "scores": {}}
    cases = [{"label": label, "command": [], "output": str(root / label)} for label in deferred.LABELS]
    plan = {"source_revision": "fixture", "output": str(root), "cases": cases}
    if dataset != "QfO":
        plan["cases"] = [{"dataset": name, "arms": cases} for name in ("OrthoBench", "ThreeKingdoms")]
    identity = write(root / "execution_plan.json", plan)
    producer = {"source_revision": "fixture", "execution_plan": identity, "job_id": "123",
                "arms": [row, {"label": "deferred500", "status": "native_running", "scores": {"VGNC F": .99}}]}
    write(root / "results.json", producer)
    if dataset != "QfO":
        write(root / dataset / "results.json", {"dataset": dataset, "arms": producer.pop("arms")})
        write(root / "results.json", producer)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    return root, snapshots, row, plan, producer


def verified_cell(dataset, score=.8):
    scores = {key: score for key in ENDPOINTS if key.startswith(dataset + ": ")}
    return {"status": "independent_content_readback_verified", "scores": scores, "products": {}}


def test_missing_copy_deferral_results_remain_missing_not_zero(tmp_path):
    cells, captures = deferred.collect_deferred_cells(tmp_path, "deferred_inferred", "QfO", tmp_path, tmp_path)
    assert [c["method_key"] for c in cells] == ["deferred_inferred:control", "deferred_inferred:deferred500"]
    assert all(c["scores"] == {} and "missing" in c["status"] for c in cells)
    assert captures == []


@pytest.mark.parametrize("dataset,kind", [("QfO", "deferred_fixed"), ("QfO", "deferred_inferred"),
                                          ("OrthoBench", "deferred_inferred"), ("ThreeKingdoms", "deferred_inferred")])
def test_only_read_back_cells_receive_scores_and_sources_stay_read_only(tmp_path, monkeypatch, dataset, kind):
    root, snapshots, _, _, _ = fixture_reports(tmp_path, dataset, kind)
    originals = {path: path.read_bytes() for path in root.rglob("*.json")}
    monkeypatch.setattr(deferred, "audit_deferred_cell", lambda *args: verified_cell(dataset))
    cells, captures = deferred.collect_deferred_cells(tmp_path, kind, dataset, root, snapshots)
    assert cells[0]["scores"] == verified_cell(dataset)["scores"]
    assert cells[1]["scores"] == {} and cells[1]["claimed_scores"] == {"VGNC F": .99}
    assert all(path.read_bytes() == data for path, data in originals.items())
    assert len(captures) == (2 if dataset == "QfO" else 3)
    assert all(item["snapshot"]["path"].startswith(str(snapshots)) for item in captures)


@pytest.mark.parametrize("mutation", ["checksum", "revision", "output", "audit"])
def test_input_or_readback_failure_is_retained_without_scores(tmp_path, monkeypatch, mutation):
    root, snapshots, _, plan, producer = fixture_reports(tmp_path)
    if mutation == "checksum":
        plan["source_revision"] = "changed"
        write(root / "execution_plan.json", plan)
    elif mutation in ("revision", "output"):
        plan["source_revision" if mutation == "revision" else "output"] = "foreign"
        producer["execution_plan"] = write(root / "execution_plan.json", plan)
        write(root / "results.json", producer)
    def audit(*args):
        raise ValueError("readback mismatch")
    monkeypatch.setattr(deferred, "audit_deferred_cell", audit)
    cells, _ = deferred.collect_deferred_cells(tmp_path, "deferred_fixed", "QfO", root, snapshots)
    assert cells[0]["status"] == "readback_failed" and cells[0]["scores"] == {}
    assert cells[0]["error"]


def test_failed_and_not_yet_reported_cells_are_not_dropped(tmp_path):
    root, snapshots, row, _, producer = fixture_reports(tmp_path)
    producer["arms"] = [{**row, "status": "failed", "scores": {"VGNC F": .99}}]
    write(root / "results.json", producer)
    cells, _ = deferred.collect_deferred_cells(tmp_path, "deferred_fixed", "QfO", root, snapshots)
    assert cells[0]["status"] == "failed" and cells[0]["scores"] == {}
    assert "not_yet_reported" in cells[1]["status"]


@pytest.mark.parametrize("label", ["control", "unknown"])
def test_duplicate_or_unknown_arm_is_rejected(tmp_path, label):
    root, snapshots, row, _, producer = fixture_reports(tmp_path)
    producer["arms"] = [row, {**row, "label": label}]
    write(root / "results.json", producer)
    with pytest.raises(ValueError, match="Duplicate or unknown"):
        deferred.collect_deferred_cells(tmp_path, "deferred_fixed", "QfO", root, snapshots)


def test_inferred_child_dataset_mismatch_is_rejected(tmp_path):
    root, snapshots, _, _, _ = fixture_reports(tmp_path, "OrthoBench", "deferred_inferred")
    write(root / "OrthoBench/results.json", {"dataset": "QfO", "arms": []})
    with pytest.raises(ValueError, match="dataset identity"):
        deferred.collect_deferred_cells(tmp_path, "deferred_inferred", "OrthoBench", root, snapshots)


def test_planned_inventory_is_strict(tmp_path):
    _, _, _, plan, _ = fixture_reports(tmp_path)
    assert deferred.planned_case(plan, "QfO", "deferred500")[0]["label"] == "deferred500"
    plan["cases"].reverse()
    with pytest.raises(ValueError, match="arm inventory"):
        deferred.planned_case(plan, "QfO", "control")


@pytest.mark.parametrize("control_status", ["failed", "native_running", "independent_content_readback_verified"])
def test_variant_requires_independent_matched_control_admission(tmp_path, monkeypatch, control_status):
    root, snapshots, row, _, producer = fixture_reports(tmp_path)
    producer["arms"] = [{**row, "status": control_status}, {**row, "label": "deferred500"}]
    write(root / "results.json", producer)
    def audit(*args):
        if args[3]["label"] == "control":
            raise ValueError("control does not reproduce")
        return verified_cell("QfO")
    monkeypatch.setattr(deferred, "audit_deferred_cell", audit)
    cells, _ = deferred.collect_deferred_cells(tmp_path, "deferred_fixed", "QfO", root, snapshots)
    assert cells[1]["status"] == "readback_failed" and cells[1]["scores"] == {}
    assert "matched-control" in cells[1]["error"]


def fixture_audit(tmp_path, monkeypatch, dataset="OrthoBench", fixed=False):
    from benchmark_tools import run_held_tree_hillclimb as held
    from benchmark_tools import run_three_kingdoms_membership_supplement as supplement
    kind = "deferred_fixed" if fixed else "deferred_inferred"
    root = tmp_path / "output"
    root.mkdir()
    candidate = {}
    for key in ("candidate_partition", "membership_constraints", "seed_partition"):
        path = tmp_path / key
        path.write_text(key)
        candidate[key] = record(path)
    tree = tmp_path / "baseline.nwk"
    tree.write_text("(A:1,B:1);\n")
    command = ["python", "replay.py", "--candidate-clusters", candidate["candidate_partition"]["path"],
        "--membership-constraints", candidate["membership_constraints"]["path"],
        "--species-tree-mode", "supplied" if fixed else "infer", "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair"]
    if fixed:
        command.extend(["--species-tree", str(tree)])
    cases = [{"label": label, "command": list(command), "output": str(root / label), "scoring_id": label}
             for label in deferred.LABELS]
    scores = {key: .8 for key in reporting.METRIC_NAMES[dataset].values()}
    settings = {"golden": {}, "expected_control_scores": dict(scores) if dataset == "QfO" else reporting.endpoint_vector(dataset, scores),
                "baseline_candidate": candidate["candidate_partition"], "baseline_constraints": candidate["membership_constraints"]}
    plan = {"output": str(root), "inputs": [], "sources": [], "runtime": {"files": []}, "tree": record(tree)}
    if dataset == "QfO":
        plan.update(settings, cases=cases, evaluation_mode="fixed_species_tree" if fixed else "inferred_parent")
    else:
        plan["cases"] = [{**settings, "dataset": name, "arms": cases} for name in ("OrthoBench", "ThreeKingdoms")]
    producer = {"execution_plan": write(tmp_path / "plan.json", plan)}
    native = {"genes": reporting.EXPECTED_GENES[dataset], "species": 78 if dataset == "QfO" else 12,
              "groups": 1, "pairs": 1, "phylogeny": str(root), "products": {}, "checked": [],
              "summary": {"reconciled_families": 2, "checkpoint_hits": 2, "species_tree_families": 2}}
    row = {"label": "control", "status": "independent_content_readback_verified",
           "inputs_reverified_unchanged": True, "candidate": candidate}
    monkeypatch.setattr(deferred, "replay_native", lambda *args: native)
    monkeypatch.setattr(deferred, "readback_scores", lambda *args: (scores, {"fixture": True}, []))
    monkeypatch.setattr(deferred, "compare_golden", lambda *args: None)
    monkeypatch.setattr(supplement, "candidate_control_matches", lambda *args: None)
    monkeypatch.setattr(held, "tree_contrast", lambda *args: {"rooted_clade_distance": 0})
    return kind, row, plan, producer, cases, native


def resolve_fixture_parent(tmp_path, plan, producer):
    plan["parent_job"] = "123"
    plan["parent_plan"] = write(tmp_path / "parent_plan.json", {"fixture": True})
    parent = {"status": "finished_all_arms_verified", "job_id": "123", "execution_plan": plan["parent_plan"],
              "arms": [{"label": label, "status": "independent_content_readback_verified",
                        "inputs_reverified_unchanged": True} for label in deferred.LABELS]}
    snapshot = write(tmp_path / "parent_snapshot.json", parent)
    cache = tmp_path / "tree_cache.txt"
    cache.write_text("immutable raw tree")
    resolution = {"status": "terminal_parent_native_inputs_resolved_and_verified", "snapshot": snapshot,
                  "checked": [snapshot, record(cache)]}
    receipt = write(tmp_path / "output/resolved_parent_inputs.json", resolution)
    producer["parent_admission"] = receipt
    return receipt, resolution, parent, cache


def test_own_inferred_qfo_requires_and_admits_resolved_parent_inputs(tmp_path, monkeypatch):
    kind, row, plan, producer, _, _ = fixture_audit(tmp_path, monkeypatch, "QfO")
    receipt, _, _, _ = resolve_fixture_parent(tmp_path, plan, producer)
    result = deferred.audit_deferred_cell(tmp_path, kind, "QfO", row, plan, producer)
    assert result["scores"] == verified_cell("QfO")["scores"]
    assert receipt in result["checked"]
    assert result["tree_mode"] == "candidate_specific_inferred"


@pytest.mark.parametrize("mutation", ["receipt", "resolution", "snapshot", "parent_running", "parent_job",
                                      "cache", "raw_reuse", "markers", "mode", "fas", "score"])
def test_own_qfo_resolution_and_control_gates(tmp_path, monkeypatch, mutation):
    kind, row, plan, producer, _, native = fixture_audit(tmp_path, monkeypatch, "QfO")
    receipt, resolution, parent, cache = resolve_fixture_parent(tmp_path, plan, producer)
    if mutation == "receipt":
        producer["parent_admission"] = write(tmp_path / "foreign.json", resolution)
    elif mutation == "resolution":
        resolution["status"] = "unresolved"
        producer["parent_admission"] = write(tmp_path / "output/resolved_parent_inputs.json", resolution)
    elif mutation == "snapshot":
        write(tmp_path / "parent_snapshot.json", {})
    elif mutation in ("parent_running", "parent_job"):
        parent["status" if mutation == "parent_running" else "job_id"] = "foreign"
        resolution["snapshot"] = write(tmp_path / "parent_snapshot.json", parent)
        producer["parent_admission"] = write(tmp_path / "output/resolved_parent_inputs.json", resolution)
    elif mutation == "cache":
        cache.write_text("changed")
    elif mutation == "raw_reuse":
        native["summary"]["checkpoint_hits"] = 1
    elif mutation == "markers":
        native["summary"]["species_tree_families"] = False
    elif mutation == "mode":
        plan["evaluation_mode"] = "fixed_species_tree"
    else:
        plan["expected_control_scores"]["FAS" if mutation == "fas" else "VGNC F"] = .9
    if mutation == "fas":
        result = deferred.audit_deferred_cell(tmp_path, kind, "QfO", row, plan, producer)
        assert result["scores"]["QfO: FAS"] == .8
    else:
        with pytest.raises(ValueError):
            deferred.audit_deferred_cell(tmp_path, kind, "QfO", row, plan, producer)


@pytest.mark.parametrize("dataset,fixed", [("OrthoBench", False), ("ThreeKingdoms", False), ("QfO", True)])
def test_audit_admits_matched_native_controls(tmp_path, monkeypatch, dataset, fixed):
    kind, row, plan, producer, _, _ = fixture_audit(tmp_path, monkeypatch, dataset, fixed)
    result = deferred.audit_deferred_cell(tmp_path, kind, dataset, row, plan, producer)
    assert result["scores"] == verified_cell(dataset)["scores"]
    assert result["tree_mode"] == ("supplied_fixed_diagnostic" if fixed else "candidate_specific_inferred")


@pytest.mark.parametrize("mutation", ["postflight", "mode", "supplied", "root", "pairs", "rooting", "policy",
                                      "unconstrained", "candidate", "constraints", "markers", "species", "genes", "score"])
def test_own_inferred_admission_rejects_protocol_or_control_corruption(tmp_path, monkeypatch, mutation):
    kind, row, plan, producer, cases, native = fixture_audit(tmp_path, monkeypatch)
    argv = cases[0]["command"]
    flags = {"mode": "--species-tree-mode", "root": "--root-rule", "pairs": "--pair-rule",
             "rooting": "--species-tree-rooting", "policy": "--satellite-membership-rule",
             "candidate": "--candidate-clusters", "constraints": "--membership-constraints"}
    if mutation in flags:
        argv[argv.index(flags[mutation]) + 1] = "foreign"
    elif mutation == "postflight":
        row["inputs_reverified_unchanged"] = False
    elif mutation == "supplied":
        argv.extend(["--species-tree", "foreign"])
    elif mutation == "unconstrained":
        argv.append("--unconstrained-membership")
    elif mutation == "markers":
        native["summary"]["species_tree_families"] = 0
    elif mutation in ("species", "genes"):
        native[mutation] -= 1
    else:
        plan["cases"][0]["expected_control_scores"] = {ENDPOINTS[0]: .9}
    with pytest.raises(ValueError):
        deferred.audit_deferred_cell(tmp_path, kind, "OrthoBench", row, plan, producer)


def test_fixed_tree_diagnostics_cannot_fill_own_inferred_metric_vector():
    cells = []
    for dataset, kind in (("QfO", "deferred_fixed"), ("OrthoBench", "deferred_inferred"),
                          ("ThreeKingdoms", "deferred_inferred")):
        cells.append({**verified_cell(dataset), "method_key": kind + ":control", "dataset": dataset})
    methods = reporting.combine_cells({}, cells)
    assert len(methods["deferred_fixed:control"]["scores"]) == 6
    assert len(methods["deferred_inferred:control"]["scores"]) == 2
    vectors = {key: method["scores"] for key, method in methods.items()}
    result = reporting.screen(vectors, "deferred_inferred:control")
    assert all(not row["complete_eight_endpoint_contrast"] for row in result["comparisons"].values())


def test_inferred_cells_form_complete_vector_without_fixed_tree_values():
    cells = [{**verified_cell(dataset), "method_key": "deferred_inferred:control", "dataset": dataset}
             for dataset in ("QfO", "OrthoBench", "ThreeKingdoms")]
    cells.append({**verified_cell("QfO", .99), "method_key": "deferred_fixed:control", "dataset": "QfO"})
    methods = reporting.combine_cells({}, cells)
    assert methods["deferred_inferred:control"]["scores"] == dict.fromkeys(ENDPOINTS, .8)


def test_comparison_report_optional_deferred_rounds_are_separate(tmp_path, monkeypatch):
    registry = tmp_path / "registry.json"
    write(registry, {"rows": [], "limitations": []})
    def collect(repo, kind, dataset, root, snapshots):
        return [{**verified_cell(dataset), "method_key": kind + ":" + label, "dataset": dataset}
                for label in deferred.LABELS], []
    monkeypatch.setattr(deferred, "collect_deferred_cells", collect)
    output = tmp_path / "comparison"
    report = reporting.build_report(tmp_path, registry, tmp_path / "membership", tmp_path / "supplement.json",
        tmp_path / "full", output, tmp_path / "fixed", tmp_path / "own_supplement", tmp_path / "own_qfo")
    methods = {m["key"]: m for m in report["methods"]}
    assert len(methods["deferred_fixed:control"]["scores"]) == 6
    assert len(methods["deferred_inferred:deferred500"]["scores"]) == 8
    assert not report["round_screening"]["deferred_fixed"]["comparisons"]["deferred500"]["complete_eight_endpoint_contrast"]
    assert report["round_screening"]["deferred_inferred"]["comparisons"]["deferred500"]["complete_eight_endpoint_contrast"]
    assert "deferred_fixed:control" in (output / "scores.tsv").read_text()
    assert "deferred_inferred:deferred500" in (output / "scores.tsv").read_text()
    assert report["defaults_promoted"] is False
