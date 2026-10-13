import json
from pathlib import Path

import pytest

from benchmark_tools import report_copy_branch_cells as reporting
from benchmark_tools import report_hillclimb as combined
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return record(path)


def candidate(root, content):
    root.mkdir(parents=True)
    seed = root / "seed.txt"
    seed.write_text(content)
    groups = root / "groups.txt"
    groups.write_text(content)
    constraints = root / "constraints.json"
    constraints.write_text(json.dumps({"groups": content}))
    return {"seed_partition": record(seed), "candidate_partition": record(groups),
            "membership_constraints": record(constraints), "expansion": {"parameters": {"fixed": True}}}


def fixture(tmp_path, monkeypatch, dataset="OrthoBench"):
    root = tmp_path / "branch"
    root.mkdir()
    old = tmp_path / "parent"
    cases, parent_rows, old_candidates = [], [], []
    score = {"OrthoBench F1": .74} if dataset == "OrthoBench" else {"Three Kingdoms F1": .87} if dataset == "ThreeKingdoms" else {
        "VGNC F": .9, "SwissTrees F": .83, "TreeFam-A F": .61, "EC": .96, "GO": .49, "FAS": .76}
    products = {key: {"path": str(old / key), "sha256": key} for key in PRODUCTS}
    for index, label in enumerate(reporting.runner.PARENT_LABELS):
        cand = candidate(old / label, "a\nb\n" if index == 0 else "a b\n")
        old_candidates.append({"candidate": cand["candidate_partition"], "constraints": cand["membership_constraints"]})
        cases.append({"label": label, "output": str(old / label), "command": ["python", "replay.py"]})
        parent_rows.append({"label": label, "status": "independent_content_readback_verified", "scores": dict(score)})
    old_plan = {"source_root": str(tmp_path / "old_source"), "source_revision": "old", "cases": cases}
    old_identity = write(old / "execution_plan.json", old_plan)
    top = {"job_id": "123", "execution_plan": old_identity, "source_revision": "old", "status": "finished_all_arms_verified"}
    child = {"status": "finished_all_arms_verified", "arms": parent_rows}
    if dataset == "QfO":
        top.update(child)
        child = top
        child_path = old / "results.json"
    else:
        child_path = old / dataset / "results.json"
        write(child_path, child)
    write(old / "results.json", top)
    scheduler = {"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0", "Elapsed": "01:00:00", "NodeList": "bizon", "AllocCPUS": "16"}
    resolution = {"status": "terminal_matched_parent_independently_admitted", "scheduler": scheduler,
        "checked": [old_identity, record(old / "results.json"), record(child_path)],
        "captures": [], "parent_rows": parent_rows}
    receipt = write(root / "resolved_parent_inputs.json", resolution)
    equivalence = {"status": "unchanged"}
    arms = [{"label": label, "profile": profile, "output": str(root / label), "scoring_id": f"short{index}",
        "command": ["python", "replay.py"]} for index, (label, profile) in enumerate(reporting.runner.ARMS)]
    plan = {"repo": str(tmp_path), "output": str(root), "dataset": dataset, "parent_root": str(old), "parent_job": "123",
        "parent_results": str(child_path), "parent_plan": old_identity, "parent_cases": cases, "parent_candidates": old_candidates,
        "source_root": str(tmp_path / "source"), "source_revision": "new", "replay_equivalence": equivalence,
        "arms": arms, "data": {}}
    rows = {}
    for index, label in enumerate(reporting.LABELS):
        ref = 0 if label in ("control", "large_only500") else 1
        cand = candidate(root / label, "a\nb\n" if ref == 0 else "a b\n")
        trace = write(root / label / "trace.json", {"events": []})
        rows[label] = {"label": label, "profile": arms[index - 2]["profile"] if index >= 2 else "default",
            "status": "independent_content_readback_verified", "scores": dict(score), "candidate": cand,
            "seed": cand["seed_partition"], "copy_split_trace": trace, "copy_split_counts": {},
            "inputs_reverified_unchanged": True, "execution_kind": "exact_candidate_structural_reuse_no_native_rerun",
            "reused_parent_label": reporting.runner.PARENT_LABELS[ref], "native_case": cases[ref],
            "native_products": products, "native_summary": {"species_tree_families": 2}, "scoring_output": cases[ref]["output"]}
    producer = {"dataset": dataset, "source_revision": "new", "execution_plan": write(root / "execution_plan.json", plan),
        "job_id": "456", "status": "finished_all_branch_arms_verified", "parent_admission": receipt,
        "arms": [rows[label] for label in reporting.LABELS[2:]], "structural_controls": {label: rows[label] for label in reporting.LABELS[:2]}}
    write(root / "results.json", producer)
    native = {"genes": combined.EXPECTED_GENES[dataset], "species": 78 if dataset == "QfO" else 12,
        "groups": 2, "pairs": 2, "products": products, "checked": [], "summary": {"species_tree_families": 2}}
    monkeypatch.setattr(reporting.runner, "verify", lambda *args: None)
    monkeypatch.setattr(reporting.runner, "cases_for", lambda *args: cases)
    monkeypatch.setattr(reporting.runner, "replay_equivalence", lambda *args: equivalence)
    monkeypatch.setattr(reporting.subprocess, "check_output", lambda *args, **kwargs: "|".join(scheduler) + "\n" + "|".join(scheduler.values()) + "\n")
    monkeypatch.setattr(reporting, "audit_deferred_cell", lambda repo, kind, dataset, row, *args: {
        "status": "independent_content_readback_verified", "scores": combined.endpoint_vector(dataset, row["scores"]),
        "products": products, "checked": []})
    monkeypatch.setattr(reporting, "replay_native", lambda *args: native)
    monkeypatch.setattr(reporting, "readback_scores", lambda repo, dataset, row, *args: (row["scores"], {}, []))
    return root, plan, producer, rows, native, resolution


def test_missing_branch_reports_do_not_invent_scores(tmp_path):
    cells, captures = reporting.collect_branch_cells(tmp_path, "OrthoBench", tmp_path, tmp_path)
    assert len(cells) == 4 and captures == []
    assert all(c["scores"] == {} and "missing" in c["status"] for c in cells)


@pytest.mark.parametrize("dataset", ["OrthoBench", "ThreeKingdoms"])
def test_exact_narrow_reuse_is_independently_scored_without_new_execution(tmp_path, monkeypatch, dataset):
    root, plan, producer, rows, _, _ = fixture(tmp_path, monkeypatch, dataset)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    before = {p: p.read_bytes() for p in root.rglob("*.json")}
    cells, captures = reporting.collect_branch_cells(tmp_path, dataset, root, snapshots)
    assert len(cells) == 4 and len(captures) == 2
    assert all(c["status"] == "independent_content_readback_verified" for c in cells)
    assert cells[2]["execution_kind"] == "exact_candidate_structural_reuse_no_native_rerun"
    assert cells[0]["scores"] == cells[2]["scores"]
    assert cells[1]["scores"] == cells[3]["scores"]
    assert all(p.read_bytes() == data for p, data in before.items())


@pytest.mark.parametrize("mutation", ["candidate", "constraints", "native_case", "scoring_output", "kind", "timing", "label", "profile", "products", "summary", "scores", "postflight"])
def test_reuse_rejects_changed_inputs_handoffs_and_claims(tmp_path, monkeypatch, mutation):
    _, plan, producer, rows, native, _ = fixture(tmp_path, monkeypatch)
    parent = reporting.admit_parent(tmp_path, plan, producer)
    row = rows["large_only500"]
    if mutation in ("candidate", "constraints"):
        key = "candidate_partition" if mutation == "candidate" else "membership_constraints"
        path = Path(row["candidate"][key]["path"])
        path.write_text("changed")
        row["candidate"][key] = record(path)
    elif mutation == "native_case":
        row["native_case"] = plan["parent_cases"][1]
    elif mutation == "scoring_output":
        row["scoring_output"] = "foreign"
    elif mutation == "kind":
        row["execution_kind"] = "new_cached_graph_branch_inference_own_tree"
    elif mutation == "timing":
        row["native_execution"] = {"returncode": 0}
    elif mutation == "label":
        row["reused_parent_label"] = "deferred500"
    elif mutation == "profile":
        row["profile"] = "phylogeny_deferred_medium"
    elif mutation == "products":
        row["native_products"] = {}
    elif mutation == "summary":
        row["native_summary"] = {}
    elif mutation == "scores":
        row["scores"] = {"OrthoBench F1": .99}
    else:
        row["inputs_reverified_unchanged"] = False
    with pytest.raises(ValueError):
        reporting.audit_branch_cell(tmp_path, row, plan, producer, parent)


@pytest.mark.parametrize("mutation", ["coverage", "species", "markers", "bool_markers"])
def test_invalid_native_universe_or_markers_are_not_admitted(tmp_path, monkeypatch, mutation):
    _, plan, producer, rows, native, _ = fixture(tmp_path, monkeypatch)
    parent = reporting.admit_parent(tmp_path, plan, producer)
    if mutation == "coverage":
        native["genes"] -= 1
    elif mutation == "species":
        native["species"] -= 1
    else:
        native["summary"]["species_tree_families"] = True if mutation == "bool_markers" else 0
    with pytest.raises(ValueError):
        reporting.audit_branch_cell(tmp_path, rows["large_only500"], plan, producer, parent)


def traced_row(tmp_path, label="medium_only500"):
    event = {"members": list(range(70)), "large_copy_eligible": False, "medium_copy_eligible": True,
        "deferred": True, "size": 70, "maximum_species_copies": 14, "dataset_species_count": 12}
    row = {"label": label, "copy_split_counts": {"medium_eligible_families": 1, "medium_deferred_families": 1, "medium_deferred_genes": 70}}
    row["copy_split_trace"] = write(tmp_path / "trace.json", {"events": [event]})
    return row, event


@pytest.mark.parametrize("mutation", ["wrong_scope", "size", "members", "bool_index", "counts", "dataset", "eligibility", "copies", "ceiling"])
def test_executed_trace_counts_and_bounds_are_verified(tmp_path, mutation):
    row, event = traced_row(tmp_path)
    assert reporting.trace_counts(row, "OrthoBench")["medium_deferred_genes"] == 70
    if mutation == "wrong_scope":
        row["label"] = "large_only500"
    elif mutation == "size":
        event["size"] = 69
    elif mutation == "members":
        event["members"] = [0] * 70
    elif mutation == "bool_index":
        event["members"][0] = False
    elif mutation == "counts":
        row["copy_split_counts"]["medium_deferred_genes"] = 69
    elif mutation == "dataset":
        event["dataset_species_count"] = 78
    elif mutation == "eligibility":
        event["large_copy_eligible"] = True
    elif mutation == "copies":
        event["maximum_species_copies"] = True
    else:
        event["members"] = list(range(501))
        event["size"] = 501
    row["copy_split_trace"] = write(tmp_path / "trace.json", {"events": [event]})
    with pytest.raises(ValueError):
        reporting.trace_counts(row, "OrthoBench")


def test_running_and_failed_arms_remain_scoreless(tmp_path, monkeypatch):
    root, _, producer, _, _, _ = fixture(tmp_path, monkeypatch)
    producer["status"] = "failed_with_explicit_gaps"
    producer["arms"] = [{"label": "large_only500", "status": "failed", "error": "native failure", "scores": {"OrthoBench F1": .99}}]
    write(root / "results.json", producer)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    cells, _ = reporting.collect_branch_cells(tmp_path, "OrthoBench", root, snapshots)
    assert cells[2]["scores"] == {} and cells[2]["error"] == "native failure"
    assert cells[3]["scores"] == {} and cells[3]["status"] == "not_run_terminal_producer_gap"


@pytest.mark.parametrize("mutation", ["checksum", "revision", "dataset", "output", "inventory", "parent"])
def test_plan_or_parent_failure_prevents_branch_score_admission(tmp_path, monkeypatch, mutation):
    root, plan, producer, _, _, _ = fixture(tmp_path, monkeypatch)
    if mutation == "parent":
        monkeypatch.setattr(reporting, "audit_deferred_cell", lambda *args: (_ for _ in ()).throw(ValueError("control mismatch")))
    elif mutation == "checksum":
        Path(producer["execution_plan"]["path"]).write_text("{}")
    else:
        if mutation == "revision":
            plan["source_revision"] = "changed"
        elif mutation == "dataset":
            plan["dataset"] = "QfO"
        elif mutation == "output":
            plan["output"] = "foreign"
        else:
            plan["arms"].reverse()
        producer["execution_plan"] = write(root / "execution_plan.json", plan)
        write(root / "results.json", producer)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    cells, _ = reporting.collect_branch_cells(tmp_path, "OrthoBench", root, snapshots)
    assert all(c["status"] == "readback_failed" and c["scores"] == {} for c in cells)


@pytest.mark.parametrize("label", ["large_only500", "foreign"])
def test_duplicate_or_unknown_arm_rejected(tmp_path, monkeypatch, label):
    root, _, producer, _, _, _ = fixture(tmp_path, monkeypatch)
    producer["arms"].append({"label": label, "status": "native_running"})
    write(root / "results.json", producer)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    with pytest.raises(ValueError, match="Duplicate or unknown"):
        reporting.collect_branch_cells(tmp_path, "OrthoBench", root, snapshots)


def test_reporting_integrates_all_four_branch_vectors_without_mean_ranking(tmp_path, monkeypatch):
    registry = tmp_path / "registry.json"
    write(registry, {"rows": [], "limitations": []})
    def collect(repo, dataset, root, snapshots):
        assert root == tmp_path / "branches" / dataset
        return [{"method_key": "copy_branch:" + label, "dataset": dataset,
            "status": "independent_content_readback_verified", "products": {},
            "scores": combined.endpoint_vector(dataset, {metric: .8 for metric in combined.METRIC_NAMES[dataset].values()})}
            for label in reporting.LABELS], []
    monkeypatch.setattr(reporting, "collect_branch_cells", collect)
    report = combined.build_report(tmp_path, registry, tmp_path / "membership", tmp_path / "supplement",
        tmp_path / "full", tmp_path / "report", copy_branch_root=tmp_path / "branches")
    assert not report["readback_failures"]
    branch = [m for m in report["methods"] if m["key"].startswith("copy_branch:")]
    assert len(branch) == 4 and all(len(m["scores"]) == 8 for m in branch)
    assert report["round_screening"]["copy_branch"]["control"] == "control"
    assert not report["round_screening"]["copy_branch"]["comparisons"]["large_only500"]["passes_point_estimate_screen"]


def qfo_fixture(tmp_path, monkeypatch):
    from benchmark_tools.run_copy_branch_ablation import native_command
    root, plan, producer, rows, native, resolution = fixture(tmp_path, monkeypatch, "QfO")
    parent = reporting.admit_parent(tmp_path, plan, producer)
    fixed = write(tmp_path / "fixed.json", {"cases": [{"output": "unused"}, {"output": str(tmp_path / "cache")}]})
    parent["parent"]["parent_plan"] = fixed
    command = ["python", "old/replay.py", "--fasta-directory", "fasta", "--candidate-clusters", "old/candidate",
        "--membership-constraints", "old/constraints", "--output-directory", "old/output", "--json", "old/metrics.json",
        "--cpu", "16", "--species-tree-mode", "infer", "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
        "--species-tree-rooting", "min_variance", "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", "cache"]
    plan["parent_cases"][0]["command"] = command
    for arm in plan["arms"]:
        folder = Path(arm["output"])
        cand = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
        arm["command"] = native_command(command, Path(plan["source_root"]), folder, cand,
            cand.with_name("phylogeny_candidate_merges.json"), tmp_path / "cache/output")
        row = rows[arm["label"]]
        old_candidate = row["candidate"]
        new = candidate(cand.parent, Path(old_candidate["candidate_partition"]["path"]).read_text())
        Path(new["candidate_partition"]["path"]).rename(cand)
        constraint = cand.with_name("phylogeny_candidate_merges.json")
        Path(new["membership_constraints"]["path"]).rename(constraint)
        new["candidate_partition"], new["membership_constraints"] = record(cand), record(constraint)
        row.update(candidate=new, execution_kind="new_cached_graph_branch_inference_own_tree", native_execution={"returncode": 0},
                   native_case=arm, scoring_output=str(folder))
    return root, plan, producer, rows, native, parent


def test_qfo_new_branch_uses_exact_native_and_scoring_handoff(tmp_path, monkeypatch):
    _, plan, producer, rows, _, parent = qfo_fixture(tmp_path, monkeypatch)
    handoffs = []
    def readback(repo, dataset, row, native, folder, participant):
        handoffs.append((folder, participant))
        return row["scores"], {}, []
    monkeypatch.setattr(reporting, "readback_scores", readback)
    for row in (rows[label] for label in reporting.LABELS[2:]):
        result = reporting.audit_branch_cell(tmp_path, row, plan, producer, parent)
        assert result["status"] == "independent_content_readback_verified" and len(result["scores"]) == 6
    assert handoffs == [(Path(arm["output"]), "ohmm_hillclimb_r1_" + arm["scoring_id"]) for arm in plan["arms"]]


@pytest.mark.parametrize("mutation", ["returncode", "kind", "case", "folder", "candidate", "protocol"])
def test_qfo_mismatched_inference_or_scoring_rejected(tmp_path, monkeypatch, mutation):
    _, plan, producer, rows, _, parent = qfo_fixture(tmp_path, monkeypatch)
    row = rows["large_only500"]
    if mutation == "returncode":
        row["native_execution"]["returncode"] = 1
    elif mutation == "kind":
        row["execution_kind"] = "exact_candidate_structural_reuse_no_native_rerun"
    elif mutation == "case":
        row["native_case"] = plan["parent_cases"][0]
    elif mutation == "folder":
        row["scoring_output"] = "foreign"
    elif mutation == "candidate":
        row["candidate"]["candidate_partition"] = rows["medium_only500"]["candidate"]["candidate_partition"]
    else:
        command = plan["arms"][0]["command"]
        command[command.index("--pair-rule") + 1] = "lca"
    with pytest.raises(ValueError):
        reporting.audit_branch_cell(tmp_path, row, plan, producer, parent)


@pytest.mark.parametrize("mutation", ["scheduler", "job", "revision", "case", "rows", "equivalence", "no_control"])
def test_parent_admission_rejects_changed_terminal_identity(tmp_path, monkeypatch, mutation):
    _, plan, producer, _, _, resolution = fixture(tmp_path, monkeypatch)
    old = Path(plan["parent_root"])
    if mutation == "scheduler":
        monkeypatch.setattr(reporting.subprocess, "check_output", lambda *args, **kwargs:
            "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|RUNNING|0:0|01:00:00|bizon|16\n")
    elif mutation == "equivalence":
        monkeypatch.setattr(reporting.runner, "replay_equivalence", lambda *args: {"status": "changed"})
    elif mutation == "case":
        plan["parent_cases"] = []
    else:
        top = json.loads((old / "results.json").read_text())
        child = json.loads(Path(plan["parent_results"]).read_text())
        if mutation == "job":
            top["job_id"] = "other"
        elif mutation == "revision":
            top["source_revision"] = "other"
        elif mutation == "rows":
            child["arms"].reverse()
        else:
            child["arms"][0]["status"] = "failed"
            resolution["parent_rows"] = child["arms"]
        write(old / "results.json", top)
        write(Path(plan["parent_results"]), child)
        resolution["checked"] = [record(old / "results.json"), record(Path(plan["parent_results"]))]
        producer["parent_admission"] = write(Path(plan["output"]) / "resolved_parent_inputs.json", resolution)
    with pytest.raises(ValueError):
        reporting.admit_parent(tmp_path, plan, producer)
