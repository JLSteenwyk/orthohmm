from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import recover_full_qfo_search as helper
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def terminal_fixture():
    identity = {"path": "parent", "sha256": "hash", "bytes": 1}
    rows = [{"label": "control100", "status": "failed", "exit_code": 0,
             "native": {"status": "fresh_native_content_verified"}, "scores": {},
             "cpu_backend_probe": {"passed": True}, "inputs_sources_runtime_reverified": True,
             "error_type": "ValueError", "error": helper.PATH_ERROR},
            {"label": "cap400", "status": "not_run_failed_control", "scores": {},
             "inputs_sources_runtime_reverified": True}]
    report = {"job_id": "24201", "source_revision": helper.REVISION, "execution_plan": identity,
              "status": "finished_with_explicit_gaps", "arms": rows}
    return report, identity


def test_specific_scoring_failure_runs_only_missing_variant():
    report, identity = terminal_fixture()
    assert helper.recovery_actions(report, identity, "24201", {"State": "FAILED"}) == ["score_existing", "run_missing"]


@pytest.mark.parametrize("field,value", [("status", "running"), ("job_id", "24190"),
    ("source_revision", "new"), ("execution_plan", {}), ("arms", [])])
def test_unrelated_or_live_parent_is_rejected(field, value):
    report, identity = terminal_fixture()
    report[field] = value
    with pytest.raises(ValueError, match="exact terminal"):
        helper.recovery_actions(report, identity, "24201", {"State": "FAILED"})


@pytest.mark.parametrize("field,value", [("error", "Out of memory"), ("error_type", "RuntimeError"),
    ("exit_code", 1), ("native", {}), ("cpu_backend_probe", {"passed": False}),
    ("inputs_sources_runtime_reverified", False), ("scores", {"VGNC F": .9})])
def test_unexpected_native_runtime_or_scoring_failure_is_not_bypassed(field, value):
    report, identity = terminal_fixture()
    report["arms"][0][field] = value
    with pytest.raises(ValueError):
        helper.recovery_actions(report, identity, "24201", {"State": "FAILED"})


def test_completed_native_variant_is_scored_not_repeated():
    report, identity = terminal_fixture()
    report["arms"][1] = {**deepcopy(report["arms"][0]), "label": "cap400"}
    assert helper.recovery_actions(report, identity, "24201", {"State": "FAILED"}) == ["score_existing"] * 2
    for row in report["arms"]:
        row["status"] = "scored_development_result"
    report["status"] = "finished_all_arms_scored"
    assert helper.recovery_actions(report, identity, "24201", {"State": "COMPLETED"}) == ["reuse_scored"] * 2


def test_skipped_variant_with_partial_output_claims_is_rejected():
    report, identity = terminal_fixture()
    report["arms"][1]["native"] = {}
    with pytest.raises(ValueError, match="Skipped variant"):
        helper.recovery_actions(report, identity, "24201", {"State": "FAILED"})


def test_scheduler_requires_exact_raw_terminal_job_and_resources():
    header = "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
    line = "24201|FAILED|1:0|01:00:00|bizon|32\n"
    assert helper.terminal_scheduler(header + line, "24201")["State"] == "FAILED"
    assert helper.terminal_scheduler(header + line.replace("FAILED|1:0", "COMPLETED|0:0"), "24201")["State"] == "COMPLETED"
    for bad in (line.replace("FAILED", "RUNNING"), line.replace("FAILED", "CANCELLED"),
                line.replace("1:0", "137:0"), line.replace("bizon", "dgx"),
                line.replace("|32", "|16"), line.replace("24201", "24190_1"), line + line):
        with pytest.raises(ValueError):
            helper.terminal_scheduler(header + bad, "24201")


def test_short_namespace_passes_real_repo_length_preflight():
    repo = Path(helper.__file__).resolve().parent.parent
    name = helper.short_id(repo, repo / "benchmarks/results/recovery_example", 0)
    helper.require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + name))
    for index in (True, -1, 2, "0", .5):
        with pytest.raises(ValueError, match="arm index"):
            helper.short_id(repo, repo, index)


def parent_fixture(tmp_path):
    source = tmp_path / "original_source"
    cases = []
    for label, cap in zip(helper.LABELS, (100, 400)):
        folder = tmp_path / "original" / label
        cases.append({"label": label, "dataset": "QfO", "output": str(folder), "fasta": "fasta",
            "inputs": [{"path": f"species{i}.fa"} for i in range(78)], "search_max_candidates_per_query": cap,
            "command": helper.native_command(source, Path("fasta"), folder, ["mafft", "fasttree"], cap, 32)})
    return {"source_revision": helper.REVISION, "cpu": 32, "source_root": str(source), "cases": cases}


def test_parent_command_inventory_preserves_native_method_and_full_inputs(tmp_path):
    parent = parent_fixture(tmp_path)
    assert helper.parent_cases(parent) == parent["cases"]
    for mutate in (lambda p: p.update(cpu=True), lambda p: p.update(source_revision="other"),
        lambda p: p["cases"][0]["command"].append("--unconstrained-membership"),
        lambda p: p["cases"][1]["inputs"].pop(), lambda p: p["cases"][1].update(fasta="subset")):
        broken = deepcopy(parent)
        mutate(broken)
        with pytest.raises(ValueError):
            helper.parent_cases(broken)


def test_native_launch_rejects_control_and_changed_variant_without_subprocess(tmp_path, monkeypatch):
    parent = parent_fixture(tmp_path)
    monkeypatch.setattr(helper.subprocess, "run", lambda *a, **k: pytest.fail("No native process may start"))
    with pytest.raises(ValueError, match="Only missing cap400"):
        helper.run_missing(parent["cases"][0], parent, {})
    bad = deepcopy(parent["cases"][1])
    bad["command"][bad["command"].index("--cpm-resolution") + 1] = "0.2"
    with pytest.raises(ValueError, match="Only missing cap400"):
        helper.run_missing(bad, parent, {})


def setup_worker(tmp_path, monkeypatch, actions):
    source = Path(helper.__file__).resolve().parent.parent
    parent_path = tmp_path / "parent.json"
    parent_path.write_text(json.dumps({"source_root": "original_source"}))
    cases, resolved = [], []
    for index, (label, action) in enumerate(zip(helper.LABELS, actions)):
        case = {"label": label, "output": str(tmp_path / "original" / label), "scoring_id": f"old{index}"}
        native = {"phylogeny": "native", "checked": []}
        resolved.append({"label": label, "action": action, "native": native, "parent_row": {
            "scores": {"GO": .5}, "metric_files": [], "conversion": {}, "command": []}})
        cases.append({"label": label, "original_case": case, "fresh_case": {**case, "fresh": True},
                      "scoring_output": str(tmp_path / label / "scoring"), "scoring_id": f"new{index}"})
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path), "source_revision": "test",
            "tool_path": "/usr/bin", "scope": "test", "parent_plan": record(parent_path), "cases": cases}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setattr(helper, "verify_inputs", lambda *a: None)
    monkeypatch.setattr(helper, "resolve", lambda *a: (resolved, [record(parent_path)]))
    monkeypatch.setattr(helper, "host_activity", lambda: "shared-host fixture")
    calls = []
    def native(case, *a):
        calls.append(("native", case["label"]))
        return {"phylogeny": "native", "checked": []}
    def scoring(repo, phylo, folder, scoring_id):
        calls.append(("score", scoring_id))
        return {"scores": {"GO": .5}, "metric_files": [], "conversion": {}, "command": []}
    def readback(repo, row, *a):
        calls.append(("readback", row["label"]))
        return {"status": "independent_content_readback_verified", "checked": []}
    monkeypatch.setattr(helper, "run_missing", native)
    monkeypatch.setattr(helper, "score_qfo", scoring)
    monkeypatch.setattr(helper, "independent_scores", readback)
    return path, calls


def test_worker_verifies_control_before_only_missing_variant(tmp_path, monkeypatch):
    path, calls = setup_worker(tmp_path, monkeypatch, ("score_existing", "run_missing"))
    result = helper.run(path)
    assert result["status"] == "finished_all_arms_independently_verified"
    assert calls == [("score", "new0"), ("readback", "control100"), ("native", "cap400"),
                     ("score", "new1"), ("readback", "cap400")]
    assert result["arms"][0]["native_case"]["output"].endswith("original/control100")
    assert "fresh_native_wall_s" not in result["arms"][0]
    assert "fresh_native_wall_s" in result["arms"][1]
    with pytest.raises(FileExistsError, match="already started"):
        helper.run(path)


def test_worker_failed_control_readback_prevents_variant_execution(tmp_path, monkeypatch):
    path, calls = setup_worker(tmp_path, monkeypatch, ("score_existing", "run_missing"))
    def fail(*a):
        raise ValueError("Independent control disagrees")
    monkeypatch.setattr(helper, "independent_scores", fail)
    result = helper.run(path)
    assert result["status"] == "failed_with_explicit_gaps"
    assert result["error"] == "Independent control disagrees"
    assert calls == [("score", "new0")]
    assert result["arms"][0]["status"] == "failed" and len(result["arms"]) == 1


@pytest.mark.parametrize("actions,expected", [
    (("score_existing", "score_existing"), [("score", "new0"), ("readback", "control100"),
                                            ("score", "new1"), ("readback", "cap400")]),
    (("reuse_scored", "reuse_scored"), [("readback", "control100"), ("readback", "cap400")])])
def test_worker_never_repeats_completed_native_or_scoring_work(tmp_path, monkeypatch, actions, expected):
    path, calls = setup_worker(tmp_path, monkeypatch, actions)
    result = helper.run(path)
    assert result["status"] == "finished_all_arms_independently_verified"
    assert calls == expected


def test_worker_live_parent_gate_failure_starts_nothing(tmp_path, monkeypatch):
    path, calls = setup_worker(tmp_path, monkeypatch, ("score_existing", "run_missing"))
    def fail(*a):
        raise ValueError("Parent still RUNNING")
    monkeypatch.setattr(helper, "resolve", fail)
    result = helper.run(path)
    assert result["status"] == "failed_with_explicit_gaps" and result["arms"] == [] and calls == []


def setup_resolution(tmp_path, monkeypatch):
    parent = parent_fixture(tmp_path)
    path = tmp_path / "parent_plan.json"
    path.write_text(json.dumps(parent))
    report, _ = terminal_fixture()
    report["execution_plan"] = record(path)
    native = {"genes": 984137, "species": 78, "products": {}, "checked": []}
    report["arms"][0]["native"]["products"] = {}
    results = tmp_path / "parent_results.json"
    results.write_text(json.dumps(report))
    for case in parent["cases"]:
        case["scoring_id"] = "x" * 200
    plan = {"parent_plan": record(path), "parent_results": str(results), "parent_job": "24201",
            "repo": "/r", "cases": [helper.planned_case(Path("/r"), tmp_path / "recovery", c, i)
                                      for i, c in enumerate(parent["cases"])]}
    monkeypatch.setattr(helper.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24201|FAILED|1:0|01:00:00|bizon|32\n")
    monkeypatch.setattr(helper, "verified_native", lambda *a: native)
    return plan, parent, results


def test_resolution_freezes_terminal_failed_attempt_and_admits_only_valid_native(tmp_path, monkeypatch):
    plan, parent, results = setup_resolution(tmp_path, monkeypatch)
    output = tmp_path / "recovery"
    output.mkdir()
    resolved, pins = helper.resolve(plan, parent, output)
    assert [r["action"] for r in resolved] == ["score_existing", "run_missing"]
    assert resolved[1]["native"] is None
    assert record(results) in pins
    assert json.loads((output / "parent_results_snapshot.json").read_text())["status"] == "finished_with_explicit_gaps"
    assert pins[0] == record(output / "resolved_parent_inputs.json")


def test_resolution_never_repeats_a_variant_with_unexpected_native_directory(tmp_path, monkeypatch):
    plan, parent, _ = setup_resolution(tmp_path, monkeypatch)
    Path(parent["cases"][1]["output"]).mkdir(parents=True)
    with pytest.raises(ValueError, match="do not repeat"):
        helper.resolve(plan, parent, tmp_path / "recovery")


def test_resolution_rejects_missing_or_reordered_cases(tmp_path, monkeypatch):
    plan, parent, _ = setup_resolution(tmp_path, monkeypatch)
    for cases in (plan["cases"][:1], list(reversed(plan["cases"]))):
        with pytest.raises(ValueError, match="both ordered"):
            helper.resolve({**plan, "cases": cases}, parent, tmp_path)


def test_resolution_requires_path_error_to_reproduce(tmp_path, monkeypatch):
    plan, parent, _ = setup_resolution(tmp_path, monkeypatch)
    plan["repo"] = "/r"
    parent["cases"][0]["scoring_id"] = "tiny"
    plan["cases"][0] = helper.planned_case(Path("/r"), tmp_path / "recovery", parent["cases"][0], 0)
    with pytest.raises(ValueError, match="does not reproduce"):
        helper.resolve(plan, parent, tmp_path / "recovery")


def test_resolution_refuses_changed_native_receipt(tmp_path, monkeypatch):
    plan, parent, _ = setup_resolution(tmp_path, monkeypatch)
    monkeypatch.setattr(helper, "verified_native", lambda *a: {"products": {"new": {}}, "checked": []})
    with pytest.raises(ValueError, match="products differ"):
        helper.resolve(plan, parent, tmp_path / "recovery")


def test_resolution_rejects_changed_native_or_scoring_namespaces(tmp_path, monkeypatch):
    plan, parent, _ = setup_resolution(tmp_path, monkeypatch)
    plan["cases"][1]["fresh_case"]["output"] = parent["cases"][1]["output"]
    with pytest.raises(ValueError, match="case inventory"):
        helper.resolve(plan, parent, tmp_path / "recovery")


def test_native_admission_rejects_wrong_universe_and_source(tmp_path, monkeypatch):
    case = parent_fixture(tmp_path)["cases"][0]
    native = {"genes": 984137, "species": 78, "metrics": {"harness": {
        "git_commit": "wrong", "exit_code": 0}}}
    monkeypatch.setattr(helper, "validate_native", lambda *a: native)
    with pytest.raises(ValueError, match="source revision"):
        helper.verified_native(case)
    native["genes"] = 984136
    with pytest.raises(ValueError, match="QfO universe"):
        helper.verified_native(case)


def test_native_admission_checks_exact_harness_sources_inputs_and_checkpoint_pins(tmp_path, monkeypatch):
    case = parent_fixture(tmp_path)["cases"][0]
    source = Path(case["command"][1]).parent.parent
    (source / "orthohmm").mkdir(parents=True)
    (source / "orthohmm/core.py").write_text("pass\n")
    (source / "benchmark_tools").mkdir()
    (source / "benchmark_tools/benchmark_production.py").write_text("pass\n")
    folder = Path(case["output"])
    paths = ["metrics.json", "output/orthohmm_phylogeny/provenance_manifest.json",
        "output/orthohmm_phylogeny/reconciliation_summary.json",
        "output/orthohmm_working_res/orthohmm_edges_clustered.txt",
        "output/orthohmm_working_res/phylogeny_candidate_superfamilies.txt",
        "output/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json"]
    for relative in paths:
        path = folder / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n")
    harness = {"git_commit": helper.REVISION, "exit_code": 0, "source_manifest": [],
               "input_manifest": [{"path": f"species{i}.fa", "bytes": 1, "sha256": "hash"} for i in range(78)]}
    case["inputs"] = harness["input_manifest"]
    for path in [source / "orthohmm/core.py", source / "benchmark_tools/benchmark_production.py"]:
        item = record(path)
        harness["source_manifest"].append({**item, "path": str(path.relative_to(source))})
    native = {"genes": 984137, "species": 78, "metrics": {"harness": harness}, "products": {}}
    monkeypatch.setattr(helper, "validate_native", lambda *a: native)
    verified = helper.verified_native(case)
    assert len(verified["checked"]) == 6
    assert verified["prediction"].endswith("orthohmm_edges_clustered.txt")
    harness["source_manifest"][0]["sha256"] = "changed"
    with pytest.raises(ValueError, match="manifests differ"):
        helper.verified_native(case)

