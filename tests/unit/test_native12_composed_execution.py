"""Prospective gates and real small tree inventories; no production authorization."""

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import native12_composed_execution as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.snapshot_runtime_trees import inventory


def store(path, value):
    path.write_text(json.dumps(value), encoding="ascii")
    return record(path)


@pytest.fixture
def history():
    refs = [dict(path=f"/fixture/review_{index}.json", sha256=str(index), bytes=1)
            for index in range(11)]
    request = dict(index=11, job_id=23985, history=refs,
                   plan=dict(path="/fixture/plan.json"), amendment=dict(path="/fixture/amendment.json"))
    runs = [dict(index=i) for i in range(13)]
    runs[11].update(cell="p1_c1_r0")
    runs[12].update(cell="p1_c1_r1", dataset="qfo_corrected")
    review = dict(schema=module.binding.composed.SCHEMA, status="native_success", index=11,
        job_id=23985, review_job_id=24034, plan=request["plan"], amendment=request["amendment"],
        terminal_reviewed=True, composed_full_review_complete=True, native_outputs_validated=True,
        primary_resources_replayed=True, shared_host_resources_reviewed=True,
        original_ordinary_full_review_success=False, current_original_inventory_equality=False,
        next_identity_authorized=False, automatic_retry=False, publication_ready=False,
        historical_failures_retained=[23986, 24033])
    return request, dict(runs=runs), review, [dict(review=ref) for ref in refs]


def test_history_gate_requires_complete_new_type_without_rewriting_authorization(history):
    module.history_gate(*history)
    assert history[2]["next_identity_authorized"] is False
    assert history[2]["original_ordinary_full_review_success"] is False


@pytest.mark.parametrize("key,value", [
    ("schema", "allocated_native_factorial_terminal_review_v1"),
    ("status", "native_failure_retained"), ("job_id", 23986), ("index", 12),
    ("review_job_id", 24033), ("plan", {}), ("amendment", {}),
    ("terminal_reviewed", False), ("composed_full_review_complete", False),
    ("native_outputs_validated", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False), ("original_ordinary_full_review_success", True),
    ("current_original_inventory_equality", True), ("next_identity_authorized", True),
    ("automatic_retry", True), ("publication_ready", True), ("historical_failures_retained", []),
])
def test_history_gate_rejects_incomplete_or_fabricated_review(history, key, value):
    history[2][key] = value
    with pytest.raises(ValueError):
        module.history_gate(*history)


@pytest.mark.parametrize("mutation", ["request_index", "request_job", "short_history",
    "short_scheduler", "reordered", "missing_final", "other_cell", "other_dataset", "bool_index"])
def test_history_gate_rejects_reordered_or_changed_recipe(history, mutation):
    request, plan, _, scheduler = history
    if mutation == "request_index":
        request["index"] = 10
    elif mutation == "request_job":
        request["job_id"] = 24034
    elif mutation == "short_history":
        request["history"].pop()
    elif mutation == "short_scheduler":
        scheduler.pop()
    elif mutation == "reordered":
        scheduler.reverse()
    elif mutation == "missing_final":
        plan["runs"].pop()
    elif mutation == "other_cell":
        plan["runs"][12]["cell"] = "p0_c1_r1"
    elif mutation == "other_dataset":
        plan["runs"][12]["dataset"] = "orthobench"
    else:
        plan["runs"][12]["index"] = True
    with pytest.raises(ValueError):
        module.history_gate(*history)


def test_history_binding_uses_full_consumer_and_original_prefix_kernel(monkeypatch, history):
    request, plan, review, scheduler = history
    execution, terminal = dict(original=True), dict(success=True)
    request_ref, review_ref = dict(request=True), dict(composed=True)
    context = (request, execution, plan, plan["runs"][11], review, {}, "group", terminal, [], {})
    calls = []
    monkeypatch.setattr(module.binding, "native_binding",
        lambda a, b: calls.append(("consumer", a, b)) or context)
    monkeypatch.setattr(module, "reviewed_history",
        lambda a, b, c, d: calls.append(("prefix", a, b, c, d)) or scheduler)
    monkeypatch.setattr(module, "check", lambda ref: None)
    returned, result = module.history_binding(request_ref, review_ref)
    assert returned is context
    assert calls[0] == ("consumer", request_ref, review_ref)
    assert calls[1] == ("prefix", request, request["amendment"], execution, plan)
    assert result["historical_prefix"] == request["history"]
    assert result["native11_scheduler"] == terminal
    assert result["next_unrun_index"] == 12
    assert result["next_identity_authorized"] is False
    assert result["original_review_translated"] is False
    assert result["production_execution_launched"] is False


@pytest.fixture
def unrun(tmp_path):
    return (dict(index=12, cell="p1_c1_r1", output_root=str(tmp_path / "run_12")),
            dict(panel_root=str(tmp_path / "panel")), tmp_path / "request_12.json")


def test_unrun_paths_are_explicit_and_absent(unrun):
    assert module.require_unrun(*unrun) == [Path(unrun[0]["output_root"]),
        Path(unrun[1]["panel_root"]) / "sessions/run_12", unrun[2]]


@pytest.mark.parametrize("which", [0, 1, 2])
@pytest.mark.parametrize("kind", ["directory", "file", "broken_symlink"])
def test_existing_attempt_or_dangling_symlink_rejects_implicit_retry(unrun, which, kind):
    paths = module.require_unrun(*unrun)
    path = paths[which]
    path.parent.mkdir(parents=True, exist_ok=True)
    if kind == "directory":
        path.mkdir()
    elif kind == "file":
        path.write_text("retained", encoding="ascii")
    else:
        path.symlink_to(path.parent / "absent_target")
    with pytest.raises(ValueError):
        module.require_unrun(*unrun)


@pytest.fixture
def runtime():
    specs = [["/fixture/os.json", "os_digest"], ["/fixture/private.json", "private_digest"]]
    rows = []
    for number, (path, digest) in enumerate(specs):
        rows.append(dict(manifest=dict(path=path, sha256=digest, bytes=1), comparison=dict(
            missing=[], changed=[], metadata_differences={}, equal=number == 1,
            added=copy.deepcopy(module.LATE_ADDITIONS) if number == 0 else [],
            expected_records=10, observed_records=12 if number == 0 else 10)))
    return dict(schema="native11_composed_runtime_gate_v1",
        status="historical_brackets_and_current_exact_additions_revalidated",
        current_original_inventory_equality=False, current_original_entries_unchanged=True,
        current_private_inventory_equality=True, continuous_runtime_integrity_established=False,
        current_additions=copy.deepcopy(module.LATE_ADDITIONS), current_inventories=rows), specs


def test_runtime_gate_allows_only_exact_bound_late_additions(runtime):
    module.runtime_gate(*runtime)


@pytest.mark.parametrize("key,value", [
    ("schema", "runtime_tree_identity_matches"), ("status", "success"),
    ("current_original_inventory_equality", True), ("current_original_entries_unchanged", False),
    ("current_private_inventory_equality", False), ("continuous_runtime_integrity_established", True),
    ("current_additions", []), ("current_inventories", []),
])
def test_runtime_gate_rejects_wrong_scope_and_false_equality(runtime, key, value):
    runtime[0][key] = value
    with pytest.raises(ValueError):
        module.runtime_gate(*runtime)


@pytest.mark.parametrize("field,value", [
    ("missing", [dict(path="/missing")]), ("changed", [dict(path="/changed")]),
    ("metadata_differences", dict(roots=True)), ("added", []), ("equal", True),
    ("expected_records", True), ("observed_records", 13), ("observed_records", True),
])
def test_runtime_gate_rejects_any_other_os_delta(runtime, field, value):
    runtime[0]["current_inventories"][0]["comparison"][field] = value
    with pytest.raises(ValueError):
        module.runtime_gate(*runtime)


@pytest.mark.parametrize("mutation", ["private_added", "private_changed", "reorder", "digest",
    "addition_hash", "addition_mode", "addition_path", "addition_bytes"])
def test_runtime_gate_rejects_changed_private_tree_or_allowlist(runtime, mutation):
    report, specs = runtime
    if mutation == "private_added":
        report["current_inventories"][1]["comparison"]["added"] = [dict(path="/new")]
    elif mutation == "private_changed":
        report["current_inventories"][1]["comparison"]["equal"] = False
    elif mutation == "reorder":
        report["current_inventories"].reverse()
    elif mutation == "digest":
        specs[0][1] = "changed"
    else:
        report["current_additions"][0][mutation.removeprefix("addition_")] = "changed"
    with pytest.raises(ValueError):
        module.runtime_gate(report, specs)


@pytest.fixture
def real_trees(tmp_path, monkeypatch):
    roots = [tmp_path / "os", tmp_path / "private"]
    specs, original = [], []
    for number, root in enumerate(roots):
        root.mkdir()
        (root / "original").write_text("unchanged", encoding="ascii")
        manifest = inventory([root])
        original.append(manifest)
        ref = store(tmp_path / f"manifest_{number}.json", manifest)
        specs.append([ref["path"], ref["sha256"]])
    for name in ("lftp", "lftpget"):
        (roots[0] / name).write_text(name, encoding="ascii")
    additions = module.compare_inventory(original[0], inventory([roots[0]]))["added"]
    monkeypatch.setattr(module, "LATE_ADDITIONS", additions)
    source = tmp_path / "worker.py"
    source.write_text("# frozen fixture source\n", encoding="ascii")
    component = tmp_path / "component.json"
    store(component, dict(component=True))
    monkeypatch.setattr(module.binding, "composed", SimpleNamespace(__file__=str(source), RUNTIME=component))
    monkeypatch.setattr(module.binding, "COMPOSED_SOURCE_SHA", record(source)["sha256"])
    report = dict(schema="native11_composed_runtime_gate_v1",
        status="historical_brackets_and_current_exact_additions_revalidated",
        current_original_inventory_equality=False, current_original_entries_unchanged=True,
        current_private_inventory_equality=True, continuous_runtime_integrity_established=False,
        current_additions=additions, source=record(source), runtime_component=record(component),
        current_inventories=[dict(manifest=record(specs[i][0]),
            comparison=module.compare_inventory(original[i], inventory([roots[i]]))) for i in range(2)])
    runtime_ref = store(tmp_path / "runtime.json", report)
    return roots, specs, runtime_ref, report


def test_fresh_real_tree_walk_preserves_original_inequality(real_trees):
    _, specs, runtime_ref, _ = real_trees
    result = module.current_trees(specs, runtime_ref)
    assert len(result) == 2
    assert [row["original_inventory_equality"] for row in result] == [False, True]
    assert all(row["prospective_inventory_equality"] is True for row in result)
    assert all(row["schema"] == module.TREE_SCHEMA for row in result)
    assert all(row["scientific_execution_authorized"] is False for row in result)


@pytest.mark.parametrize("mutation", ["extra_os", "missing_os", "changed_os", "changed_private",
    "private_added", "late_file_changed", "manifest_changed", "worker_changed", "component_changed"])
def test_real_tree_walk_rejects_new_mutation(real_trees, mutation):
    roots, specs, runtime_ref, report = real_trees
    if mutation == "extra_os":
        (roots[0] / "extra").write_text("unexpected", encoding="ascii")
    elif mutation == "missing_os":
        (roots[0] / "original").unlink()
    elif mutation == "changed_os":
        (roots[0] / "original").write_text("changed", encoding="ascii")
    elif mutation == "changed_private":
        (roots[1] / "original").write_text("changed", encoding="ascii")
    elif mutation == "private_added":
        (roots[1] / "extra").write_text("unexpected", encoding="ascii")
    elif mutation == "late_file_changed":
        (roots[0] / "lftp").write_text("changed", encoding="ascii")
    elif mutation == "manifest_changed":
        Path(specs[0][0]).write_text("{}", encoding="ascii")
    elif mutation == "worker_changed":
        Path(report["source"]["path"]).write_text("changed", encoding="ascii")
    else:
        Path(report["runtime_component"]["path"]).write_text("changed", encoding="ascii")
    with pytest.raises(ValueError):
        module.current_trees(specs, runtime_ref)


def test_lookup_checker_hook_keeps_original_binding_and_fresh_report(tmp_path, monkeypatch, real_trees):
    _, specs, runtime_ref, _ = real_trees
    lookup_binding = store(tmp_path / "binding.json", dict(runtime_specs=specs))
    baseline = dict(path="/fixture/baseline.json", sha256="baseline")
    lookup = store(tmp_path / "lookup.json", dict(binding=lookup_binding, baseline=baseline))
    plan = dict(runtime_lookup=lookup, baseline=baseline)
    observed = {}

    class StubLookupChecker:
        def __init__(self, path, digest, output, *, tree_checker):
            observed.update(path=path, digest=digest, hook=tree_checker)
            output.mkdir()
            self.count = 0

        def __call__(self, specifications):
            self.count += 1
            return dict(runtime=observed["hook"](specifications), lookup="stub")

    monkeypatch.setattr(module, "RuntimeChecker", StubLookupChecker)
    output = tmp_path / "checks"
    checker = module.ProspectiveRuntimeChecker(plan, runtime_ref, output)
    assert checker.count == 0
    result = checker(specs)
    assert observed["path"] == Path(lookup["path"])
    assert observed["digest"] == lookup["sha256"]
    assert result["lookup"] == "stub"
    assert checker.count == 1
    receipt = json.loads((output / "tree_check_01.json").read_text())
    assert receipt["runtime_basis"] == runtime_ref
    assert receipt["inventories"] == result["runtime"]
    assert receipt["original_os_inventory_equality"] is False
    assert receipt["current_identity_revalidated"] is True
    assert receipt["next_identity_authorized"] is False
    assert receipt["check_wall_s"] > 0
    with pytest.raises(ValueError):
        checker.check_trees(list(reversed(specs)))


@pytest.mark.parametrize("drift", ["runtime", "lookup"])
def test_original_lookup_checker_and_wrapper_preserve_postflight_failure(tmp_path, real_trees, drift):
    from benchmark_tools.measure_native_scaling_run import run_checked

    roots, specs, runtime_ref, _ = real_trees
    expected = dict(python="fixture version", executable="/fixture/python", cwd="/fixture/core",
        requested=[], paths=[], modules={}, mapped_files=[], meta_path=[], path_hooks=[], editable={},
        files=[], dont_write_bytecode=True, coverage=dict(missing=[], changed=[], all_covered=True),
        scientific_origin=dict(root="/fixture/core"))
    prior = store(tmp_path / "prior_lookup.json", expected)
    baseline = store(tmp_path / "baseline.json", dict(tool_entrypoints=dict(
        orthohmm_python=dict(absolute_path="/fixture/python"))))
    lookup_binding = store(tmp_path / "lookup_binding.json", dict(runtime_specs=specs))
    lookup = store(tmp_path / "lookup_receipt.json", dict(
        status="native_lookup_repeated_identity_match", binding=lookup_binding, baseline=baseline,
        source=record(Path(module.__file__).with_name("inspect_native_python_lookup.py")),
        interpreters={name: dict(reports=[prior]) for name in ("orthohmm", "orthofinder")}))
    checker = module.ProspectiveRuntimeChecker(dict(runtime_lookup=lookup, baseline=baseline),
        runtime_ref, tmp_path / "actual_checker")
    calls, launches = [], []
    observed = copy.deepcopy(expected)

    def probe_stub(command, **kwargs):
        calls.append(command)
        destination = Path(command[-1])
        destination.mkdir()
        for name in ("orthohmm", "orthofinder"):
            store(destination / f"{name}.json", observed)
        return SimpleNamespace(returncode=0)

    checker.checker.runner = probe_stub

    def measurement(_):
        launches.append(True)
        if drift == "runtime":
            (roots[0] / "lftp").write_text("changed after measurement", encoding="ascii")
        else:
            observed["modules"] = dict(unexpected="/outside/runtime")
        return dict(status="command_exited_zero", native=dict(exit_code=0))

    result = run_checked(specs, tmp_path / "wrapper", measurement, checker)
    assert launches == [True]
    assert result["status"] == "runtime_changed_or_unverifiable"
    assert result["measurement"]["native"]["exit_code"] == 0
    assert result["scientific_results_admitted"] is False
    assert result["before"]["status"] == "runtime_and_lookup_checked"
    assert result["before"]["runtime"][0]["original_inventory_equality"] is False
    assert len(calls) == (1 if drift == "runtime" else 2)
    assert not (checker.output / "checked_02.json").exists()
