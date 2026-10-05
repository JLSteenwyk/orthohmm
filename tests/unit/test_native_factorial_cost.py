"""Control/input/resource semantics; no production native inference in tests."""

from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_native_factorial_cost as cost
from benchmark_tools.native_factorial_adapter import factors
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def plan_fixture():
    panel = cost.ROOT / "benchmarks/results/test_factorial_cost_uncreated"
    runs = []
    for i, (dataset, cell) in enumerate(cost.IDENTITIES):
        inputs = [dict(path=f"/input/{dataset}/s{j:02d}.fa", bytes=9, sha256="a" * 64)
                  for j in range(12 if i < 6 else 78)]
        names = [Path(r["path"]).name for r in inputs]
        runs.append(dict(index=i, dataset=dataset, cell=cell, repeat=0, inputs=inputs,
            input_directory=f"/input/{dataset}", native_order=names, input_creation_order=names,
            genes=251378 if i < 6 else 984137, proteomes=len(inputs), output_root=str(panel / f"run_{i:02d}")))
    return dict(schema="native_factorial_cost_plan_v1", root=str(cost.ROOT), panel_root=str(panel),
        core_commit=cost.CORE_COMMIT, execution_scope=cost.SCOPE, automatic_retry=False, runs=runs,
        resources=dict(native_cpu_ids=list(range(32)), slurm_slots=64, memory_bytes=cost.MEMORY,
            timeout_s=85800, sample_period_s=1., host_period_s=30., minimum_available_memory_bytes=cost.MEMORY),
        helper_sources=[dict(path=str(cost.ROOT / "benchmark_tools" / name)) for name in (
            "run_native_factorial_cost.py", "prepare_native_factorial_cost.py", "native_factorial_adapter.py", "run_native_factorial_cost.sh")])


def test_exact_thirteen_identities_and_no_completed_cells():
    plan = plan_fixture()
    assert cost.validate_plan(plan) == plan["runs"]
    assert len(cost.IDENTITIES) == 13
    assert ("orthobench", "p1_c0_r0") not in cost.IDENTITIES
    assert ("orthobench", "p1_c1_r1") not in cost.IDENTITIES
    assert ("qfo_corrected", "p1_c0_r0") not in cost.IDENTITIES


@pytest.mark.parametrize("key,value", [("schema", "old_scaling_panel"), ("root", "/tmp"),
    ("execution_scope", "isolated_controlled"), ("automatic_retry", True), ("core_commit", "changed")])
def test_plan_scope_mutation(key, value):
    plan = plan_fixture()
    plan[key] = value
    with pytest.raises(ValueError):
        cost.validate_plan(plan)


@pytest.mark.parametrize("key,value", [("memory_bytes", 96 * 1024**3), ("native_cpu_ids", list(range(1, 33))),
    ("slurm_slots", 32), ("timeout_s", 3600), ("sample_period_s", 5.), ("host_period_s", 60.),
    ("minimum_available_memory_bytes", 0)])
def test_resource_mutation(key, value):
    plan = plan_fixture()
    plan["resources"][key] = value
    with pytest.raises(ValueError):
        cost.validate_plan(plan)


@pytest.mark.parametrize("key,value", [("index", True), ("cell", "p1_c0_r0"), ("dataset", "other"),
    ("repeat", 1), ("genes", 10), ("proteomes", 4), ("output_root", "/tmp/output")])
def test_identity_mutation(key, value):
    plan = plan_fixture()
    plan["runs"][0][key] = value
    with pytest.raises(ValueError):
        cost.validate_plan(plan)


def test_duplicate_or_unbound_sources_and_inputs():
    plan = plan_fixture()
    plan["helper_sources"].pop()
    with pytest.raises(ValueError):
        cost.validate_plan(plan)
    plan = plan_fixture()
    plan["runs"][0]["inputs"][1] = plan["runs"][0]["inputs"][0]
    with pytest.raises(ValueError):
        cost.validate_plan(plan)


def request_fixture():
    ref = dict(path="/plan.json", bytes=123, sha256="b" * 64)
    return dict(schema="native_factorial_cost_request_v1", execution_authorized=True, job_id=123,
        plan=ref, scheduler_command=str(cost.SCRIPT), allocation_cwd=str(cost.ROOT), index=0, history=[]), ref


def test_request_exact_job_plan_history():
    request, ref = request_fixture()
    cost.validate_request(request, ref, 123)
    for key, value in (("job_id", 124), ("index", True), ("index", 13), ("execution_authorized", False),
                       ("history", [dict(job_id=12)]), ("allocation_cwd", "/tmp"), ("scheduler_command", "/other.sh")):
        changed = dict(request, **{key: value})
        with pytest.raises(ValueError):
            cost.validate_request(changed, ref, 123)


@pytest.mark.parametrize("p", (0, 1))
@pytest.mark.parametrize("c", (0, 1))
@pytest.mark.parametrize("r", (0, 1))
def test_all_factor_native_arguments_and_stages(p, c, r):
    run = dict(cell=f"p{p}_c{c}_r{r}", output_root="/native/run_00", genes=8)
    module = SimpleNamespace(StopStep=SimpleNamespace(infer="infer"), SubstitutionMatrix=SimpleNamespace(blosum62="BLOSUM62"))
    baseline = dict(tool_entrypoints={n: dict(absolute_path="/tools/" + n) for n in ("mafft", "FastTree")})
    args = cost.native_kwargs(module, run, baseline)
    assert args["accuracy_profile"] == "high_sensitivity"
    assert args["start"] is None and args["stop"] == "infer"
    assert args["cpu"] == 32 and args["threads_per_worker"] == 4
    assert args["evalue_threshold"] == .0001 and args["substitution_matrix"] == "BLOSUM62"
    assert args["phylogeny"] == ("reconcile" if r else "off")
    assert args["phylogeny_candidates"] == ("satellite_v2" if c else "seed")
    flags = factors(run["cell"])
    stages = {"search", "edge_thresholds", "network_edges", "clustering", "refinement", "orthogroup_materialization"}
    for enabled, stage in ((p, "profile_expansion"), (c, "phylogeny_candidates"), (r, "phylogeny")):
        if enabled:
            stages.add(stage)
    metrics = dict(status="complete", metadata=dict(native_factorial=flags), counts=dict(genes=8,
        phylogeny_checkpoint_hits=0, phylogeny_species_tree_checkpoint_hit=False), stages=dict.fromkeys(stages))
    cost.verify_metrics(metrics, run, flags)
    wrong = deepcopy(metrics)
    wrong["stages"].pop("search")
    with pytest.raises(ValueError):
        cost.verify_metrics(wrong, run, flags)
    if r:
        wrong = deepcopy(metrics)
        wrong["counts"]["phylogeny_checkpoint_hits"] = 1
        with pytest.raises(ValueError):
            cost.verify_metrics(wrong, run, flags)


def preparation_fixture(tmp_path):
    core = tmp_path / "core"
    (core / "orthohmm").mkdir(parents=True)
    enumerator = core / "orthohmm/files.py"
    enumerator.write_text("import glob, os\ndef fetch_fasta_files(directory):\n    return [os.path.basename(p) for p in glob.glob(directory + '/*.fa')]\n")
    pin = record(enumerator)
    baseline = dict(core_root=str(core), core_sources=[dict(absolute_path=pin["path"], bytes=pin["bytes"], sha256=pin["sha256"])])
    source = tmp_path / "original"
    source.mkdir()
    for i in range(4):
        (source / f"s{i}.fa").write_text(f">g{i}a description\nAAAA\n>g{i}b\nCCCC\n")
    root = tmp_path / "run_00"
    root.mkdir()
    names = [f"s{i}.fa" for i in range(4)]
    run = dict(output_root=str(root), inputs=[record(source / n) for n in names],
               input_creation_order=names, native_order=[p.name for p in source.iterdir()], genes=8, proteomes=4)
    return run, baseline


def test_actual_fresh_copy_ownership_and_enum(tmp_path):
    run, baseline = preparation_fixture(tmp_path)
    result = cost.prepare_inputs(run, baseline)
    assert result["genes"] == 8
    assert result["per_species_counts"] == dict.fromkeys(run["native_order"], 2)
    assert result["input_snapshot"]["datasets"][0]["native_order"] == run["native_order"]
    assert result["cold_cache_claim"] is False
    for pin in run["inputs"]:
        assert Path(pin["path"]).read_bytes() == (Path(run["output_root"]) / "input" / Path(pin["path"]).name).read_bytes()
    with pytest.raises(FileExistsError):
        cost.prepare_inputs(run, baseline)
    changed = Path(run["output_root"]) / "input/s0.fa"
    changed.write_text(changed.read_text() + "\n")
    with pytest.raises(ValueError, match="bytes changed"):
        cost.prepared_inputs(run, baseline)


def test_duplicate_ownership_fails_and_preserves_copy(tmp_path):
    run, baseline = preparation_fixture(tmp_path)
    original = Path(run["inputs"][1]["path"])
    original.write_text(">g0a\nAAAA\n>other\nCCCC\n")
    run["inputs"][1] = record(original)
    with pytest.raises(ValueError, match="Duplicate gene ownership"):
        cost.prepare_inputs(run, baseline)
    assert (Path(run["output_root"]) / "input/s0.fa").exists()


def test_enumeration_mismatch_does_not_recopy_or_infer(tmp_path):
    run, baseline = preparation_fixture(tmp_path)
    run["native_order"] = list(reversed(run["native_order"]))
    with pytest.raises(ValueError, match="enumeration differs"):
        cost.prepare_inputs(run, baseline)
    assert not (Path(run["output_root"]) / "native").exists()
    assert len(list((Path(run["output_root"]) / "input").iterdir())) == 4


@pytest.mark.parametrize("memory,comment,matched,success", [
    (cost.MEMORY, "request_sha", True, True),
    (cost.MEMORY - 1024, "request_sha", True, False),
    (cost.MEMORY, "other_request", True, False),
    (cost.MEMORY, "request_sha", False, False)])
def test_release_capacity_attribution_and_scheduler_request(tmp_path, monkeypatch, memory, comment, matched, success):
    from benchmark_tools import probe_host_counters, observe_threadripper_process_identity
    from benchmark_tools import review_threadripper_process_policy, slurm_resource_snapshot, verify_threadripper_controller
    plan_ref = dict(path="/plan.json", bytes=1, sha256="plan_sha")
    request_ref = dict(path="/request.json", bytes=1, sha256="request_sha")
    policy_ref = dict(path="/policy.json", bytes=1, sha256="policy_sha")
    process_ref = dict(path="/process.json", bytes=1, sha256="process_sha")
    request = dict(job_id=123)
    plan = dict(helper_sources=[])
    policy = dict(process_policy=process_ref)
    refs = {"/request.json": request, "/policy.json": policy, "/process.json": dict(schema="typed")}
    monkeypatch.setattr(cost, "read", lambda ref: refs[ref["path"]])
    monkeypatch.setattr(cost, "check", lambda ref: None)
    monkeypatch.setattr(cost, "prepared_inputs", lambda *args: None)
    monkeypatch.setattr(probe_host_counters, "snapshot", lambda: dict(cpu="diagnostic"))
    monkeypatch.setattr(observe_threadripper_process_identity, "enriched_snapshot", lambda: dict(sequence=2))
    monkeypatch.setattr(review_threadripper_process_policy, "review", lambda *args, **kw: dict(process_policy_matched=matched))
    monkeypatch.setattr(slurm_resource_snapshot, "scoped_path", lambda *args: Path("/slurm/job_123/step_0/user/task_0"))
    class Budget:
        def __init__(self, *args, **kwargs):
            assert kwargs["allocation_mode"] == "shared"
        def __call__(self, directory):
            cost.save(directory / "release_budget.json", dict(test_only=True))
            return dict(allocation=dict(fields=dict(Comment=comment)))
    monkeypatch.setattr(verify_threadripper_controller, "ReleaseBudgetGuard", Budget)
    original_read = Path.read_text
    monkeypatch.setattr(Path, "read_text", lambda path, *a, **kw:
        f"MemAvailable: {memory // 1024} kB\n" if str(path) == "/proc/meminfo" else original_read(path, *a, **kw))
    (tmp_path / "ready.json").write_text(json.dumps(dict(cgroup="test")))
    (tmp_path / "host_processes.jsonl").write_text(json.dumps(dict(snapshot=dict(sequence=1), observer_pid=99)) + "\n")
    guard = cost.ReleaseGuard(request_ref, plan_ref, plan, dict(index=0), {}, policy_ref)
    if success:
        guard(tmp_path)
        preflight = json.loads((tmp_path / "environment_preflight.json").read_text())
        assert preflight["background_competition_recorded"] is True
        assert preflight["uncontended_timing"] is False
        assert preflight["plan_sha256"] == "plan_sha"
    else:
        with pytest.raises(ValueError):
            guard(tmp_path)
        assert not (tmp_path / "environment_preflight.json").exists()
    assert (tmp_path / "launch_environment_observation.json").exists()


@pytest.mark.parametrize("raw", ("MemAvailable: -1 kB", "MemAvailable: 200 MB", "MemTotal: 5 kB", "MemAvailable: 12"))
def test_invalid_capacity(raw):
    with pytest.raises(ValueError):
        cost.available_memory(raw)


def test_capacity_is_not_a_foreign_cpu_gate():
    assert cost.available_memory("MemTotal: 9999 kB\nMemAvailable: 1048576 kB\n") == 1024**3
    source = (cost.ROOT / "benchmark_tools/run_native_factorial_cost.py").read_text()
    assert "capacity < MEMORY or not verdict[\"process_policy_matched\"]" in source
    assert "maximum_foreign_average_cores" not in source
    assert "scancel" not in source and "os.kill" not in source


def test_launch_bootstrap_uses_separate_request_scope():
    script = cost.SCRIPT.read_text()
    assert "#SBATCH --cpus-per-task=64" in script and "#SBATCH --mem=128G" in script
    assert "#SBATCH --time=1-02:00:00" in script and "#SBATCH --no-requeue" in script
    assert "#SBATCH --exclusive" not in script
    assert "unset PYTHONHOME PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT" in script
    assert "scheduler-comment" in script
    assert "run_native_factorial_cost.py" in script and "run_threadripper_scaling.py" not in script


def test_actual_prepared_plan_and_policy_bytes():
    import hashlib
    import subprocess
    path = cost.ROOT / "benchmark_tools/results/native_factorial_cost_plan_20261004.json"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "c697d9c9b08df1c0ad40807e6db6044ffd597dc3dc997debf60fbcf492d45784"
    plan = json.loads(path.read_text())
    assert len(cost.validate_plan(plan)) == 13
    assert plan["scientific_execution_authorized"] is False
    assert plan["source_commit"] == "3c275a025b6e3743beffe08db6d48a54454e4732"
    for ref in plan["helper_sources"] + plan["evidence"]:
        raw = Path(ref["path"]).read_bytes()
        if Path(ref["path"]).name in {"prepare_native_factorial_cost.py", "run_native_factorial_cost.py"}:
            # The failed attempt keeps its historical generator, not a rewritten pin.
            relative = Path(ref["path"]).relative_to(cost.ROOT).as_posix()
            raw = subprocess.check_output(["git", "-C", str(cost.ROOT), "show", plan["source_commit"] + ":" + relative])
        assert len(raw) == ref["bytes"]
        assert hashlib.sha256(raw).hexdigest() == ref["sha256"]
    for index in (0, 6):
        run = plan["runs"][index]
        names = {Path(ref["path"]).name for ref in run["inputs"]}
        assert names == set(run["native_order"]) == set(run["input_creation_order"])
        assert run["native_order"] == plan["filename_enumeration_probes"][run["dataset"]]["datasets"][0]["native_order"]
    policy_path = cost.ROOT / "benchmark_tools/results/native_factorial_cost_policy_20261004.json"
    assert hashlib.sha256(policy_path.read_bytes()).hexdigest() == "5a7ef6f0b0b4f201ca41376a1bc3c2c652ac916e5b6fbf72f803765151954c31"
    policy = json.loads(policy_path.read_text())
    assert policy["plan_sha256"] == "c697d9c9b08df1c0ad40807e6db6044ffd597dc3dc997debf60fbcf492d45784"
    assert policy["foreign_cpu_role"] == policy["native_pressure_role"] == policy["preflight_pressure_role"] == "diagnostic_only"
    assert policy["minimum_available_memory_bytes"] == cost.MEMORY


def test_exact_packaging_only_runtime_delta():
    from benchmark_tools.refresh_native_factorial_lookup import delta, CHANGES
    before = dict(roots=["/unchanged"], records=[])
    after = deepcopy(before)
    for name, (a, b) in CHANGES.items():
        row = dict(path=str(cost.ROOT / "benchmark_tools" / name), mode=436, kind="file", bytes=10, sha256=a)
        before["records"].append(row)
        after["records"].append(dict(row, sha256=b, bytes=20))
    assert len(delta(before, after)) == 3
    for mutation in ("add", "hash", "mode", "roots", "other"):
        wrong = deepcopy(after)
        if mutation == "add":
            wrong["records"].append(dict(path="/unexpected", mode=436, kind="file", bytes=4, sha256="bad"))
        elif mutation in ("hash", "mode"):
            wrong["records"][0]["sha256" if mutation == "hash" else "mode"] = "bad" if mutation == "hash" else 0
        elif mutation == "roots":
            wrong["roots"] = ["/other"]
        else:
            wrong["unreviewed_metadata"] = True
        with pytest.raises(ValueError):
            delta(before, wrong)


def test_actual_22426_failed_before_materialization_or_inference():
    import hashlib
    root = cost.ROOT / "benchmarks/results/native_factorial_cost_v1_20261004"
    report = root / "sessions/run_00/result.json"
    result = json.loads(report.read_text())
    assert result["job_id"] == 22426 and result["index"] == 0
    assert result["status"] == "factorial_attempt_failed_retained"
    assert result["wrapper"]["status"] == "verified_wrapper_failed"
    assert result["wrapper"]["error"] == "Runtime inventory changed"
    assert result["automatic_retry"] is False
    assert result["next_identity_authorized"] is False
    assert not (root / "run_00/input").exists()
    assert not (root / "run_00/native").exists()
    assert not (root / "run_00/measurement").exists()
    assert hashlib.sha256(Path(result["plan"]["path"]).read_bytes()).hexdigest() == result["plan"]["sha256"]


def test_expired_controller_accounting_identity():
    raw = "22426|FAILED|1:0|64|128G|bizon|gpu|1560|2026-10-04T18:52:24|2026-10-04T18:52:59|2026-10-04T18:53:27|orthohmm_factorial_cost\n"
    assert cost.terminal_accounting(raw, 22426)["State"] == "FAILED"
    for a, b in (("22426", "22426_0"), ("FAILED", "RUNNING"), ("64|128G", "32|128G"),
                 ("128G", "96G"), ("1560", "1440"), ("bizon", "dgx")):
        with pytest.raises(ValueError):
            cost.terminal_accounting(raw.replace(a, b), 22426)


def test_transient_controller_error_does_not_imply_terminal(monkeypatch):
    import subprocess
    calls = []
    def failed(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 1, "", "Connection refused")
    monkeypatch.setattr(cost.subprocess, "run", failed)
    with pytest.raises(ValueError, match="not an expired"):
        cost.verify_terminal(22426)
    assert len(calls) == 1


def test_actual_repaired_plan_preserves_factors_and_failed_attempt():
    import hashlib
    path = cost.ROOT / "benchmark_tools/results/native_factorial_cost_plan_repaired_20261004.json"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "61756d74d1778fd268c9c4fe4f85b32cc130121f9bf02793f470b202bc4c2b1a"
    plan = json.loads(path.read_text())
    cost.validate_plan(plan)
    old = json.loads(Path(plan["supersedes_plan"]["path"]).read_text())
    for i, run in enumerate(plan["runs"]):
        for key in ("index", "dataset", "cell", "repeat", "inputs", "genes", "proteomes", "native_order", "input_creation_order"):
            assert run[key] == old["runs"][i][key]
        assert run["output_root"] != old["runs"][i]["output_root"]
    for ref in plan["helper_sources"] + plan["evidence"]:
        source = Path(ref["path"])
        if source == cost.ROOT / "benchmark_tools/run_native_factorial_cost.py":
            # Verify the historical executor, not the repaired current source.
            source = cost.ROOT / "benchmark_tools/results/native_factorial_failed_source_22427/run_native_factorial_cost.py"
        raw = source.read_bytes()
        assert len(raw) == ref["bytes"] and hashlib.sha256(raw).hexdigest() == ref["sha256"]
    failure = json.loads(Path(plan["retained_preflight_failure"]["path"]).read_text())
    assert failure["job_id"] == 22426 and failure["status"] == "factorial_attempt_failed_retained"
    lookup = json.loads(Path(plan["runtime_lookup"]["path"]).read_text())
    prior = json.loads(Path(old["runtime_lookup"]["path"]).read_text())
    assert lookup["baseline"] == prior["baseline"] == plan["baseline"]
    a, b = (json.loads(Path(r["binding"]["path"]).read_text()) for r in (prior, lookup))
    for key in ("baseline", "controller_python", "command_plan", "baseline_paths", "retired_roots"):
        assert a[key] == b[key]
    assert a["runtime_specs"][1:] == b["runtime_specs"][1:]
    delta = json.loads(Path(lookup["current_validation"]["delta"]["path"]).read_text())
    assert len(delta["changed"]) == 3 and delta["added"] == delta["removed"] == []
    for side in ("before", "after"):
        assert all(row["status"] == "runtime_tree_identity_matches" for row in lookup["current_validation"]["runtime_checks"][side])


@pytest.mark.parametrize("status", ["native_factorial_completed_pending_output_review", "native_factorial_failed"])
def test_terminal_native_receipt_updates_owned_running_file(tmp_path, status):
    path = tmp_path / "native_execution.json"
    initial = dict(status="native_factorial_running", index=0, cell="p0_c0_r0")
    cost.save(path, initial)
    terminal = dict(initial, status=status, counts=dict(genes=2))
    cost.update_native_receipt(path, initial, terminal)
    assert json.loads(path.read_text()) == terminal
    assert not path.with_name(path.name + ".terminal.pending").exists()
    with pytest.raises(ValueError, match="changed"):
        cost.update_native_receipt(path, initial, terminal)


def test_receipt_update_does_not_overwrite_collision_or_foreign_change(tmp_path):
    path = tmp_path / "native_execution.json"
    initial = dict(status="native_factorial_running")
    cost.save(path, initial)
    pending = path.with_name(path.name + ".terminal.pending")
    pending.write_text("foreign evidence")
    with pytest.raises(FileExistsError):
        cost.update_native_receipt(path, initial, dict(status="complete"))
    assert json.loads(path.read_text()) == initial
    assert pending.read_text() == "foreign evidence"
    path.write_text(json.dumps(dict(status="externally_changed")))
    with pytest.raises(ValueError, match="changed"):
        cost.update_native_receipt(path, initial, dict(status="complete"))


def test_receipt_update_rejects_symlink(tmp_path):
    actual = tmp_path / "actual.json"
    initial = dict(status="native_factorial_running")
    actual.write_text(json.dumps(initial))
    link = tmp_path / "native_execution.json"
    link.symlink_to(actual)
    with pytest.raises(ValueError, match="changed"):
        cost.update_native_receipt(link, initial, dict(status="complete"))
    assert json.loads(actual.read_text()) == initial
