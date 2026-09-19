import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_root_context_controls as module
from benchmark_tools.run_root_context_controls import ORDER, summarize

JOB = 123
MANAGER = "/user.slice/user-1000.slice/user@1000.service"


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def scheduler_text():
    path = Path(module.__file__).parent / "results/full_node_control_scheduler_21918_20260919.txt"
    value = path.read_text().replace("JobId=21918", f"JobId={JOB}").replace("00:45:00", "00:15:00")
    value = value.replace("full_node_controls_recipe_v1", module.RECIPE)
    return value.replace("Command=benchmarks/work/run_dgx_full_node_controls_dd912c6.sh",
        "Command=" + str(module.ROOT / module.RECIPE / "benchmark_tools/run_dgx_root_context_controls.sh"))


@pytest.fixture
def context():
    plan = dict(launcher_python="/pinned/python", runtime_manifests=[dict(path="/runtime/one", sha256="one"),
                                                                 dict(path="/runtime/two", sha256="two")])
    base = module.ROOT / module.RECIPE / "benchmark_tools"
    recipe = dict(records=[dict(kind="file", path=str(base / name), sha256=sha) for name, sha in (
        ("results/ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md", module.PROTOCOL_SHA),
        ("run_verified_slurm_measurement.py", "wrapper"))])
    rows = [dict(block=block, mode=mode, index=index, job_id=JOB, scientific_timings_admitted=False,
                 status="root_context_workload_validated", positive_control_response=True if mode == "user-contended" else None)
        for index, (block, mode) in enumerate((b, m) for b, modes in enumerate(ORDER) for m in modes)]
    panel = dict(status="root_context_controls_completed", trials=rows, job_id=JOB, manager=MANAGER,
        summary=summarize(rows), protocol_sha256=module.PROTOCOL_SHA, scientific_timings_admitted=False, publication_ready=False)
    expected = [dict(row, records=count, status="runtime_tree_identity_matches", scientific_execution_authorized=False)
        for row, count in zip(plan["runtime_manifests"], (26673, 10066))]
    expected.append(dict(path=str(module.ROOT / (module.RECIPE + ".json")), sha256="recipe",
        records=2, status="runtime_tree_identity_matches", scientific_execution_authorized=False))
    verification = dict(before=expected, after=copy.deepcopy(expected), measurement=panel,
        status=panel["status"], source_sha256="wrapper", scientific_results_admitted=False,
        before_check_wall_s=1, after_check_wall_s=1)
    launch = dict(job_id=JOB, recipe_sha256="recipe", protocol_sha256=module.PROTOCOL_SHA,
        runtime_plan_sha256=module.PLAN_SHA, executable=plan["launcher_python"], host="spark-7ff0",
        manager=MANAGER, source_directory=str(module.ROOT / module.RECIPE), scientific_timings_admitted=False)
    return plan, recipe, verification, launch, panel


def bind(context):
    plan, recipe, verification, launch, panel = context
    module.bind(plan, recipe, "recipe", verification, launch, panel, JOB)


def test_valid_bindings_and_scheduler(context):
    bind(context)
    assert module.scheduler(scheduler_text(), JOB)["NumCPUs"] == "20"


def test_session_deployment_is_distinct_and_old_identity_rejected(context):
    plan, recipe, verification, launch, panel = context
    for row in recipe["records"]:
        row["path"] = row["path"].replace("recipe_v1", "recipe_v2")
    for key in ("before", "after"):
        verification[key][-1]["path"] = verification[key][-1]["path"].replace("recipe_v1", "recipe_v2")
    launch["source_directory"] = launch["source_directory"].replace("recipe_v1", "recipe_v2")
    module.bind(plan, recipe, "recipe", verification, launch, panel, JOB, "v2")
    text = scheduler_text().replace("recipe_v1", "recipe_v2").replace("run_dgx_root_context_controls.sh", "run_dgx_root_context_session.sh")
    assert module.scheduler(text, JOB, "v2")["TimeLimit"] == "00:15:00"
    with pytest.raises(ValueError):
        module.scheduler(text, JOB)
    with pytest.raises(ValueError):
        module.scheduler(scheduler_text(), JOB, "v2")


@pytest.mark.parametrize("fault", ["recipe", "before", "after", "summary", "launch", "protocol", "duration", "order", "continued"])
def test_binding_drift_rejected(context, fault):
    plan, recipe, verification, launch, panel = context
    if fault == "recipe":
        recipe["records"].append(copy.deepcopy(recipe["records"][0]))
    elif fault in ("before", "after"):
        verification[fault][0]["records"] = 1
    elif fault == "summary":
        panel["summary"]["all_workloads_valid"] = False
    elif fault == "launch":
        launch["executable"] = "/unverified/python"
    elif fault == "protocol":
        recipe["records"][0]["sha256"] = "changed"
    elif fault == "duration":
        verification["after_check_wall_s"] = float("nan")
    elif fault == "order":
        panel["trials"].reverse()
    else:
        panel["trials"][0]["status"] = "failed"
        panel["summary"] = summarize(panel["trials"])
    with pytest.raises(ValueError):
        bind(context)


@pytest.mark.parametrize("old,new", [("NumCPUs=20", "NumCPUs=19"), ("00:15:00", "00:45:00"),
    ("OverSubscribe=NO", "OverSubscribe=YES"), ("Restarts=0", "Restarts=1"),
    ("JobState=COMPLETED", "JobState=RUNNING"), ("run_dgx_root_context_controls.sh", "other.sh")])
def test_scheduler_drift_rejected(old, new):
    with pytest.raises(ValueError):
        module.scheduler(scheduler_text().replace(old, new), JOB)


@pytest.fixture
def archive(tmp_path, monkeypatch, context):
    plan, recipe, verification, launch, panel = context
    output = tmp_path / module.OUTPUT
    save(output / "verification.json", verification)
    save(output / "launch.json", launch)
    save(output / "measurement/result.json", panel)
    for index, row in enumerate(panel["trials"]):
        path = output / "measurement" / f"trial_{index:02d}"
        save(path / "panel_trial.json", row)
        save(output / "measurement" / f"progress_{index:02d}.json",
             dict(trials=panel["trials"][:index+1], job_id=JOB, failed_index=None))
        if row["mode"] == "user-contended":
            save(path / "controller.json", dict(removal=[dict(monotonic_ns=index*100+90)]))
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text(scheduler_text())
    recipe_path = tmp_path / "recipe.json"
    save(recipe_path, recipe)
    monkeypatch.setattr(module, "recipe_evidence", lambda *args: [])
    monkeypatch.setattr(module, "read_pinned", lambda path, sha: plan if sha == module.PLAN_SHA else recipe)

    def replay(path, remote, source, python, mode, job, index, manager):
        assert (source, python, job, manager) == (module.ROOT / module.RECIPE, plan["launcher_python"], JOB, MANAGER)
        result = {k: v for k, v in panel["trials"][index].items() if k != "block"}
        points = [dict(host=[dict(started_monotonic_ns=index*100+t-1)], root_context=dict(
            boot_before="boot", ticks=100, identities_before={"/": [1, 2]},
            host_after=dict(finished_ns=index*100+t))) for t in (10, 80)]
        return dict(trial=result, measurement=dict(lineage=dict(measured=dict(points=points))))

    monkeypatch.setattr(module, "replay", replay)
    return tmp_path, recipe_path, scheduler, panel


def audit(archive):
    root, recipe, scheduler, _ = archive
    return module.audit(root, recipe, "recipe", scheduler, JOB)


def test_complete_archive_orchestration(archive):
    result = audit(archive)
    assert result["validated_trials"] == 12
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("unexpected_work", [False, True])
def test_stopped_panel_retains_failure_and_unrun_conditions(archive, unexpected_work):
    root, _, scheduler, panel = archive
    output = root / module.OUTPUT
    for index, row in enumerate(panel["trials"]):
        if index == 3:
            row.update(status="failed", error="unconfirmed cleanup")
        elif index > 3:
            row.update(status="not_run_after_failure", failed_index=3)
        path = output / "measurement" / f"trial_{index:02d}"
        if index > 3 and (path / "controller.json").exists():
            (path / "controller.json").unlink()
        save(path / "panel_trial.json", row)
        save(output / "measurement" / f"progress_{index:02d}.json",
             dict(trials=panel["trials"][:index+1], job_id=JOB, failed_index=3 if index >= 3 else None))
    panel.update(status="root_context_controls_stopped", summary=summarize(panel["trials"]))
    save(output / "measurement/result.json", panel)
    verification = json.loads((output / "verification.json").read_text())
    verification.update(status=panel["status"], measurement=panel)
    save(output / "verification.json", verification)
    scheduler.write_text(scheduler.read_text().replace("ExitCode=0:0", "ExitCode=1:0")
                         .replace("JobState=COMPLETED", "JobState=FAILED"))
    if unexpected_work:
        (output / "measurement/trial_04/workload").mkdir()
        with pytest.raises(ValueError, match="Unrun condition"):
            audit(archive)
    else:
        result = audit(archive)
        assert result["validated_trials"] == 3
        assert result["trials"][3]["status"] == "retained_failure"
        assert all(row["status"] == "retained_unrun" for row in result["trials"][4:])


@pytest.mark.parametrize("fault", ["checkpoint", "extra", "row", "exit", "overlap"])
def test_archive_inconsistency_rejected(archive, fault):
    root, _, scheduler, panel = archive
    measurement = root / module.OUTPUT / "measurement"
    if fault == "checkpoint":
        save(measurement / "progress_03.json", {})
    elif fault == "extra":
        (measurement / "trial_12").mkdir()
    elif fault == "row":
        save(measurement / "trial_00/panel_trial.json", {})
    elif fault == "exit":
        scheduler.write_text(scheduler.read_text().replace("ExitCode=0:0", "ExitCode=1:0"))
    else:
        save(measurement / "trial_03/controller.json", dict(removal=[dict(monotonic_ns=10000)]))
    with pytest.raises(ValueError):
        audit(archive)
