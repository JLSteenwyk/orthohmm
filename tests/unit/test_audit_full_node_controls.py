import copy
from pathlib import Path

import pytest

from benchmark_tools import audit_full_node_controls as module


@pytest.fixture
def binding():
    results = Path(module.__file__).parent / "results"
    plan = module.read_pinned(results / "dgx_dual_native_plan_20260919.json", module.PLAN_SHA)
    recipe = module.read_pinned(results / "full_node_controls_recipe_20260919.json", module.RECIPE_SHA)
    rows = [dict(block=block, mode=mode, status="failed", error="fixture")
            for block, modes in enumerate(module.ORDER) for mode in modes]
    panel = dict(job_id=module.JOB, protocol_sha256=module.PROTOCOL_SHA, trials=rows,
                 summary=module.summarize(rows), scientific_timings_admitted=False, publication_ready=False)
    runtime = [dict(row, records=count, status="runtime_tree_identity_matches", scientific_execution_authorized=False)
               for row, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(module.ROOT / (module.RECIPE + ".json")), sha256=module.RECIPE_SHA,
                        records=len(recipe["records"]), status="runtime_tree_identity_matches",
                        scientific_execution_authorized=False))
    wrapper = next(row["sha256"] for row in recipe["records"]
                   if row["path"].endswith("/run_verified_slurm_measurement.py"))
    verification = dict(before=runtime, after=copy.deepcopy(runtime), measurement=panel,
        status="full_node_controls_completed", source_sha256=wrapper, scientific_results_admitted=False)
    launch = dict(job_id=module.JOB, recipe_sha256=module.RECIPE_SHA, protocol_sha256=module.PROTOCOL_SHA,
                  plan_sha256=module.PLAN_SHA, executable=plan["launcher_python"], host="spark-7ff0")
    return plan, recipe, verification, launch, panel


def test_binding_retains_failed_trials(binding):
    module.bind(*binding)


@pytest.mark.parametrize("fault", ["runtime", "source", "launch", "order", "summary", "admission"])
def test_binding_rejects_inconsistency(binding, fault):
    plan, recipe, verification, launch, panel = binding
    if fault == "runtime":
        verification["after"][0]["records"] += 1
    elif fault == "source":
        verification["source_sha256"] = "0"*64
    elif fault == "launch":
        launch["job_id"] += 1
    elif fault == "order":
        panel["trials"].reverse()
    elif fault == "summary":
        panel["summary"]["all_workloads_valid"] = True
    else:
        panel["scientific_timings_admitted"] = True
    with pytest.raises(ValueError):
        module.bind(*binding)


def test_live_scheduler_blocks_archive_access(tmp_path):
    with pytest.raises(ValueError, match="terminal"):
        module.scheduler(f"JobId={module.JOB} JobState=RUNNING")


def test_allocation_check_rejects_missing_fields(monkeypatch):
    monkeypatch.setattr(module, "terminal_record", lambda text, job: text)
    with pytest.raises(ValueError, match="allocation"):
        module.scheduler(f"JobId={module.JOB} JobState=COMPLETED")
