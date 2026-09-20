"""Bind each overhead arm to its frozen native task, runtime and report bytes."""

from copy import deepcopy
import json
import math
from pathlib import Path
import re

from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_overhead import ROOT, build
from benchmark_tools.run_root_context_overhead import PLAN_SHA
from benchmark_tools.verify_lineage_native_provenance import same

RECIPE_ROOT = ROOT / "root_context_overhead_recipe_v1"
RECIPE_PATH = ROOT / "root_context_overhead_recipe_v1.json"
SUBMISSION_SCRIPT = RECIPE_ROOT / "benchmark_tools/run_dgx_root_context_overhead.sh"


def load_context(results, recipe_path, recipe_sha):
    if not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Invalid recipe hash")
    plan = read_pinned(results / "dgx_root_context_overhead_plan_20260919.json", PLAN_SHA)
    expected = build(results / "dgx_root_context_native_plan_20260919.json",
                     results / "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md")
    if not same(plan, expected):
        raise ValueError("Frozen overhead plan derivation differs")
    recipe = read_pinned(recipe_path, recipe_sha)
    paths = [row["path"] for row in recipe["records"]]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate recipe record")
    return dict(plan=plan, recipe=recipe, recipe_sha256=recipe_sha)


def scheduler_identity(text, job, require_completed=True):
    if type(job) is not int or job <= 0:
        raise ValueError("Invalid job identity")
    fields = {}
    for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text):
        if key in fields:
            raise ValueError("Duplicate scheduler field")
        fields[key] = value
    if terminal_record(text, job) is None:
        raise ValueError("Require terminal overhead allocation")
    expected = dict(JobId=str(job), Restarts="0", Requeue="0", NodeList="spark-7ff0", Partition="spark",
        OverSubscribe="NO", MinMemoryNode="96G", NumNodes="1", NumCPUs="20", NumTasks="1",
        TimeLimit="05:00:00", Command=str(SUBMISSION_SCRIPT), WorkDir=str(RECIPE_ROOT))
    expected["CPUs/Task"] = "20"
    if require_completed:
        expected.update(JobState="COMPLETED", ExitCode="0:0")
    if any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError("Overhead scheduler identity, allocation or outcome differs")
    return fields


def verify(context, index, directory, receipt, scheduler, job, require_completed=True):
    scheduler_identity(scheduler, job, require_completed)
    if type(index) is not int or index not in range(18):
        raise ValueError("Invalid overhead task index")
    plan, recipe = context["plan"], context["recipe"]
    task = plan["runs"][index]
    remote = Path(task["run"]["measurement_directory"]).parent
    paths = [directory / name for name in ("preparation.json", "verification.json", "measurement/lineage_report.json")]
    supplementary = directory / "measurement/root_context_report.json"
    if task["arm"] == "root_context":
        paths.append(supplementary)
    elif supplementary.exists() or supplementary.is_symlink():
        raise ValueError("Supplementary report in lineage arm")
    if directory.is_symlink() or any(path.is_symlink() for path in paths):
        raise ValueError("Symlinked task evidence")
    evidence = [record(path) for path in paths]
    preparation, verification, measured, *root_report = [json.loads(path.read_text()) for path in paths]
    expected = {key: task[key] for key in ("index", "block", "pair", "arm", "method")}
    expected.update(task=task, job_id=job, scientific_timings_admitted=False,
        status="measurement_completed", wrapper_status="command_exited_zero", native_wall_s=measured["native_wall_s"],
        verification=dict(evidence[1], path=str(remote / "verification.json")),
        lineage_report=dict(evidence[2], path=str(remote / "measurement/lineage_report.json")))
    if root_report:
        expected["root_context_report"] = dict(evidence[3], path=str(remote / "measurement/root_context_report.json"))
        if not same(root_report[0]["lineage_report"], dict(evidence[2], path="lineage_report.json")):
            raise ValueError("Supplementary report refers to different lineage bytes")
    if not same(receipt, expected):
        raise ValueError("Task receipt, arm or report bytes differ")
    if type(measured["job_id"]) is not int or measured["job_id"] != job:
        raise ValueError("Measured job differs")
    if not measured["points"] or any(("root_context" in point) != (task["arm"] == "root_context") for point in measured["points"]):
        raise ValueError("Observation fields do not match collector arm")
    native = measured["native"]
    if type(native["exit_code"]) is not int or native["exit_code"] != 0 or native["timed_out"] is not False:
        raise ValueError("Native command failed or timed out")
    run = deepcopy(task["run"])
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(remote / "native.time.tsv"))
    argv = time_command(run["native_argv"], run["gnu_time"]["output"])
    order = plan["order"]
    originals = order["inputs_in_native_order"]
    copies = []
    if run["native_method"] == "orthofinder_full":
        copies = [dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name))
                  for item in sorted(originals, key=lambda item: Path(item["path"]).name)]
    expected = dict(run=run, measured_argv=argv, original_inputs=originals, copied_inputs=copies,
        expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False)
    if any(not same(preparation[k], v) for k, v in expected.items()):
        raise ValueError("Prepared native command or input inventory differs")
    if len(plan["runtime_manifests"]) != 2:
        raise ValueError("Unexpected runtime inventory")
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(RECIPE_PATH), sha256=context["recipe_sha256"], records=len(recipe["records"]),
                       scientific_execution_authorized=False, status="runtime_tree_identity_matches"))
    before = dict(native_order=order["native_order"], original_inputs=originals, runtime=runtime)
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = copies
    wrappers = [row["sha256"] for row in recipe["records"] if row["kind"] == "file"
                and row["path"] == str(RECIPE_ROOT / "benchmark_tools/run_verified_slurm_measurement.py")]
    if len(wrappers) != 1:
        raise ValueError("Missing or duplicate wrapper")
    expected = dict(status="command_exited_zero", before=before, after=after, measurement=measured,
                    source_sha256=wrappers[0], scientific_results_admitted=False)
    if any(not same(verification[k], v) for k, v in expected.items()):
        raise ValueError("Runtime, wrapper or embedded measurement differs")
    for value in (preparation["preparation_wall_s"], verification["before_check_wall_s"],
                  verification["after_check_wall_s"], measured["native_wall_s"]):
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid measurement or verification duration")
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
        plan["launcher_python"], "-B", str(RECIPE_ROOT / "benchmark_tools/measure_native_lineage_step.py"),
        "--worker", run["measurement_directory"]]
    if not same(measured["launched"], launched):
        raise ValueError("Native worker command differs")
    for item in evidence:
        check(item)
    return dict(status="root_context_overhead_task_bound", index=index, arm=task["arm"], job_id=job,
        run=run, measured_argv=argv, measurement=measured, evidence=evidence, scientific_timings_admitted=False,
        limitations=["Requires separate complete-panel, recipe archive and waiting-session checks.",
            "Report byte binding and arm fields do not replace raw collector replay or native-output validation.",
            "Before/after runtime checks cannot exclude temporary changes during execution."])
