"""Execute a pinned corrected-input primary command; admission stays separate."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.prepare_qfo_corrected_primary import primary_commands, BASELINE_SHA

METHODS = ("orthohmm_high_sensitivity", "orthofinder_full")


def verify_child_resolution(runtime, env, python, cwd):
    inspector = runtime["source"]
    check(inspector)
    check(runtime["python_executable"])
    if record(python) != runtime["python_executable"]:
        raise ValueError("Changed OrthoFinder Python executable")
    code = (
        "import importlib.util,json,sys; "
        "spec=importlib.util.spec_from_file_location('runtime_inspector',sys.argv[1]); "
        "module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); "
        "print(json.dumps(module.inspect()))"
    )
    observed = subprocess.run(
        [str(python), "-c", code, inspector["path"]],
        env=env, cwd=cwd, capture_output=True, text=True, check=True, timeout=240,
    )
    observed = json.loads(observed.stdout)
    for field in ("python_executable", "package_sources", "packages"):
        if observed[field] != runtime[field]:
            raise ValueError(f"Changed OrthoFinder runtime: {field}")
    for name, expected in runtime["child_tools"].items():
        tool = observed["child_tools"].get(name, {})
        if tool.get("status") != "resolved" or tool.get("file") != expected["file"] or tool.get("exit_code") != 0:
            raise ValueError(f"Changed OrthoFinder child-tool resolution: {name}")


def verify(plan):
    if plan["status"] != "corrected_primary_commands_frozen_unrun" or plan["execution_authorized"] is not False:
        raise ValueError("Unexpected command-plan state")
    baseline = read_frozen(Path(plan["baseline"]["path"]), BASELINE_SHA)
    expected = primary_commands(Path(plan["input_directory"]), Path(plan["output_root"]), baseline)
    if plan["methods"] != expected or plan["resource_plan"]["cpus"] != 32:
        raise ValueError("Changed native commands or CPU allocation")
    for item in [plan["source"], plan["baseline"], plan["runtime"], *plan["helper_sources"], *plan["inputs"]]:
        check(item)
    directory = Path(plan["input_directory"])
    stage = json.loads((directory / "staging_manifest.json").read_text())
    inputs = stage["input_fastas"]
    if {p.name for p in directory.iterdir()} != {Path(r["path"]).name for r in inputs} | {"staging_manifest.json"}:
        raise ValueError("Changed input directory contents")
    for item in inputs:
        check(item)
    verify_environment(baseline)
    env, resolved = execution_environment(baseline)
    runtime = json.loads(Path(plan["runtime"]["path"]).read_text())
    for tool in runtime["child_tools"].values():
        check(tool["file"])
    for item in runtime["package_sources"]:
        check(item)
    if resolved != plan["resolved_outer_executables"]:
        raise ValueError("Changed executable resolution")
    # Preserve the virtual-environment entrypoint, not its resolved system symlink.
    python = Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python"
    verify_child_resolution(runtime, env, python, baseline["core_root"])
    return env, inputs


def run(manifest, expected_sha, index):
    if type(index) is not int or index not in (0, 1):
        raise ValueError("Require primary task index0or1")
    if os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled32CPU allocation")
    plan = read_frozen(manifest, expected_sha)
    env, inputs = verify(plan)
    method = METHODS[index]
    config = plan["methods"][method]
    directory = Path(config["output"])
    execution = Path(plan["output_root"]) / "execution" / method
    if directory.exists() or execution.exists():
        raise FileExistsError("Existing native or execution output; no implicit resume")
    execution.mkdir(parents=True, exist_ok=False)
    directory.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "method": method, "index": index, "manifest": record(manifest),
              "source": record(__file__), "job_id": os.environ["SLURM_JOB_ID"],
              "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
              "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
              "node": os.uname().nodename, "native_argv": config["native_argv"], "cwd": config["cwd"],
              "accuracy_admitted": False, "native_outputs_validated": False,
              "limitation": "Shared-host inference; not dedicated matched timing. Conversion and scoring require independent admission."}
    status = execution / "status.json"
    status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    try:
        if method == "orthofinder_full":
            copied = Path(config["copy_inputs_to"])
            copied.mkdir(exist_ok=False)
            for item in inputs:
                source = Path(item["path"])
                if source.parent != Path(config["copy_inputs_from"]):
                    raise ValueError("Unexpected input copy source")
                target = copied / source.name
                shutil.copy2(source, target)
                observed = record(target)
                if any(observed[k] != item[k] for k in ("bytes", "sha256")):
                    raise ValueError("Input copy differs")
            report["copied_inputs"] = [record(p) for p in sorted(copied.iterdir())]
        report.update(status="running", started_epoch=time.time())
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (execution / "native.log").open("xb") as log:
            result = subprocess.run(["/usr/bin/time", "-v", "-o", str(execution / "time.txt"),
                                     *config["native_argv"]], cwd=config["cwd"], env=env,
                                    stdout=log, stderr=subprocess.STDOUT)
        report.update(exit_code=result.returncode, finished_epoch=time.time())
        verify(plan)
        check(report["manifest"])
        if method == "orthofinder_full":
            for item in report["copied_inputs"]:
                check(item)
        report.update(status="process_succeeded_pending_native_admission" if result.returncode == 0 else "native_process_failed",
                      outputs=[record(p) for p in sorted(directory.rglob("*")) if p.is_file()],
                      log=record(execution / "native.log"), timing=record(execution / "time.txt"))
        if result.returncode:
            raise RuntimeError(f"Native process failed: {result.returncode}")
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--index", type=int, choices=(0, 1), required=True)
    args = parser.parse_args()
    run(args.manifest.resolve(), args.manifest_sha256, args.index)
