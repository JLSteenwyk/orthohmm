"""Freeze and run fresh SonicParanoid on the audited staged Three Kingdoms panel."""

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
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_sonic import command, validate_runtime
from benchmark_tools.run_qfo_corrected_sonic import environment, compare_resolution, INVENTORIES, RUNTIME_SHA
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree
from benchmark_tools.audit_three_kingdoms_method_inputs import SOURCE_SHA

ENTRYPOINT = "/home/bizon/anaconda3/bin/sonicparanoid"


def validate_plan(plan):
    root = Path(plan["output_root"])
    if (not root.is_absolute() or plan["status"] != "matched_three_kingdoms_sonic_frozen_unrun"
            or plan["native_argv"] != command(ENTRYPOINT, root)
            or plan["resources"] != {"node": "bizon", "cpus": 32, "memory_gib": 192, "hours": 72}
            or plan["accuracy_admitted"] is not False or plan["reuse_search"] is not False):
        raise ValueError("Changed matched-input Sonic plan")
    source = read_frozen(Path(plan["input_source_audit"]["path"]), SOURCE_SHA)
    expected = [row["files"]["staged"] for row in source["inputs"]]
    if (plan["inputs"] != expected or len(expected) != 12
            or source["counts"]["proteins"] != 443217
            or source["status"] != "retained_three_kingdoms_lineage_and_reference_verified"):
        raise ValueError("Changed Three Kingdoms input universe")
    directory = Path(expected[0]["path"]).parent
    if {str(p) for p in directory.glob("*.fasta")} != {r["path"] for r in expected}:
        raise ValueError("Changed staged FASTA membership")
    return root


def prepare(repo, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require new output and plan")
    source_path = repo / "benchmark_tools/results/three_kingdoms_sources_20260918.json"
    source = read_frozen(source_path, SOURCE_SHA)
    runtime_path = repo / "benchmark_tools/results/qfo_corrected_sonic_runtime_20260918.json"
    runtime = read_frozen(runtime_path, RUNTIME_SHA)
    validate_runtime(runtime)
    inputs = [row["files"]["staged"] for row in source["inputs"]]
    records = [record(source_path), record(runtime_path), record(ENTRYPOINT), runtime["python"],
               *[t["file"] for t in runtime["tools"].values()], *inputs]
    for path in ("benchmark_tools/normalize_three_kingdoms_orthogroups.py", "three_kingdoms/score_against_busco.py",
                 "three_kingdoms/busco/reference_orthogroups.txt", "benchmark_tools/inspect_sonicparanoid_runtime.py",
                 "benchmark_tools/prepare_qfo_corrected_sonic.py", "benchmark_tools/run_qfo_corrected_sonic.py",
                 "benchmark_tools/snapshot_runtime_trees.py", "benchmark_tools/audit_three_kingdoms_method_inputs.py",
                 "benchmark_tools/results/three_kingdoms_method_inputs_20260918.json"):
        records.append(record(repo / path))
    inventories = []
    for name, digest in INVENTORIES:
        path = repo / "benchmarks/work" / name
        read_frozen(path, digest)
        inventories.append(record(path))
    records.extend(inventories)
    for item in records:
        check(item)
    plan = {"status": "matched_three_kingdoms_sonic_frozen_unrun", "accuracy_admitted": False,
            "source": record(__file__), "input_source_audit": record(source_path), "inputs": inputs,
            "output_root": str(output), "native_argv": command(ENTRYPOINT, output), "reuse_search": False,
            "runtime": record(runtime_path), "runtime_inventories": inventories, "checked_records": records,
            "resources": {"node": "bizon", "cpus": 32, "memory_gib": 192, "hours": 72},
            "conversion": "Existing SonicParanoid ortholog_groups.tsv normalizer; group co-membership semantics.",
            "scoring": "Existing BUSCO-reference gene-pair micro P/R/F1, same 255-group reference; no retuning.",
            "remaining_gates": ["Independent terminal/native snapshot and output validation.",
                                "Frozen conversion and scoring execution with independent arithmetic audit.",
                                "Retain historical raw-input result; update comparative table only after admission."],
            "limitations": ["Contemporary matched-input run, not a pure causal test of stop-marker removal.",
                            "Current runtime identity does not prove historical April dependency versions.",
                            "Shared-host resource data are descriptive, not dedicated timing."]}
    validate_plan(plan)
    with destination.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return plan


def verify(plan_path, digest):
    plan = read_frozen(plan_path, digest)
    root = validate_plan(plan)
    if record(__file__)["sha256"] != plan["source"]["sha256"]:
        raise ValueError("Launcher source differs from frozen plan")
    for item in [plan["source"], *plan["checked_records"]]:
        check(item)
    runtime = read_frozen(Path(plan["runtime"]["path"]), RUNTIME_SHA)
    validate_runtime(runtime)
    package = verify_tree(runtime["package_inventory"])
    inventories = [verify_tree(read_frozen(Path(item["path"]), item["sha256"])) for item in plan["runtime_inventories"]]
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unreviewed loader preload")
    env = environment(plan, runtime)
    code = ("import importlib.util,json,sys; s=importlib.util.spec_from_file_location('probe',sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); s.loader.exec_module(m); print(json.dumps(m.inspect()))")
    done = subprocess.run([runtime["python"]["path"], "-B", "-s", "-c", code, runtime["source"]["path"]],
                          env=env, cwd=plan_path.parent, capture_output=True, text=True, check=True, timeout=240)
    compare_resolution(runtime, json.loads(done.stdout))
    if Path(env["PYTHONPYCACHEPREFIX"]).exists():
        raise ValueError("Isolated bytecode cache prefix must remain absent")
    return plan, root, env, {"package": package, "inventories": inventories}


def run(plan_path, digest, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require scheduled 32-CPU allocation")
    plan, root, env, before = verify(plan_path, digest)
    if root.exists():
        raise FileExistsError("No implicit restart or overwrite")
    if check_only:
        return {"status": "preflight_passed_no_inference", "runtime": before}
    root.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "plan": record(plan_path), "source": record(__file__),
              "job_id": os.environ["SLURM_JOB_ID"], "native_argv": plan["native_argv"], "environment": env,
              "node": os.uname().nodename, "runtime_before": before, "accuracy_admitted": False}
    status = root / "execution.json"
    try:
        directory = root / "input"
        directory.mkdir()
        copies = []
        for item in plan["inputs"]:
            check(item)
            target = directory / Path(item["path"]).name
            shutil.copy2(item["path"], target)
            copied = record(target)
            if any(copied[k] != item[k] for k in ("bytes", "sha256")):
                raise ValueError("Input copy differs")
            copies.append(copied)
        report.update(status="running", copied_inputs=copies, started_epoch=time.time())
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (root / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(root / "time.txt"), *plan["native_argv"]],
                                  cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
        report.update(exit_code=done.returncode, finished_epoch=time.time())
        _, _, _, report["runtime_after"] = verify(plan_path, digest)
        for item in copies:
            check(item)
        if {p.name for p in directory.glob("*.fasta")} != {Path(r["path"]).name for r in copies}:
            raise ValueError("Copied input membership changed")
        report.update(outputs=[record(p) for p in sorted((root / "output").rglob("*")) if p.is_file()],
                      native_log=record(root / "native.log"), timing=record(root / "time.txt"))
        done.check_returncode()
        report["status"] = "process_succeeded_pending_native_admission"
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--repo", type=Path, required=True)
    prep.add_argument("--output-root", type=Path, required=True)
    prep.add_argument("--plan", type=Path, required=True)
    execute = sub.add_parser("run")
    execute.add_argument("--plan", type=Path, required=True)
    execute.add_argument("--sha256", required=True)
    execute.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.output_root.resolve(), args.plan.resolve()) if args.action == "prepare"
              else run(args.plan.resolve(), args.sha256, args.check_only))
    print(json.dumps({"status": result["status"]}))
