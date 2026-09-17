"""Execute a pinned biological application run; accuracy admission is separate."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.snapshot_orthohmm_input_order import record, snapshot

METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "sonicparanoid")
DIRECTORIES = ("high_sensitivity", "satellite_v2", "orthofinder", "sonicparanoid")


def pinned(row):
    path = Path(row["path"])
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise ValueError("Changed pinned file: " + str(path))
    return json.loads(data)


def select(spec, index):
    if spec.get("execution_authorized") is not True or spec.get("purpose") != "wgd_application":
        raise ValueError("Application execution is not authorized")
    plan = pinned(spec["command_plan"])
    if tuple(r["method"] for r in plan["runs"]) != METHODS or not 0 <= index < 4:
        raise ValueError("Unexpected method panel or index")
    root = Path(plan["output_root"]) / DIRECTORIES[index]
    if not root.is_absolute() or root.resolve() != root:
        raise ValueError("Output path is not absolute or traverses a symlink")
    return plan, plan["runs"][index], root


def check_inputs(plan, spec, runtime):
    inputs = pinned(plan["inputs"])
    order = pinned(spec["input_order"])
    directory = Path(inputs["inputs"][0]["path"]).parent
    if set(directory.iterdir()) != {Path(r["path"]) for r in inputs["inputs"]}:
        raise ValueError("Changed input directory membership")
    for row in [inputs["reference"], inputs["cohort"], inputs["protocol"], plan["contrast_protocol"]]:
        if record(row["path"]) != row:
            raise ValueError("Changed application reference or protocol")
    observed = snapshot(plan["core_root"], {"datasets": [{"proteomes": 4,
                        "input_directory": str(directory), "inputs": inputs["inputs"]}]}, runtime)
    if observed["datasets"] != order["datasets"] or observed["source"] != order["source"]:
        raise ValueError("Changed native input order or source")
    return inputs


def check_all(plan, spec):
    checks = []
    runtime = None
    for index, row in enumerate(spec["runtime_manifests"]):
        manifest = pinned(row)
        checks.append(verify(manifest))
        if index == 0:
            runtime = manifest
    if runtime is None:
        raise ValueError("Missing runtime inventory")
    inputs = check_inputs(plan, spec, runtime)
    for row in plan["entrypoints"]:
        if record(row["path"]) != row:
            raise ValueError("Changed native entrypoint")
    return inputs, checks


def check_copies(directory, originals):
    expected = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in originals}
    observed = {p.name: record(p) for p in directory.glob("*.fasta")}
    if {name: (row["bytes"], row["sha256"]) for name, row in observed.items()} != expected:
        raise ValueError("Changed comparator input copies")
    return list(observed.values())


def execute(argv, cwd, log, timeout):
    started = time.monotonic()
    with log.open("xb") as handle:
        proc = subprocess.Popen(argv, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT,
                                start_new_session=True)
        timed_out = False
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            # Only this command's new process group is terminated.
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            time.sleep(1)
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
    return {"exit_code": proc.returncode, "timed_out": timed_out,
            "elapsed_monotonic_s": time.monotonic() - started}


def run(spec, index, check_only=False):
    plan, selected, root = select(spec, index)
    if root.exists():
        raise FileExistsError(root)
    cache = root.with_name(root.name + "_pycache")
    if not sys.dont_write_bytecode or not sys.flags.no_user_site or sys.pycache_prefix != str(cache) or cache.exists():
        raise ValueError("Require disabled writes and an absent isolated bytecode prefix at startup")
    if os.environ.get("PYTHONHASHSEED") != "0":
        raise ValueError("Require fixed hash seed at startup")
    if not check_only and (platform.node() != "bizon" or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
                           or os.environ.get("SLURM_MEM_PER_NODE") != "131072"):
        raise ValueError("Require original host Slurm task with32CPUs/128GiB")
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unreviewed loader preload")
    for key in spec["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(spec["environment"])
    inputs, checks = check_all(plan, spec)
    if check_only:
        return {"status": "preflight_passed_no_inference", "runtime_checks": checks}
    root.mkdir(parents=True, exist_ok=False)
    report = {"method": selected, "job_id": os.environ["SLURM_JOB_ID"],
              "hostname": platform.node(), "environment": spec["environment"],
              "runtime_before": checks, "status": "prepared", "output_admitted": False,
              "limitations": ["Uncontrolled original-host runtime; not a scaling measurement.",
                              "Zero exit does not establish native output correctness."]}
    with (root / "preparation.json").open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    try:
        if "copy_inputs_to" in selected:
            destination = Path(selected["copy_inputs_to"])
            if destination != root / "input":
                raise ValueError("Unexpected comparator copy destination")
            destination.mkdir()
            for item in inputs["inputs"]:
                source = Path(item["path"])
                with source.open("rb") as src, (destination / source.name).open("xb") as dst:
                    shutil.copyfileobj(src, dst)
            report["copied_inputs"] = check_copies(destination, inputs["inputs"])
        # Freeze the immediate pre-inference enumeration again after preparation.
        check_inputs(plan, spec, pinned(spec["runtime_manifests"][0]))
        argv = ["/usr/bin/time", "-v", "-o", str(root / "native.time.log"), *selected["argv"]]
        report["executed_argv"] = argv
        report["cwd"] = str(root)
        report["native"] = execute(argv, root, root / "native.log", spec["timeout_s"])
        report["status"] = "native_exited_zero" if report["native"]["exit_code"] == 0 and not report["native"]["timed_out"] else "native_failed"
    except Exception as error:
        report.update(status="execution_error", error=repr(error))
    finally:
        try:
            _, report["runtime_after"] = check_all(plan, spec)
            if "copy_inputs_to" in selected:
                check_copies(Path(selected["copy_inputs_to"]), inputs["inputs"])
        except Exception as error:
            report.update(status="postflight_failed", postflight_error=repr(error))
        with (root / "execution.json").open("x") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    spec = pinned({"path": str(args.spec), "sha256": args.sha256})
    result = run(spec, args.index, args.check_only)
    print(json.dumps({"status": result["status"]}))
    raise SystemExit(0 if result["status"] in {"native_exited_zero", "preflight_passed_no_inference"} else 1)
