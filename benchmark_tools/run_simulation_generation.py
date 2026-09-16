"""Execute one immutable panel generation task, preserving failed artifacts."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record


def verify_file(path, record):
    actual = file_record(path, path.parent)
    if actual["bytes"] != record["bytes"] or actual["sha256"] != record["sha256"]:
        raise ValueError(f"Artifact changed: {path}")


def child_path(base, relative):
    path = (base / relative).resolve()
    if not path.is_relative_to(base.resolve()):
        raise ValueError("Artifact path escapes expected root")
    return path


def preflight(manifest_path, expected_hash, panel, label):
    raw = manifest_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_hash:
        raise ValueError("Manifest checksum mismatch")
    manifest = json.loads(raw)
    if manifest["status"] != "materialized_not_executed" or len(manifest["simulation_runs"]) != 40 or len(manifest["datasets"]) != 70:
        raise ValueError("Not the complete frozen generation panel")
    matches = [r for r in manifest["simulation_runs"] if r["label"] == label]
    if len(matches) != 1:
        raise ValueError("Requested run not unique in manifest")
    run = matches[0]
    source = Path(manifest["source"])
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != manifest["source_commit"]:
        raise ValueError("Simulator source revision changed")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD"], check=True, capture_output=True)
    for record in manifest["native_sources"] + manifest["native_defaults"]:
        verify_file(child_path(source, record["path"]), record)
    for record in manifest["workflow_sources"]:
        verify_file(Path(record["absolute_path"]), record)
    for record in manifest.get("extra_inputs", []):
        verify_file(child_path(panel, record["path"]), record)
    for row in manifest["simulation_runs"]:
        for record in row["parameters"].values():
            verify_file(child_path(panel, record["file"]["path"]), record["file"])
    freeze = manifest["environment_freeze"]
    verify_file(child_path(panel, freeze["path"]), freeze)
    query = "import importlib.metadata as m,json,sys; names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); print(json.dumps({'python':sys.version,'inventory':''.join(f'{n}=={m.version(n)}\\n' for n in names)}))"
    current = json.loads(subprocess.check_output([manifest["python"], "-c", query], text=True))
    if current["python"] != manifest["python_version"] or current["inventory"] != child_path(panel, freeze["path"]).read_text():
        raise ValueError("Simulation interpreter/package inventory changed")
    native = child_path(panel, run["native_output"])
    if native.exists():
        raise FileExistsError("Native output already exists; no automatic restart")
    for command in run["commands"]:
        argv = command["argv"]
        if argv[0] != manifest["python"]:
            raise ValueError("Unexpected interpreter in stage command")
        target = child_path(panel, argv[argv.index("--output") + 1])
        if target.exists():
            raise FileExistsError(f"Stage output already exists: {target}")
    return manifest, run


def execute(run, environment, evidence, provenance=None):
    evidence.mkdir(parents=True, exist_ok=False)
    status = {"schema_version": 1, "label": run["label"], "status": "running",
              "started_epoch": time.time(), "stages": [], "method_inference_executed": False,
              "provenance": provenance or {}, "output_inventory_recorded": False}
    path = evidence / "status.json"
    try:
        for stage in run["commands"]:
            name = stage["stage"]
            if not name.isidentifier():
                raise ValueError("Invalid stage name")
            record = {"stage": name, "command": stage["argv"], "status": "running", "started_epoch": time.time()}
            status["stages"].append(record)
            path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
            command = ["/usr/bin/time", "-v", "-o", str(evidence / f"{name}.time.log"), *stage["argv"]]
            with (evidence / f"{name}.log").open("w") as log:
                completed = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False)
            record.update(exit_code=completed.returncode, finished_epoch=time.time(),
                          status="complete" if completed.returncode == 0 else "failed")
            if completed.returncode:
                raise RuntimeError(f"Generation stage {name} failed: {completed.returncode}")
        status["status"] = "complete"
    except Exception as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        status["finished_epoch"] = time.time()
        path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    manifest, run = preflight(args.manifest, args.manifest_sha256, args.panel.resolve(), args.run)
    if args.check_only:
        print(f"Verified generation preflight: {args.run}; no execution")
        return
    if args.evidence is None:
        parser.error("--evidence required unless --check-only")
    env = os.environ.copy()
    env.update(manifest["environment"])
    provenance = {"manifest": file_record(args.manifest.resolve(), args.manifest.resolve().parent),
                  "runner": file_record(Path(__file__).resolve(), Path(__file__).resolve().parent)}
    result = execute(run, env, args.evidence.resolve(), provenance)
    try:
        result["outputs"] = []
        targets = {Path(command["argv"][command["argv"].index("--output") + 1]) for command in run["commands"]}
        for target in sorted(targets):
            files = [p for p in sorted(target.rglob("*")) if p.is_file()]
            if not files:
                raise ValueError(f"Expected nonempty stage output: {target}")
            result["outputs"].extend(file_record(p, args.panel.resolve()) for p in files)
        result["output_inventory_recorded"] = True
    except Exception as error:
        result.update(status="output_inventory_failed", error=str(error))
        raise
    finally:
        (args.evidence / "status.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
