"""Execute frozen simulation inference; native conversion and scoring stay separate."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_generation import child_path, verify_file
from benchmark_tools.verify_simulation_histories import compare_histories, verify_generation

METHOD_HASH = "4f0717fa8b0c34d5a6982a3193772d64bc39a3d75b8fe23195d9ce1b21eee186"
GENERATION_HASH = "ee31ea38d3b5c047abf80f04636f06838959f941d6704649a216bb93958343b2"


def read_frozen(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"Frozen manifest changed: {path}")
    return json.loads(raw)


def execution_environment(manifest):
    env = os.environ.copy()
    env.update(manifest["environment_overrides"])
    env["PATH"] = os.pathsep.join([*manifest["prepend_path"], env.get("PATH", "")])
    resolved = {}
    for name in ("mafft", "FastTree", "diamond"):
        found = shutil.which(name, path=env["PATH"])
        expected = manifest["tool_entrypoints"][name]
        if found is None or Path(found).resolve() != Path(expected["absolute_path"]).resolve():
            raise ValueError(f"Unexpected executable resolution: {name}")
        verify_file(Path(found), expected)
        resolved[name] = found
    return env, resolved


def verify_environment(manifest):
    for record in (manifest["core_sources"] + manifest["adapter_sources"] +
                   manifest["orthofinder_distribution"] + list(manifest["tool_entrypoints"].values())):
        verify_file(Path(record["absolute_path"]), record)
    commit = subprocess.check_output(["git", "-C", manifest["core_root"], "rev-parse", "HEAD"], text=True).strip()
    if commit != manifest["core_commit"]:
        raise ValueError("Frozen core revision changed")
    query = "import importlib.metadata as m,json,sys; names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in names}}))"
    interpreters = {"orthohmm": manifest["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
                    "orthofinder": str(Path(manifest["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python")}
    for name, executable in interpreters.items():
        actual = json.loads(subprocess.check_output([executable, "-c", query], text=True))
        if actual != manifest["environments"][name]:
            raise ValueError(f"Package inventory changed: {name}")


def verify_inputs(dataset, generation, panel, generation_hash=GENERATION_HASH):
    runs = {r["label"]: r for r in generation["simulation_runs"]}
    parent = dataset["parent"]
    native = verify_generation(panel, runs[parent], generation_hash)
    check = next(c for c in generation["history_equivalence_checks"] if parent in (c["first"], c["second"]))
    other = check["second"] if parent == check["first"] else check["first"]
    paired = verify_generation(panel, runs[other], generation_hash)
    if not compare_histories(native, paired)["matched"]:
        raise ValueError("Required matched biological histories differ")
    status = json.loads((panel / "execution" / parent / "status.json").read_text())
    inventory = {r["path"]: r for r in status["outputs"]}
    inputs, truth = Path(dataset["input"]), Path(dataset["truth"])
    if not truth.exists():
        derived = panel / "derived" / str(dataset["seed"]) / "manifest.json"
        key = derived.relative_to(panel).as_posix()
        if key not in inventory:
            raise ValueError("Missing truth without verified derived-condition evidence")
        condition = json.loads(derived.read_text())["conditions"][dataset["condition"]]
        if condition["status"] != "inapplicable":
            raise ValueError("Missing truth for applicable condition")
        return {"status": "inapplicable", "reason": condition, "inputs": []}
    paths = [truth, *sorted(inputs.glob("*"))]
    prefix = inputs.resolve().relative_to(panel.resolve()).as_posix() + "/"
    expected_inputs = {name for name in inventory if name.startswith(prefix)}
    actual_inputs = {p.resolve().relative_to(panel.resolve()).as_posix() for p in paths[1:]}
    if actual_inputs != expected_inputs:
        raise ValueError("Generated species input file set changed")
    if len(paths) < 3 or any(not p.is_file() for p in paths):
        raise ValueError("Expected truth and at least two species input files")
    records = []
    for path in paths:
        key = path.resolve().relative_to(panel.resolve()).as_posix()
        if key not in inventory:
            raise ValueError(f"Unrecorded inference input: {path}")
        verify_file(path, inventory[key])
        records.append(dict(inventory[key], absolute_path=str(path)))
    return {"status": "ready", "truth": records[0], "inputs": records[1:],
            "matched_history_pair": check}


def copy_inputs(source, destination, expected):
    destination.mkdir(parents=True, exist_ok=False)
    for record in expected:
        path = Path(record["absolute_path"])
        if path.parent.resolve() != source.resolve():
            raise ValueError("Input copy source mismatch")
        verify_file(path, record)
        target = destination / path.name
        shutil.copy2(path, target)
        verify_file(target, record)


def execute(dataset, order, env, evidence, verified, provenance):
    if evidence.exists() or any(Path(m["output"]).exists() or ("metrics" in m and Path(m["metrics"]).exists())
                                for m in dataset["methods"].values()):
        raise FileExistsError("Existing inference artifacts; no automatic restart")
    evidence.mkdir(parents=True)
    status = {"schema_version": 1, "dataset": dataset["label"], "status": "running",
              "started_epoch": time.time(), "verified_inputs": verified, "provenance": provenance,
              "methods": {}, "accuracy_evaluated": False, "native_outputs_validated": False}
    path = evidence / "status.json"

    def save():
        temporary = evidence / "status.tmp"
        temporary.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)

    if verified["status"] == "inapplicable":
        status.update(status="inapplicable", finished_epoch=time.time())
        save()
        return status
    save()
    for name in order:
        method = dataset["methods"][name]
        record = {"status": "running", "argv": method["argv"], "started_epoch": time.time()}
        status["methods"][name] = record
        save()
        try:
            for item in verified["inputs"]:
                verify_file(Path(item["absolute_path"]), item)
            if "copy_inputs_to" in method:
                copy_inputs(Path(method["copy_inputs_from"]), Path(method["copy_inputs_to"]), verified["inputs"])
            else:
                Path(method["output"]).parent.mkdir(parents=True, exist_ok=True)
            command = ["/usr/bin/time", "-v", "-o", str(evidence / f"{name}.time.log"), *method["argv"]]
            with (evidence / f"{name}.log").open("w") as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
            record.update(exit_code=result.returncode, status="process_succeeded" if result.returncode == 0 else "failed")
            if result.returncode == 0:
                files = [p for p in sorted(Path(method["output"]).rglob("*")) if p.is_file()]
                if "metrics" in method:
                    files.append(Path(method["metrics"]))
                if not files:
                    raise ValueError("Process exited successfully without output files")
                record["outputs"] = [dict(file_record(p, p.parent), absolute_path=str(p)) for p in files]
        except Exception as error:
            record.update(status="failed", error=str(error))
        record["finished_epoch"] = time.time()
        save()
    status.update(status="finished_pending_native_validation", finished_epoch=time.time(),
                  failed_methods=[n for n, r in status["methods"].items() if r["status"] == "failed"])
    save()
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", default=METHOD_HASH)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    manifest = read_frozen(args.manifest, args.manifest_sha256)
    generation_hash = manifest["generation_manifest"]["sha256"]
    generation = read_frozen(Path(manifest["generation_manifest"]["absolute_path"]), generation_hash)
    matches = [d for d in manifest["datasets"] if d["label"] == args.dataset]
    if len(matches) != 1:
        raise ValueError("Unknown or duplicate dataset")
    dataset = matches[0]
    evidence = args.evidence if args.evidence is not None else Path(dataset["methods"]["orthohmm_high_sensitivity"]["output"]).parent / "execution"
    verify_environment(manifest)
    env, resolved = execution_environment(manifest)
    verified = verify_inputs(dataset, generation, args.panel.resolve(), generation_hash)
    if args.check_only:
        print(f"Verified {args.dataset}: {verified['status']}; no inference")
        return
    provenance = {"method_manifest_sha256": args.manifest_sha256, "generation_manifest_sha256": generation_hash,
                  "path_resolved_executables": resolved, "resolution_scope": "PATH lookup, not process execution tracing",
                  "sources": [file_record(Path(__file__).with_name(n), Path(__file__).parent) for n in
                              ("run_simulation_methods.py", "verify_simulation_histories.py", "run_simulation_generation.py", "benchmark_production.py")],
                  "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                  "timing_limitation": "Shared-machine inference timings are not controlled scaling measurements"}
    result = execute(dataset, manifest["execution_order"], env, evidence.resolve(), verified, provenance)
    if result.get("failed_methods"):
        raise SystemExit("One or more methods failed; other method outcomes preserved")


if __name__ == "__main__":
    main()
