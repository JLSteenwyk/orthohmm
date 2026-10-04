"""Rebind exactly three known packaging changes; never launch/retry inference."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.snapshot_runtime_trees import inventory
from benchmark_tools.check_threadripper_runtime import check_manifests
from benchmark_tools.inspect_native_python_lookup import compare_lookup

OLD_SHA = "7e3ac02b46cf1ead426ec19889fc5bd23e6a0e01dbd3737bf374c3ae3180ed1d"
CHANGES = {
    "bundle_publication_handoff.py": ("9be5629233687af860856d297a698b7df0705d44e208c911f1937d0687f7ef6e", "0b45067395262d168733c772e8f2a40ab6b069112a319eb031cbcc0e6b044a00"),
    "bundle_publication_review.py": ("f102f5f70e2c11d78e3ddfa0e0778c21372eebc4c889bb6347154ae341b31de2", "4b9646bf7fd88c750bf9df4e7fad22c4fac6d5f9810cf211a639c15009a07999"),
    "bundle_publication_source.py": ("b199863d32a2370d8ce0b20dfbeb2e082d3486a0fc89a0dafd1b568b2bb0cc32", "cf5f2c49ffe8871de6ec5dbf471d16eaa949cb338218ff1f18b49e3e571f3880")}


def delta(before, after):
    a = {r["path"]: r for r in before["records"]}
    b = {r["path"]: r for r in after["records"]}
    wanted = {str(ROOT / "benchmark_tools" / name): pair for name, pair in CHANGES.items()}
    changed = {p for p in a.keys() & b.keys() if a[p] != b[p]}
    if (set(a) != set(b) or changed != set(wanted) or before.keys() != after.keys()
            or any(before[k] != after[k] for k in before if k != "records")):
        raise ValueError("Runtime changed beyond exactly three packaging helpers")
    rows = []
    for path, (old_sha, new_sha) in sorted(wanted.items()):
        if (a[path]["sha256"] != old_sha or b[path]["sha256"] != new_sha
                or any(a[path][k] != b[path][k] for k in ("path", "mode", "kind"))):
            raise ValueError("Unexpected packaging helper identity")
        rows.append(dict(before=a[path], after=b[path]))
    return rows


def refresh(work, output):
    if work.exists() or output.exists():
        raise FileExistsError("Require fresh lookup validation paths")
    old_path = ROOT / "benchmark_tools/results/threadripper_private_lookup_preparation_sync_20261004.json"
    old_ref = record(old_path)
    if old_ref["sha256"] != OLD_SHA:
        raise ValueError("Previous lookup receipt differs")
    def read(ref):
        check(ref)
        data = json.loads(Path(ref["path"]).read_text())
        check(ref)
        return data
    old = read(old_ref)
    prior = read(old["binding"])
    first_path, first_sha = prior["runtime_specs"][0]
    first_ref = record(first_path)
    if first_ref["sha256"] != first_sha:
        raise ValueError("Historical runtime manifest differs")
    before = read(first_ref)
    work.mkdir(parents=True, exist_ok=False)
    print("Checking exact packaging-only runtime differences", flush=True)
    after = inventory(before["roots"])
    changes = delta(before, after)
    save(work / "native_os_helpers.json", after)
    binding = dict(prior)
    current_ref = record(work / "native_os_helpers.json")
    binding["runtime_specs"] = [[current_ref["path"], current_ref["sha256"]], *prior["runtime_specs"][1:]]
    binding["runtime_manifests"] = [current_ref, *prior["runtime_manifests"][1:]]
    binding["supersedes"] = old["binding"]
    save(work / "binding.json", binding)
    binding_ref = record(work / "binding.json")
    save(work / "helper_delta.json", dict(previous=first_ref, current=current_ref, changed=changes,
        added=[], removed=[], scientific_private_os_drift=False, source=record(__file__)))
    inspector = ROOT / "benchmark_tools/inspect_native_python_lookup.py"
    if record(inspector) != old["source"]:
        raise ValueError("Native lookup inspector changed")
    checks = {}
    started = time.monotonic()
    checks["before"] = check_manifests(binding["runtime_specs"])
    checks["before_wall_s"] = time.monotonic() - started
    command = [binding["controller_python"]["path"], "-B", str(inspector),
        "--baseline", old["baseline"]["path"], "--baseline-sha256", old["baseline"]["sha256"],
        "--binding", binding_ref["path"], "--binding-sha256", binding_ref["sha256"], "--output", str(work / "lookup_current")]
    env = os.environ.copy()
    for key in ("PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        env.pop(key, None)
    env.update(PYTHONNOUSERSITE="1", PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1",
               PYTHONPYCACHEPREFIX=str(work / "absent_controller_cache"),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (work / "lookup_process.log").open("x") as log:
        process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, timeout=240)
    save(work / "lookup_process.json", dict(command=command, exit_code=process.returncode))
    if process.returncode or (work / "absent_controller_cache").exists():
        raise ValueError("Native lookup process failed or wrote bytecode")
    interpreters = {}
    for name in ("orthohmm", "orthofinder"):
        prior_ref = old["interpreters"][name]["reports"][-1]
        observed_ref = record(work / "lookup_current" / (name + ".json"))
        observed = read(observed_ref)
        interpreters[name] = dict(reports=[prior_ref, observed_ref], modules=len(observed["modules"]),
                                  comparison=compare_lookup(read(prior_ref), observed))
    started = time.monotonic()
    checks["after"] = check_manifests(binding["runtime_specs"])
    checks["after_wall_s"] = time.monotonic() - started
    report = dict(old, binding=binding_ref, interpreters=interpreters, supersedes=old_ref,
        current_validation=dict(source=record(__file__), runtime_checks=checks,
            delta=record(work / "helper_delta.json"), process=record(work / "lookup_process.json"), log=record(work / "lookup_process.log")))
    save(work / "lookup.json", report)
    save(output, report)
    print(json.dumps(dict(status=report["status"], lookup=record(output), changed_helpers=len(changes)), indent=2), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    refresh(args.work.absolute(), args.output.absolute())
