"""One unchanged frozen Leiden attempt on the failed cap400 final graph."""

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import timed_stage

PROTOCOL = "hillclimb_cap400_final_graph_recovery_20261016.json"
FAILURE = "benchmark_tools/results/cap400_final_graph_failure_20261016.json"
PARENT = "benchmarks/results/full_qfo_cap400_continuation_r1_v1_20261014"
REVISION = "2ccc411520d9e63678178a50207ba1fca87bfbb3"
PAYLOAD_FILES = ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")


def gene_table(path):
    names = path.read_text().splitlines()
    if (not names or len(names) > np.iinfo(np.int32).max
            or len(names) != len(set(names))
            or any(not n or any(c.isspace() for c in n) for n in names)):
        raise ValueError("Invalid complete gene table")
    return names


def reconstruct(abc, names_path, payload, count):
    if type(count) is not int or count < 1:
        raise ValueError("Require a positive exact edge count")
    names = gene_table(names_path)
    ids = {name: index for index, name in enumerate(names)}
    payload.mkdir(exist_ok=False)
    shutil.copyfile(names_path, payload / "gene_names.txt")
    arrays = [np.lib.format.open_memmap(payload / filename, mode="w+", dtype=dtype,
                                      shape=(count,))
              for filename, dtype in zip(PAYLOAD_FILES[1:], ("int32", "int32", "float64"))]
    observed = 0
    with abc.open() as stream:
        for index, line in enumerate(stream):
            parts = line.removesuffix("\n").split("\t")
            if index >= count or not line.endswith("\n") or len(parts) != 3:
                raise ValueError("Malformed edge row or excess edge count")
            a, b, token = parts
            weight = float(token)
            if a not in ids or b not in ids or not math.isfinite(weight) or str(weight) != token:
                raise ValueError("Unknown endpoint or noncanonical finite float64 token")
            arrays[0][index], arrays[1][index], arrays[2][index] = ids[a], ids[b], weight
            observed += 1
    if observed != count:
        raise ValueError("Retained graph has fewer rows than authenticated count")
    for array in arrays:
        array.flush()
    return {"genes": len(names), "edges": observed,
            "files": [record(payload / name) for name in PAYLOAD_FILES]}


def verify_rows(abc, names_path, payload, count):
    """Read every row back using inverse endpoint lookup and float bit patterns."""
    names = gene_table(names_path)
    if (payload / "gene_names.txt").read_bytes() != names_path.read_bytes():
        raise ValueError("Complete gene table or its order changed")
    arrays = [np.load(payload / name, mmap_mode="r", allow_pickle=False)
              for name in PAYLOAD_FILES[1:]]
    for array, dtype in zip(arrays, ("int32", "int32", "float64")):
        if array.shape != (count,) or array.dtype != np.dtype(dtype):
            raise ValueError("Reconstructed array shape or dtype differs")
    observed = 0
    with abc.open("rb") as stream:
        for index, raw in enumerate(stream):
            fields = raw.rstrip(b"\n").split(b"\t")
            if index >= count or len(fields) != 3 or not raw.endswith(b"\n"):
                raise ValueError("Independent graph row inventory differs")
            source, target = int(arrays[0][index]), int(arrays[1][index])
            if not (0 <= source < len(names) and 0 <= target < len(names)):
                raise ValueError("Independent graph endpoint outside complete gene table")
            weight = float(fields[2])
            if (names[source].encode() != fields[0] or names[target].encode() != fields[1]
                    or not math.isfinite(weight) or str(weight).encode() != fields[2]
                    or arrays[2][index].tobytes() != struct.pack("=d", weight)):
                raise ValueError("Independent graph endpoint/order/weight readback differs")
            observed += 1
    if observed != count:
        raise ValueError("Independent graph row count differs")
    return {"status": "all_saved_rows_and_float64_bits_verified", "genes": len(names),
            "edges": observed, "isolates_preserved_by_complete_gene_table": True,
            "edge_rows_deduplicated_or_reordered": False}


def require_failure(plan, parent, metrics, receipt, scheduler):
    if (plan.get("scientific_revision") != REVISION
            or parent.get("job_id") != "24245"
            or parent.get("status") != "failed_or_unverified_control_continuation"
            or parent.get("scores") != {} or parent.get("execution", {}).get("returncode") != 1
            or metrics.get("status") != "failed"
            or metrics.get("harness", {}).get("exit_code") != 1
            or metrics.get("harness", {}).get("git_commit") != REVISION
            or metrics.get("harness", {}).get("git_dirty") is not False
            or receipt.get("status") != "terminal_cap400_final_leiden_failure_and_retained_checkpoint_verified"
            or receipt.get("optimizer_attempted") is not False
            or receipt.get("scores") != {}):
        raise ValueError("Require authenticated original final-clustering failure")
    expected = {"JobIDRaw": "24245", "State": "FAILED", "ExitCode": "1:0",
                "Elapsed": "2-04:30:54", "NodeList": "bizon", "AllocCPUS": "32"}
    if scheduler != expected or receipt.get("scheduler") != expected:
        raise ValueError("Original scheduler state differs")
    case = plan["case"]["fresh_case"]
    if (parent["execution"]["command"] != case["command"]
            or metrics["harness"]["output_line_counts"]["edges"] != 33585180
            or receipt.get("final_graph_edges_reported_by_failed_harness") != 33585180):
        raise ValueError("Original command or final graph count differs")


def authenticate(repo, failure_pin):
    check(failure_pin)
    receipt = json.loads(Path(failure_pin["path"]).read_text())
    for pin in (receipt["execution_plan"], receipt["parent_result"], receipt["failed_native_metrics"],
                receipt["recovery_protocol"], *receipt["retained_native_products"]):
        check(pin)
    plan = json.loads(Path(receipt["execution_plan"]["path"]).read_text())
    parent = json.loads(Path(receipt["parent_result"]["path"]).read_text())
    metrics = json.loads(Path(receipt["failed_native_metrics"]["path"]).read_text())
    if (Path(plan["repo"]).resolve() != repo or Path(plan["output"]) != repo / PARENT
            or Path(receipt["execution_plan"]["path"]) != repo / PARENT / "execution_plan.json"):
        raise ValueError("Failure evidence belongs to a different run")
    accounting = subprocess.check_output(["sacct", "-j", "24245", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    if len(rows) != 1:
        raise ValueError("Parent scheduler record is not unique")
    require_failure(plan, parent, metrics, receipt, rows[0])
    for pin in [*plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]:
        check(pin)
    if runtime() != plan["runtime"]:
        raise ValueError("Scientific runtime changed")
    verify_input_inventory(plan["case"]["fresh_case"])
    science = Path(plan["scientific_source"])
    for item in metrics["harness"]["source_manifest"]:
        observed = record(science / item["path"])
        if (observed["bytes"], observed["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Frozen scientific source differs from failed harness")
    working = repo / PARENT / "cap400/native/output/orthohmm_working_res"
    abc = working / "orthohmm_edges.txt"
    expected = next(p for p in receipt["retained_native_products"] if p["path"] == str(abc))
    check(expected)
    checkpoint = working / "high_sensitivity_checkpoint"
    audited = audit(checkpoint, receipt["checkpoint_audit"]["manifest"]["sha256"])
    original_audit = receipt["checkpoint_audit"]
    check(original_audit["auditor"])
    if (any(audited[key] != value for key, value in original_audit.items() if key != "auditor")
            or any(audited["auditor"][key] != original_audit["auditor"][key]
                   for key in ("bytes", "sha256"))):
        raise ValueError("Original numeric checkpoint audit changed")
    if (audited["summary"]["genes"], audited["summary"]["species"]) != (984137, 78):
        raise ValueError("Wrong corrected QfO universe")
    return plan, receipt, abc, checkpoint / "gene_names.txt"


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require a frozen runner and new output namespace")
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    pins = [record(source / "benchmark_tools" / name) for name in (
        "probe_cap400_final_graph.py", "probe_cap400_final_graph.slurm", PROTOCOL)]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(source)
        blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(relative)])
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Runner differs from committed source")
    failure = record(repo / FAILURE)
    if record(source / FAILURE)["sha256"] != failure["sha256"]:
        raise ValueError("Failure receipt differs from frozen committed evidence")
    parent, receipt, abc, names = authenticate(repo, failure)
    if record(source / "benchmark_tools" / PROTOCOL)["sha256"] != receipt["recovery_protocol"]["sha256"]:
        raise ValueError("Recovery protocol differs from precommitted failure receipt")
    sources = [record(p) for p in sorted((source / "benchmark_tools").glob("*.py"))]
    sources += [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    result = {"status": "prepared_unrun_final_graph_probe", "repo": str(repo), "source_root": str(source),
        "source_revision": revision, "output": str(output), "failure": failure,
        "sources": sources + pins, "abc": record(abc), "gene_names": record(names),
        "scientific_source": parent["scientific_source"], "scientific_revision": REVISION,
        "runtime": parent["runtime"], "tool_path": parent["tool_path"],
        "edges": 33585180, "genes": 984137, "species": 78,
        "metadata": {"cpm_resolution": 0.1, "include_isolates": True, "seed": 4,
                     "output_directory": str(output / "native")},
        "scope": "One cached final-graph optimizer probe; no initial HMM rerun or accuracy admission."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(path):
    plan = json.loads(path.read_text())
    source, output, repo = (Path(plan[k]) for k in ("source_root", "output", "repo"))
    if (Path(__file__).resolve() != source / "benchmark_tools/probe_cap400_final_graph.py"
            or (output / "report.json").exists()):
        raise ValueError("Wrong runner or attempted namespace; no optimizer retries")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != (
            "1", "131072", "bizon", ""):
        raise ValueError("Require CPU-only one-CPU 128-GiB bizon allocation")
    report = {"status": "verifying_retained_final_graph", "execution_plan": record(path),
        "job_id": os.environ.get("SLURM_JOB_ID"), "optimizer_attempted": False,
        "scores": {}, "defaults_promoted": False, "goal_completed": False,
        "initial_hmm_search_rerun": False, "host_before": host_activity(), "scope": plan["scope"]}
    write_json(output / "report.json", report)
    pins = [record(path), plan["failure"], plan["abc"], plan["gene_names"], *plan["sources"]]
    try:
        for pin in pins:
            check(pin)
        parent, receipt, abc, names = authenticate(repo, plan["failure"])
        checkpoint_manifest = receipt["checkpoint_audit"]["manifest"]
        manifest = json.loads(Path(checkpoint_manifest["path"]).read_text())
        pins += [*parent["inputs"], *parent["sources"], *parent["runtime"]["files"],
                 *receipt["retained_native_products"], receipt["execution_plan"],
                 receipt["parent_result"], receipt["failed_native_metrics"], checkpoint_manifest]
        pins += [{"path": str(Path(checkpoint_manifest["path"]).parent / name), **identity}
                 for name, identity in manifest["files"].items()]
        if (parent["scientific_source"] != plan["scientific_source"] or runtime() != plan["runtime"]
                or plan["metadata"] != {"cpm_resolution": .1, "include_isolates": True,
                    "seed": 4, "output_directory": str(output / "native")}
                or (plan["edges"], plan["genes"], plan["species"]) != (33585180, 984137, 78)):
            raise ValueError("Prepared original scientific settings changed")
        payload = output / "payload"
        report["reconstruction"] = reconstruct(abc, names, payload, plan["edges"])
        report["independent_rows"] = verify_rows(abc, names, payload, plan["edges"])
        write_json(payload / "metadata.json", plan["metadata"])
        pins += report["reconstruction"]["files"] + [record(payload / "metadata.json")]
        for pin in pins:
            check(pin)
        (output / "native/orthohmm_working_res").mkdir(parents=True)
        report.update(status="original_frozen_final_optimizer_running", optimizer_attempted=True)
        write_json(output / "report.json", report)
        science = Path(plan["scientific_source"])
        env = {**os.environ, "PYTHONPATH": str(science), "PATH": plan["tool_path"],
            "PYTHONHASHSEED": "0", "CUDA_VISIBLE_DEVICES": "", "PYTHONFAULTHANDLER": "1",
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        command = [sys.executable, "-m", "orthohmm.leiden_worker", str(payload)]
        report["execution"] = timed_stage(command, output, "optimizer", science, env)
        if report["execution"]["returncode"]:
            raise RuntimeError("Original frozen final optimizer failed; no automatic retry")
        partition = output / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
        groups = read_partition(partition, set(gene_table(names)))
        report["partition"] = {"product": record(partition), "groups": len(groups), "genes": plan["genes"]}
        report["status"] = "final_graph_partition_verified_downstream_pending"
    except Exception as error:
        report.update(status="failed_final_graph_probe", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in pins:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during optimizer probe")
            report["postflight_pins_verified"] = True
        except Exception as error:
            report.update(status="failed_final_graph_probe", postflight_error=str(error))
        report["host_after"] = host_activity()
        write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, required=True, type=Path)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_unrun_final_graph_probe",
        "final_graph_partition_verified_downstream_pending"} else 1)
