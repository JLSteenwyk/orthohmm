"""Explicit recovered high-CPM handoff to frozen inferred phylogeny, unscored."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

JOB = "22386"
COMMIT = "c3f1968b913a83ba47e7fb021347bdd83e2f78fc"
EXECUTOR = "benchmarks/work/qfo_cpm_helper_candidate_admission_executor_20261001"
ADMISSION = "benchmarks/work/qfo_cpm_helper_candidate_admission_20260930/status.json"
ADMISSION_SHA = "55c696a037d13b48d47dc752fcb6b69b5b8bba6a7b990acfb08b0d8fbf6e07c0"
SOURCE_SHA = "cc9745c1f1cb8ddf4c595c6cc85b1cd7fec8fce1d68e8c1bd6e4cd657f68786f"
READBACK = "benchmark_tools/results/qfo_cpm_helper_candidate_admission_readback_20261001.json"
READBACK_SHA = "8a5e832febbccf7d4d55a5d2b16e4454a54f83fc3c5570d9c16dfed6b7f12e3a"
PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_PHYLOGENY_PROTOCOL_20261001.md"
OUTPUT = "benchmarks/results/qfo_cpm_helper_phylogeny_v1/cpm_high"
PYTHON = "/home/bizon/anaconda3/bin/python"


def pinned(root, relative, sha):
    item = record(root / relative)
    if item["sha256"] != sha:
        raise ValueError("Recovered phylogeny input identity changed: " + relative)
    return json.loads(Path(item["path"]).read_bytes()), item


def completed(accounting):
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == JOB]
    if len(rows) != 1 or tuple(rows[0][key] for key in
            ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
                JOB, "COMPLETED", "0:0", "2", "64G", "bizon"):
        raise ValueError("Require completed recovered candidate admission " + JOB)
    return rows[0]


def select_arm(root, report):
    if (report["status"] != "cpm_helper_recovered_candidates_admitted_unscored"
            or report["candidate_admitted"] is not True or report["arm"] != "cpm_high"
            or type(report["index"]) is not int or report["index"] != 1
            or report["seed_handoff"] != "explicit_helper_runtime_seed_amendment"
            or any(report[key] is not False for key in
                   ("accuracy_evaluated", "downstream_admitted", "publication_ready"))
            or report["verification"] != dict(genes=984137, seed_groups=390845,
                candidate_groups=346866, reconstructed_merges=43979)):
        raise ValueError("Require exact unscored recovered high-CPM admission")
    arm = report["candidate_arm"]
    directory = root / "benchmarks/results/qfo_cpm_helper_recovered_candidates_v1/candidate/orthohmm_working_res"
    expected = {"candidate_partition": directory / "orthohmm_edges_clustered.txt",
                "membership_constraints": directory / "phylogeny_candidate_merges.json",
                "seed_partition": root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1/orthogroups_profiles_refined.txt"}
    if (arm["candidate_expansion"] is not True
            or any(arm[key]["path"] != str(path) for key, path in expected.items())
            or arm["expansion"]["profile"] != "satellite_v2"
            or arm["expansion"]["membership_policy"] != "high_confidence_pair"):
        raise ValueError("Wrong recovered CPM seed/candidate/constraint handoff")
    return {"label": "cpm_high", "partition": arm["candidate_partition"],
            "constraints": arm["membership_constraints"], "seed_partition": arm["seed_partition"]}


def verify_candidates(root):
    accounting = subprocess.check_output(["sacct", "-j", JOB, "-X", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting)
    readback, readback_record = pinned(root, READBACK, READBACK_SHA)
    report, admission_record = pinned(root, ADMISSION, ADMISSION_SHA)
    arm = select_arm(root, report)
    if (readback["status"] != "recovered_candidate_admission_independently_read_back_unscored"
            or readback["admission"] != admission_record or readback["candidate_arm"] != report["candidate_arm"]
            or readback["verification"] != report["verification"] or readback["candidate_admitted"] is not True
            or any(readback[key] is not False for key in
                   ("accuracy_evaluated", "downstream_admitted", "publication_ready"))
            or readback["source_revision"] != COMMIT):
        raise ValueError("Recovered candidate readback disagrees")
    executor = root / EXECUTOR
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Recovered candidate admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    source = record(executor / "benchmark_tools/admit_helper_cpm_candidates.py")
    if source["sha256"] != SOURCE_SHA or source != report["source"]:
        raise ValueError("Wrong recovered candidate admission source")
    records = [admission_record, readback_record, readback["submission"], readback["log"], readback["time"],
               *report["checked_records"], arm["partition"], arm["constraints"], arm["seed_partition"]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting recovered candidate provenance")
        unique[item["path"]] = item
        check(item)
    for item in readback["source_git_bindings"]:
        path = Path(item["path"])
        if item["git_revision"] != COMMIT or executor not in path.parents:
            raise ValueError("Wrong recovered admission Git binding")
        actual = record(path)
        if actual != {key: item[key] for key in ("path", "bytes", "sha256")}:
            raise ValueError("Recovered admission Git-bound file changed")
        blob = subprocess.check_output(["git", "-C", str(executor), "show", f"{COMMIT}:{path.relative_to(executor)}"])
        if blob != path.read_bytes():
            raise ValueError("Recovered admission source differs from frozen Git blob")
    return {"arm": arm, "admission": report, "admission_record": admission_record,
            "admission_executor": str(executor), "scheduler": scheduler, "accounting": accounting,
            "readback_record": readback_record, "checked_records": list(unique.values())}


def verify_sources(root):
    verified = verify_candidates(root)
    from benchmark_tools.run_qfo_parameter_phylogeny import verify_baseline

    manifest, original, launcher, prepared, environment, records = verify_baseline(root)
    return {**verified, "manifest": manifest, "original": original, "launcher": str(launcher),
            "prepared": str(prepared), "environment": environment,
            "checked_records": [*verified["checked_records"], *records]}


def run(root, protocol_sha):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_MEM_PER_NODE") != "196608" or os.environ.get("SLURM_ARRAY_TASK_ID")
            or Path(sys.executable).resolve() != Path(PYTHON).resolve()):
        raise ValueError("Require standalone 32-CPU/192-GiB bizon task in the frozen full runtime")
    output = root / OUTPUT
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed recovered phylogeny protocol")
    verified = verify_sources(root)
    from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
    from benchmark_tools.run_qfo_factorial_cell import native_command
    from benchmark_tools.run_simulation_methods import read_frozen, execution_environment, execute

    launcher, prepared = Path(verified["launcher"]), Path(verified["prepared"])
    cell = variant_cell(verified["original"], verified["arm"], output)
    argv, equivalence = native_command(cell, launcher, prepared)
    env, resolved = execution_environment(verified["environment"])
    env.update(verified["manifest"]["environment_overrides"], PYTHONPATH=str(launcher))
    helpers = [record(path) for path in sorted(Path(__file__).parent.glob("*.py"))]
    provenance = {"source": record(__file__), "protocol": protocol, "helpers": helpers,
        "verified": verified, "cell": cell, "executed_argv": argv,
        "launcher_source_equivalence": equivalence, "resolved_tools": resolved,
        "cwd": str(launcher), "job_id": os.environ["SLURM_JOB_ID"],
        "seed_handoff": "explicit_helper_runtime_seed_amendment",
        "scope": "Unscored recovered high-CPM candidates; inferred phylogeny with validated raw-tree checkpoint reuse; incremental shared-host execution"}
    output.mkdir(parents=True, exist_ok=False)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    postflight = {"status": "running", "accuracy_evaluated": False,
                  "native_outputs_validated": False, "publication_ready": False}
    cwd = Path.cwd()
    try:
        fresh = output / "fresh_candidate_admission.json"
        admission = verified["admission"]
        command = [sys.executable, "-B", str(Path(verified["admission_executor"]) / "benchmark_tools/admit_helper_cpm_candidates.py"),
                   "--root", str(root), "--preparation-sha256", admission["preparation"]["sha256"],
                   "--protocol-sha256", admission["protocol"]["sha256"], "--output", str(fresh)]
        postflight["admission_command"] = command
        with (output / "admission.log").open("x") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        fresh_record = record(fresh)
        if read_frozen(fresh, fresh_record["sha256"]) != admission:
            raise ValueError("Fresh recovered candidate admission disagrees")
        postflight["fresh_candidate_admission"] = fresh_record
        method = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
        inputs = {"status": "ready", "inputs": [{**item, "absolute_path": item["path"]}
                  for item in verified["manifest"]["input_fastas"]]}
        os.chdir(launcher)
        execution = execute({"label": cell["label"], "methods": {cell["label"]: method}}, [cell["label"]],
                            env, output / "execution", inputs, provenance)
        os.chdir(cwd)
        if execution["failed_methods"]:
            raise RuntimeError("Recovered CPM phylogeny failed; preserve outputs without retry")
        if verify_sources(root) != verified:
            raise ValueError("Recovered phylogeny inputs/runtime changed during execution")
        for item in [provenance["source"], protocol, fresh_record, *helpers, *verified["checked_records"]]:
            check(item)
        for pair in equivalence:
            check(pair["prepared"])
            check(pair["executed"])
        postflight.update(status="complete_pending_native_validation", cell=cell)
    except BaseException as error:
        postflight.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        os.chdir(cwd)
        (output / "postflight.json").write_text(json.dumps(postflight, indent=2, sort_keys=True) + "\n")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.protocol_sha256)
