"""Run fresh native OrthoMCL inference only from independently admitted BPO inputs."""

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
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.stage_orthomcl_native_inputs import stage
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_orthomcl
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species

ADMITTER = "914c3fe6f8618c3d90452e7a5b6d7f03f939ccd2"
SOURCES_SHA = "9ccd23429301d754f80fc72ae4ba4b28ef5036116e0670c5aebdb85a41b00ae6"


def admitted_inputs(admission, root):
    if (admission["status"] != "corrected_orthomcl_bpo_checkpoint_admitted"
            or admission["accuracy_admitted"] is not False or admission["publication_ready"] is not False
            or admission["validation"]["content"]["input_proteins"] != 984137):
        raise ValueError("Require independently admitted corrected BPO")
    scheduler = admission["scheduler"]
    expected = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"}
    if any(scheduler.get(key) != value for key, value in expected.items()):
        raise ValueError("Wrong BPO preparation accounting")
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    paths = {"bpo": base / "bpo_preparation/checkpoint/all.bpo",
             "offsets": base / "bpo_preparation/checkpoint/indexes/all_bpo.idx",
             "ranges": base / "bpo_preparation/checkpoint/indexes/all_bpo.se",
             "species": base / "work/all.gg"}
    if {r["path"] for r in admission["native_inputs"]} != {str(paths[k]) for k in ("bpo", "offsets", "ranges")} or len(admission["native_inputs"]) != 3:
        raise ValueError("Wrong admitted BPO/index inventory")
    result = {}
    for key, path in paths.items():
        records = [r for r in admission["checked_records"] if r["path"] == str(path)]
        if not records or any(r != records[0] for r in records) or records[0]["bytes"] <= 0:
            raise ValueError("Missing/conflicting admitted native input: " + key)
        if key != "species" and records[0] not in admission["native_inputs"]:
            raise ValueError("Native input record differs from admitted artifact")
        result[key] = records[0]
    return result


def group_counts(path, owners):
    seen, labels = set(), set()
    for label, genes in iter_orthomcl(path):
        members = set(genes)
        if not label or label in labels or not members or len(members) != len(genes) or not members <= owners.keys() or members & seen:
            raise ValueError("Malformed/unknown/repeated native group membership")
        seen.update(members)
        labels.add(label)
    if not labels:
        raise ValueError("Empty native groups")
    return {"groups": len(labels), "grouped_proteins": len(seen), "ungrouped_proteins": len(owners) - len(seen)}


def command(tool, inputs):
    return [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(tool),
            str(Path(__file__).with_name("run_orthomcl_perl_script.pl")), str(tool / "orthomcl.pl"),
            "--mode", "4", "--bpo_file", str(inputs / "all.bpo"), "--gg_file", str(inputs / "all.gg")]


def pair_cache_records(directory, species=78):
    expected = {f"pair_{i}.storable" for i in range(species * (species - 1) // 2)}
    paths = list(directory.glob("pair_*"))
    if {p.name for p in paths} != expected or any(not p.is_file() or p.is_symlink() or p.stat().st_size == 0 for p in paths):
        raise ValueError("Incomplete or unexpected native species-pair cache inventory")
    return [record(p) for p in sorted(paths)]


def preflight(root, admission_path, admission_sha, admission_job):
    python_runtime = verify_runtime(root)
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    if any(scheduler.get(k) != v for k, v in {"NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"}.items()):
        raise ValueError("Wrong BPO admission allocation")
    admission = read_frozen(admission_path, admission_sha)
    inputs = admitted_inputs(admission, root)
    executor = root / "benchmarks/work/publication_qfo_corrected_bpo_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
        raise ValueError("Changed BPO admission executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_bpo.py"):
        raise ValueError("Wrong BPO admission source")
    results = root / "benchmark_tools/results"
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    sources = read_frozen(source_path, SOURCES_SHA)
    tool = root / "benchmarks/results/qfo_corrected_orthomcl_v1/native_tool"
    if (sources["status"] != "corrected_orthomcl_native_sources_prepared_unrun"
            or sources["threads"] != 180 or sources["pair_workers_planned"] != 64
            or sources["pair_parallel_patch"] is not True or sources["tool_directory"] != str(tool)
            or {p.name for p in tool.iterdir()} != {"orthomcl.pl", "orthomcl_module.pm"}):
        raise ValueError("Changed native source configuration or existing native run")
    manifests = [read_frozen(results / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = [record(admission_path), admission["source"], *admission["checked_records"],
               *admission["validation"]["outputs"], *admission["native_inputs"], record(source_path),
               *sources["checked_records"], *sources["originals"], *sources["configured_sources"],
               record(__file__), record(Path(__file__).with_name("stage_orthomcl_native_inputs.py")),
               record(Path(__file__).with_name("run_orthomcl_perl_script.pl"))]
    for item in checked:
        check(item)
    base = tool.parent
    if (base / "inference_execution").exists():
        raise FileExistsError("Native inference already attempted; no implicit resume")
    if shutil.disk_usage(base).free < sum(r["bytes"] for r in inputs.values()) + 10 * 1024**3:
        raise ValueError("Insufficient staging space plus 10-GiB reserve; graph output needs additional space")
    return admission, inputs, tool, manifests, checked, python_runtime, scheduler, accounting


def run(root, admission_path, admission_sha, admission_job):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "180"
            or os.environ.get("SLURM_MEM_PER_NODE") != "921600" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled 180-CPU 900-GiB allocation on bizon")
    admission, inputs, tool, manifests, checked, runtime, scheduler, accounting = preflight(
        root, admission_path, admission_sha, admission_job)
    execution = tool.parent / "inference_execution"
    execution.mkdir(exist_ok=False)
    env = {**environment(), "ORTHOMCL_PAIR_WORKERS": "64"}
    report = {"status": "staging", "source": record(__file__), "checked_records": checked,
              "runtime_before": runtime, "admission_scheduler": scheduler, "admission_accounting": accounting,
              "job_id": os.environ["SLURM_JOB_ID"], "node": os.uname().nodename,
              "started_epoch": time.time(), "query_coverage": admission["query_coverage"],
              "environment": env, "accuracy_admitted": False, "publication_ready": False}
    try:
        staged = stage(inputs, execution / "inputs", 984137, 78, admission["validation"]["index_validation"])
        report["staging"] = record(execution / "inputs/staging.json")
        mtimes = {k: Path(r["path"]).stat().st_mtime_ns for k, r in staged["staged"].items()}
        argv = command(tool, execution / "inputs")
        report.update(status="running_native", command=argv, cwd=str(execution), native_started_epoch=time.time())
        (execution / "status.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (execution / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(execution / "native.time.txt"), *argv],
                                  cwd=execution, env=env, stdout=log, stderr=subprocess.STDOUT)
        report.update(native_exit_code=done.returncode, native_finished_epoch=time.time(),
                      native_log=record(execution / "native.log"), native_timing=record(execution / "native.time.txt"))
        if done.returncode:
            raise ValueError("Native OrthoMCL inference failed")
        native = list(tool.glob("*/all_orthomcl.out"))
        if len(native) != 1:
            raise ValueError("Require one final native group file")
        counts = group_counts(native[0], load_species(execution / "inputs/all.gg"))
        for item in [*checked, *staged["staged"].values()]:
            check(item)
        if mtimes != {k: Path(r["path"]).stat().st_mtime_ns for k, r in staged["staged"].items()}:
            raise ValueError("Native inference rewrote staged inputs or indexes")
        for manifest in manifests:
            verify(manifest)
        report.update(status="corrected_orthomcl_native_exited_zero_pending_admission", content=counts,
                      runtime_after=verify_runtime(root), staged_inputs=staged["staged"], input_mtimes_ns=mtimes,
                      pair_caches=pair_cache_records(execution / "inputs"),
                      native_groups=record(native[0]), native_outputs=[record(p) for p in sorted(native[0].parent.rglob("*")) if p.is_file()])
        report["limitations"] = ["Native final groups, not pre-clustering graph edges; independent terminal/output admission remains required.",
                                  "Shared-host native-stage time is not matched end-to-end timing.",
                                  "Small-fixture parallel parity does not prove all 64-worker schedules.",
                                  "Source failed-query diagnostics remain retained; conversion/inference do not repair missing search hits."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        (execution / "status.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", required=True, type=int)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only:
        preflight(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job)
        result = {"status": "preflight_passed_no_inference"}
    else:
        result = run(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job)
    print(json.dumps({"status": result["status"]}))
