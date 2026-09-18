"""Run fresh corrected FastOMA inference after tree-bound input staging."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inspect_fastoma_runtime import docker_state
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_fastoma import ASSETS_SHA
from benchmark_tools.prepare_qfo_corrected_fastoma_assets import JAVA, IMAGE_DIGEST, image_identity
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree
from benchmark_tools.verify_ygob_validation import require_completed_job

STAGER = "82ef08e8a7d9511360ab52b62504871ba62753fc"
RUNTIME_SHA = "c194ac79775f17d1067bdb5fefdb406e74edab6a61c976f8203cdb04023dbb99"


def command(root, inputs, output, workflow):
    return ["/home/bizon/bin/nextflow", "-C", str(root / "benchmark_tools/fastoma_corrected_execution.config"),
            "run", str(workflow / "FastOMA.nf"), "-ansi-log", "false",
            "--input_folder", str(inputs), "--output_folder", str(output / "output"),
            "--omamer_db", str(root / "qfo_benchmark/data/LUCA.h5"),
            "--fasta_header_id_transformer", "UniProt", "--force_pairwise_ortholog_generation", "true",
            "--max_cpus", "180", "--max_memory", "700.GB", "-work-dir", str(output / "work")]


def environment():
    return {"HOME": "/home/bizon", "USER": "bizon", "LOGNAME": "bizon",
            "PATH": "/usr/bin:/bin:/usr/sbin:/sbin:/home/bizon/bin", "LANG": "C", "LC_ALL": "C",
            "JAVA_HOME": str(JAVA), "JAVA_CMD": str(JAVA / "bin/java"),
            "NXF_OFFLINE": "true", "NXF_VER": "22.10.8", "NXF_OPTS": "-Xmx1g",
            "NXF_HOME": "/home/bizon/.nextflow"}


def validate_stage(stage, root):
    if (stage["status"] != "corrected_fastoma_inputs_staged_pending_launch_freeze"
            or stage["execution_authorized"] is not False or stage["accuracy_evaluated"] is not False
            or stage["publication_ready"] is not False
            or stage["input_proteins"] != 984137 or stage["input_species"] != 78):
        raise ValueError("Require corrected, unrun FastOMA staging")
    directory = root / "benchmarks/work/qfo_corrected_fastoma_inputs_20260918"
    if stage["input_directory"] != str(directory) or len(stage["copies"]) != 79:
        raise ValueError("Changed staged directory or copy inventory")
    expected = set()
    for row in stage["copies"]:
        source, staged = row["source"], row["staged"]
        if source not in stage["checked_records"]:
            raise ValueError("Copy source not admitted")
        if any(source[k] != staged[k] for k in ("bytes", "sha256")):
            raise ValueError("Staged content differs from source")
        target = (directory / "proteome" / (Path(source["path"]).stem + ".fa")
                  if Path(source["path"]).suffix == ".fasta" else directory / "species_tree.nwk")
        if staged["path"] != str(target) or target in expected:
            raise ValueError("Unexpected or duplicate staged filename")
        expected.add(target)
    if directory / "species_tree.nwk" not in expected or len(expected) != 79:
        raise ValueError("Require species tree and 78 proteomes")
    observed = {p for p in directory.rglob("*") if p.is_file()}
    if observed != expected or any(p.is_symlink() for p in directory.rglob("*")):
        raise ValueError("Changed staged membership or unexpected symlink")
    return directory


def verify_runtime(runtime, env):
    verify_tree(runtime["runtime_trees"])
    for item in runtime["binaries"]:
        check(item)
    if docker_state() != runtime["docker"]:
        raise ValueError("Docker runtime changed")
    observed = subprocess.check_output(["/usr/bin/docker", "image", "inspect", IMAGE_DIGEST], env=env, text=True)
    if image_identity(json.loads(observed)) != runtime["image"]:
        raise ValueError("FastOMA image changed")


def preflight(root, stage_path, stage_sha, stage_job):
    accounting = subprocess.check_output(["sacct", "-j", str(stage_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, stage_job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Unexpected preparation allocation")
    stage = read_frozen(stage_path, stage_sha)
    executor = root / "benchmarks/work/publication_qfo_corrected_fastoma_prepare_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != STAGER:
        raise ValueError("Staging executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if stage["source"] != record(executor / "benchmark_tools/prepare_qfo_corrected_fastoma.py"):
        raise ValueError("Wrong staging source")
    inputs = validate_stage(stage, root)
    results = root / "benchmark_tools/results"
    assets_path = results / "qfo_corrected_fastoma_assets_20260918.json"
    assets = read_frozen(assets_path, ASSETS_SHA)
    runtime_path = results / "qfo_corrected_fastoma_runtime_20260918.json"
    runtime = read_frozen(runtime_path, RUNTIME_SHA)
    checked = [record(stage_path), stage["source"], record(assets_path), record(runtime_path),
               *stage["checked_records"], *[r["staged"] for r in stage["copies"]],
               *assets["sources"], assets["omamer_database"], record("/usr/bin/time")]
    for item in checked:
        check(item)
    env = environment()
    verify_runtime(runtime, env)
    output = root / "benchmarks/results/qfo_corrected_fastoma_v1"
    return {"stage": record(stage_path), "stage_scheduler": scheduler, "stage_accounting": accounting,
            "checked_records": checked, "runtime": record(runtime_path), "environment": env,
            "native_argv": command(root, inputs, output, Path(assets["workflow_root"])),
            "cwd": str(output / "run"), "output_root": str(output)}, runtime


def run(root, stage_path, stage_sha, stage_job, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "180"
            or os.environ.get("SLURM_MEM_PER_NODE") != "737280" or not os.environ.get("SLURM_JOB_ID")
            or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled bizon allocation: 180 CPUs and 720 GiB")
    if any(os.environ.get(k) for k in ("DOCKER_HOST", "DOCKER_CONTEXT")):
        raise ValueError("Docker connection overrides not allowed")
    report, runtime = preflight(root, stage_path, stage_sha, stage_job)
    output = Path(report["output_root"])
    if output.exists():
        raise FileExistsError("Existing FastOMA output; no implicit resume")
    if check_only:
        return {"status": "fastoma_preflight_passed_no_inference", **report}
    output.mkdir(parents=True, exist_ok=False)
    directory = Path(report["cwd"])
    directory.mkdir()
    status = output / "execution.json"
    report.update(status="preparing", source=record(__file__), job_id=os.environ["SLURM_JOB_ID"],
                  node=os.uname().nodename, accuracy_admitted=False, native_outputs_validated=False,
                  publication_ready=False, limitation="Supplied corrected OrthoFinder tree; shared-host inference is not matched timing.")
    try:
        report.update(status="running", started_epoch=time.time())
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (output / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *report["native_argv"]],
                                  cwd=directory, env=report["environment"], stdout=log, stderr=subprocess.STDOUT)
        report.update(exit_code=done.returncode, finished_epoch=time.time())
        after, _ = preflight(root, stage_path, stage_sha, stage_job)
        if any(after[k] != report[k] for k in after):
            raise ValueError("Preflight identity changed during inference")
        report.update(outputs=[record(p) for p in sorted((output / "output").rglob("*")) if p.is_file()],
                      log=record(output / "native.log"), timing=record(output / "time.txt"),
                      trace=record(directory / "trace.txt"), nextflow_log=record(directory / ".nextflow.log"))
        if done.returncode:
            raise RuntimeError(f"FastOMA process failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_native_admission"
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "stage"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--stage-sha256", required=True)
    parser.add_argument("--stage-job", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    run(args.root.resolve(), args.stage.resolve(), args.stage_sha256, args.stage_job, args.check_only)
