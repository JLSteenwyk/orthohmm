"""Run native OrthoMCL in a fresh directory from admitted recovered BPO inputs."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import write_sources, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.recovered_orthomcl_inputs import verify_inputs
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_qfo_corrected_orthomcl import SOURCES_SHA, command, group_counts, pair_cache_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.stage_orthomcl_native_inputs import stage
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime


def preflight(root, admission, digest, job, executor, commit):
    runtime = verify_runtime(root)
    evidence = verify_inputs(root, admission, digest, job, executor, commit)
    results = root / "benchmark_tools/results"
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    sources = read_frozen(source_path, SOURCES_SHA)
    manifests = [read_frozen(results / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = [*evidence["checked_records"], record(source_path), *sources["originals"],
               *sources["checked_records"], record(__file__),
               *[record(p) for p in sorted(Path(__file__).parent.glob("*.py"))],
               record(Path(__file__).with_name("run_orthomcl_perl_script.pl")),
               record(Path(__file__).with_name("validate_orthomcl_bpo_indexes.pl"))]
    for item in checked:
        check(item)
    output = root / "benchmarks/results/qfo_blast_recovery_native_v1"
    if output.exists() or output.is_symlink():
        raise FileExistsError("Recovered native inference already attempted; no implicit resume")
    if shutil.disk_usage(output.parent).free < sum(r["bytes"] for r in evidence["inputs"].values()) + 10 * 1024**3:
        raise ValueError("Insufficient staging space plus 10-GiB reserve; graph output needs additional space")
    return evidence, manifests, checked, runtime, output


def run(root, admission, digest, job, executor, commit):
    job_id = os.environ.get("SLURM_JOB_ID", "")
    if (not job_id.isascii() or not job_id.isdigit() or int(job_id) < 1
            or os.environ.get("SLURM_CPUS_PER_TASK") != "180"
            or os.environ.get("SLURM_MEM_PER_NODE") != "921600"
            or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled 180-CPU 900-GiB native inference on bizon")
    evidence, manifests, checked, runtime, output = preflight(root, admission, digest, job, executor, commit)
    output.mkdir(exist_ok=False)
    tool, inputs = output / "native_tool", output / "inputs"
    env = {**environment(), "ORTHOMCL_PAIR_WORKERS": "64"}
    report = dict(status="recovered_native_staging", source=record(__file__), job_id=job_id,
                  node="bizon", admission=evidence["admission"], admission_job=evidence["admission_job"],
                  admission_scheduler=evidence["scheduler"], admission_accounting=evidence["accounting"],
                  checked_records=checked, runtime_before=runtime, environment=env,
                  query_coverage=evidence["query_coverage"], started_epoch=time.time(),
                  accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False)
    save_status(output / "report.json", report)
    try:
        sources = write_sources(tool, inputs, 180)
        checked.extend([*sources["originals"], *sources["configured_sources"]])
        report["native_sources"] = sources
        staged = stage(evidence["inputs"], inputs, 984137, 78, evidence["index_validation"])
        report["staging"] = record(inputs / "staging.json")
        mtimes = {k: Path(r["path"]).stat().st_mtime_ns for k, r in staged["staged"].items()}
        argv = command(tool, inputs)
        report.update(status="recovered_native_running", command=argv, cwd=str(output),
                      native_started_epoch=time.time())
        save_status(output / "report.json", report)
        with (output / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "native.time.txt"), *argv],
                                  cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report.update(native_exit_code=done.returncode, native_finished_epoch=time.time(),
                      native_log=record(output / "native.log"), native_timing=record(output / "native.time.txt"))
        if done.returncode:
            raise ValueError("Recovered native OrthoMCL inference failed")
        groups = list(tool.glob("*/all_orthomcl.out"))
        if len(groups) != 1:
            raise ValueError("Require one final native group file")
        counts = group_counts(groups[0], load_species(inputs / "all.gg"))
        for item in [*checked, *staged["staged"].values()]:
            check(item)
        if mtimes != {k: Path(r["path"]).stat().st_mtime_ns for k, r in staged["staged"].items()}:
            raise ValueError("Native inference rewrote staged inputs or indexes")
        for manifest in manifests:
            verify(manifest)
        caches = pair_cache_records(inputs)
        runtime_after = verify_runtime(root)
        report.update(status="recovered_native_exited_zero_pending_admission", content=counts,
                      runtime_after=runtime_after, staged_inputs=staged["staged"], input_mtimes_ns=mtimes,
                      pair_caches=caches, native_groups=record(groups[0]),
                      native_outputs=[record(p) for p in sorted(groups[0].parent.rglob("*")) if p.is_file()],
                      limitations=[
                          "Final native groups, not graph edges; separate terminal/output admission remains required.",
                          "Failed queries remain in provenance; no missing BLAST hits are repaired by inference.",
                          "Shared-host native-stage accounting is not matched end-to-end timing.",
                          "Existing pair-parallel patch with 64 workers; no scientific default changes."])
    except BaseException as error:
        report.update(status="recovered_native_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        save_status(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "executor"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.admission.resolve(), args.admission_sha256,
                 args.admission_job, args.executor.resolve(), args.commit)
    print(result["status"])
