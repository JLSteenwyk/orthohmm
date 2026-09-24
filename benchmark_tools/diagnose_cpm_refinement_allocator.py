"""One allocator-debug refinement diagnostic; never admits scientific outputs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_failed_recovery_refinement import coverage, record
from benchmark_tools.run_blast_recovery_batch import save_status

PROTOCOL_SHA = "7ceae6cad1f6806529dead92ee0c978f28df2bdc87b5a379ead4028e84d27577"
READBACK_SHA = "aef33d04d1bec4dc469e7afc15c196e54174db89c7035650f8beb565942b9aa6"
RUNNER_SHA = "a277ab01e7fbcbfaa15f7a63a092c78183baaabbbcc25cd23640f6a750eabbe2"


def check(item):
    # Preserve the recorded path spelling; saved graph inputs include symlinks.
    if record(item["path"]) != item:
        raise ValueError("Diagnostic input identity changed: " + item["path"])


def environment(launcher):
    env = os.environ.copy()
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(key, None)
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
        PYTHONMALLOC="debug", PYTHONFAULTHANDLER="1")
    return env


def validate_child(child, reference, destination, names):
    if {k: v for k, v in child.items() if k != "output"} != {k: v for k, v in reference.items() if k != "output"}:
        raise ValueError("Diagnostic scientific metadata differs")
    output = record(destination)
    if child["output"] != output or any(output[key] != reference["output"][key] for key in ("bytes", "sha256")):
        raise ValueError("Diagnostic refinement bytes differ")
    return coverage(destination, names, reference["groups"])


def run(root, commit):
    if (os.environ.get("SLURM_CPUS_PER_TASK") != "1" or os.environ.get("SLURM_MEM_PER_NODE") != "65536"
            or os.environ.get("SLURMD_NODENAME") != "bizon" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require one-CPU 64-GiB bizon diagnostic allocation")
    executor = Path(__file__).resolve().parent.parent
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Diagnostic executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    output = root / "benchmarks/results/qfo_cpm_refinement_allocator_diagnostic_v1"
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    protocol = record(results / "QFO_RECOVERY_ADMISSION_FAILURE_22155.md")
    readback = record(results / "qfo_failed_refinement_readback_22155.json")
    native = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    runner = record(root / "benchmarks/work/cpm_checkpoint_recovery_v1_20260923/benchmark_tools/run_cpm_checkpoint_recovery.py")
    if (protocol["sha256"], readback["sha256"], runner["sha256"]) != (PROTOCOL_SHA, READBACK_SHA, RUNNER_SHA):
        raise ValueError("Diagnostic protocol or source changed")
    diagnostic = json.loads(Path(readback["path"]).read_text())
    native_report = json.loads((native / "status.json").read_text())
    reference = json.loads((native / "refinement.json").read_text())
    records = [record(__file__), protocol, readback, runner, *diagnostic["checked_records"],
        *native_report["checked_records"], *reference["modules"],
        *[record(native / "payload" / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")],
        record(native / "orthogroups_profiles.txt")]
    for item in records:
        check(item)
    from benchmark_tools.verify_qfo_replay_launcher import verify
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime = verify(core, launcher, runtime_path)
    if runtime != native_report["runtime_after"]:
        raise ValueError("Diagnostic runtime changed")
    output.mkdir()
    (output / "payload").symlink_to(native / "payload", target_is_directory=True)
    (output / "orthogroups_profiles.txt").symlink_to(native / "orthogroups_profiles.txt")
    command = [sys.executable, "-B", runner["path"], "--root", str(root), "--output", str(output), "--mode", "repeat-refinement"]
    report = dict(status="allocator_diagnostic_running", job_id=os.environ["SLURM_JOB_ID"],
        executor_commit=commit, command=command, checked_records=records, runtime_before=runtime,
        diagnostic_overrides=dict(PYTHONMALLOC="debug", PYTHONFAULTHANDLER="1"),
        child_attempts=0, seed_admitted=False, accuracy_evaluated=False, publication_ready=False)
    save_status(output / "status.json", report)
    try:
        report["child_attempts"] = 1
        save_status(output / "status.json", report)
        started = time.monotonic()
        with (output / "child.log").open("x") as log:
            child = subprocess.run(command, cwd=launcher, env=environment(launcher), stdout=log, stderr=subprocess.STDOUT)
        report.update(returncode=child.returncode, wall_s=time.monotonic() - started, child_log=record(output / "child.log"))
        if child.returncode:
            raise RuntimeError(f"Diagnostic child exited {child.returncode}; no retry")
        child_record = record(output / "refinement_repeat.json")
        child_report = json.loads(Path(child_record["path"]).read_text())
        report["coverage"] = validate_child(child_report, reference, output / "refinement_repeat.txt",
                                             (native / "payload/gene_names.txt").read_text().splitlines())
        report["child_report"] = child_record
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during diagnostic")
        for item in records:
            check(item)
        report.update(status="allocator_diagnostic_completed_not_admitted",
            limitations=["One diagnostic execution, not proof of memory safety or a repaired runtime.",
                         "No optimizer execution; original admission failure and blocked candidate remain unchanged."])
    except BaseException as error:
        report.update(status="allocator_diagnostic_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.commit)
