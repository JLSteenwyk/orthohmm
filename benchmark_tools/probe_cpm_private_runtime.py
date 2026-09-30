"""One unchanged-runner refinement control in the reconstructed private runtime."""

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import resource
import signal
import subprocess
import time

from benchmark_tools.probe_cpm_partition_parser import check, record
from benchmark_tools.audit_failed_recovery_refinement import coverage

RESTORED = "benchmarks/work/publication_restored_archive_ob_20260929"
NATIVE = "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
RUNNER = "benchmarks/work/cpm_checkpoint_recovery_v1_20260923/benchmark_tools/run_cpm_checkpoint_recovery.py"
RUNNER_SHA = "a277ab01e7fbcbfaa15f7a63a092c78183baaabbbcc25cd23640f6a750eabbe2"
PROTOCOL = "benchmark_tools/results/QFO_CPM_PRIVATE_RUNTIME_PROTOCOL_20260930.md"
HELPER_PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_ENVIRONMENT_PROTOCOL_20260930.md"
PINS = {
    "benchmark_tools/results/restored_archive_full_ob_result_22377.json":
        "43e882d1795996c567b6348d33d547052b708910a61a8bd5cf7da0ce059ba4da",
    RESTORED + "/plan.json": "e80067f7a5947bc357831e0d65c7d766d09899f004cfc6a31d1274a8dfe15448",
    RESTORED + "/independent_admission/package_audits.json":
        "639cf9ea0a42f997374e040687c642733c641a057d9bea4688ab725e825af4e0",
    RESTORED + "/run/inference/pyvenv.cfg":
        "3491a0965a1ca9f94d3cf0fc66e0550387946144ef108bad0a1286e7d571b779",
    "benchmark_tools/results/qfo_cpm_refinement_gc_20260928.json":
        "5087f465598b0234402791bc24187b193b18ddb1607493aeb4104e2230f15fa8",
    NATIVE + "/refinement.json": "12c0f5066768b829504b1266069b69975acff9961c416a0df2b77e9699ef818a",
}


def package_records(audit, site):
    records, seen = [], set()
    for package in audit["packages"]:
        if package["matched_files"] != len(package["matched"]):
            raise ValueError("Incomplete retained package inventory")
        for row in package["matched"]:
            member = PurePosixPath(row["member"])
            if (not member.parts or member.is_absolute() or ".." in member.parts
                    or "\\" in row["member"] or row["member"] in seen):
                raise ValueError("Unsafe or duplicate package member")
            seen.add(row["member"])
            path = site / str(member)
            if path.is_symlink() or not path.resolve().is_relative_to(site.resolve()):
                raise ValueError("Indirect private package member")
            records.append(dict(path=str(path.resolve()), bytes=row["bytes"], sha256=row["sha256"]))
    return records


def bound_base_records(plan, base):
    base = base.resolve()
    records = []
    for row in plan["checked_records"]:
        path = Path(row["path"])
        if not path.is_relative_to(base):
            continue
        canonical = path.resolve()
        if not canonical.is_relative_to(base):
            raise ValueError("Restored base member escapes its runtime")
        records.append(dict(row, path=str(canonical)))
    return unique_records(records)


def environment(launcher):
    env = dict(os.environ)
    removed = ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")
    for key in removed:
        env.pop(key, None)
    overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
        PYTHONMALLOC="debug", PYTHONFAULTHANDLER="1", OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.update(overrides)
    return env, overrides, list(removed)


def child_limits():
    resource.setrlimit(resource.RLIMIT_AS, (64 * 1024**3, 64 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (300, 300))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})


def execute(command, cwd, env, stdout, stderr, timeout):
    with stdout.open("xb") as out, stderr.open("xb") as err:
        child = subprocess.Popen(command, cwd=cwd, env=env, stdout=out, stderr=err,
                                 start_new_session=True, preexec_fn=child_limits)
        try:
            return child.wait(timeout=timeout), False
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
            return child.returncode, True


def validate_result(result, reference, destination):
    if ({key: value for key, value in result.items() if key != "output"}
            != {key: value for key, value in reference.items() if key != "output"}):
        raise ValueError("Private-runtime child metadata differs")
    actual = record(destination)
    if (any(result["output"][key] != actual[key] for key in ("bytes", "sha256"))
            or any(actual[key] != reference["output"][key] for key in ("bytes", "sha256"))):
        raise ValueError("Private-runtime partition differs")
    return actual


def unique_records(records):
    indexed = {}
    for row in records:
        previous = indexed.setdefault(row["path"], row)
        if previous != row:
            raise ValueError("Conflicting evidence identities")
    return list(indexed.values())


def prepared_runtime(path, expected_sha):
    identity = record(path)
    if identity["sha256"] != expected_sha:
        raise ValueError("Changed helper environment receipt")
    prepared = json.loads(path.read_bytes())
    if (prepared["status"] != "helper_complete_private_environment_prepared"
            or prepared["refinement_attempts"] != 0
            or any(prepared[key] is not False for key in
                   ("seed_admitted", "accuracy_evaluated", "publication_ready"))):
        raise ValueError("Helper environment not independently prepared")
    prefix = Path(prepared["prefix"])
    if prepared["interpreter"] != record(prefix / "bin/python"):
        raise ValueError("Changed prepared interpreter")
    check([*prepared["checked_records"], prepared["import_report"], prepared["base_python"]])
    return prepared, identity


def _run(root, output, protocol_sha, prepared_path=None, prepared_sha=None):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    pins = [record(root / name) for name in PINS]
    if [row["sha256"] for row in pins] != list(PINS.values()):
        raise ValueError("Changed prerequisite evidence")
    protocol = record(root / (HELPER_PROTOCOL if prepared_path is not None else PROTOCOL))
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed private-runtime protocol")
    restored = json.loads((root / next(iter(PINS))).read_bytes())
    if (restored.get("job_id") != 22377 or restored.get("reproduction_equal") is not True
            or restored.get("status") != "restored_archive_orthobench_independently_audited"):
        raise ValueError("Require independently admitted private runtime")
    plan = json.loads((root / RESTORED / "plan.json").read_bytes())
    package_audit = json.loads((root / RESTORED / "independent_admission/package_audits.json").read_bytes())
    prior = json.loads((root / "benchmark_tools/results/qfo_cpm_refinement_gc_20260928.json").read_bytes())
    reference = json.loads((root / NATIVE / "refinement.json").read_bytes())
    prepared, prepared_identity = (prepared_runtime(prepared_path, prepared_sha)
                                  if prepared_path is not None else (None, None))
    prefix = Path(prepared["prefix"]) if prepared else root / RESTORED / "run/inference"
    python = prefix / "bin/python"
    restored_python = root / RESTORED / "run/inference/bin/python"
    base = restored_python.resolve().parent.parent
    base_records = bound_base_records(plan, base)
    if not base_records or not any(row["path"] == str(restored_python.resolve()) for row in base_records):
        raise ValueError("Private interpreter is not bound to restored base")
    if prepared and (Path(prepared["base_prefix"]).resolve() != base
                     or prepared["base_python"] != record(restored_python)):
        raise ValueError("Prepared helper base differs from restored base")
    records = [*pins, protocol, record(__file__),
               record(root / "benchmark_tools/probe_cpm_partition_parser.py"),
               record(root / "benchmark_tools/audit_failed_recovery_refinement.py"),
               *prior["checked_records"],
               *reference["modules"], *base_records]
    if prepared:
        records.extend([prepared_identity, *prepared["checked_records"], prepared["import_report"]])
    else:
        records.extend(package_records(package_audit["inference"], prefix / "lib/python3.10/site-packages"))
    runner = record(root / RUNNER)
    if runner["sha256"] != RUNNER_SHA:
        raise ValueError("Changed frozen refinement runner")
    records.append(runner)
    records = unique_records(records)
    check(records)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    env, overrides, removed = environment(launcher)
    output.mkdir(parents=True)
    (output / "payload").symlink_to(root / NATIVE / "payload", target_is_directory=True)
    (output / "orthogroups_profiles.txt").symlink_to(root / NATIVE / "orthogroups_profiles.txt")
    # The original runner and its fresh set(names) readback execute unchanged.
    command = [str(python), "-B", runner["path"], "--root", str(root),
               "--output", str(output), "--mode", "repeat-refinement"]
    report = dict(status="private_runtime_control_started", attempts=1, refinement_attempts=0, command=command,
        cwd=str(launcher), environment_overrides=overrides, removed_environment_keys=removed,
        checked_records=records, seed_admitted=False, accuracy_evaluated=False, publication_ready=False,
        limits=dict(cpu_seconds=300, wall_seconds=360, address_space_bytes=64 * 1024**3, affinity_cpus=1),
        limitations=["One reconstructed-runtime observation, not a repair or causal runtime attribution.",
            "No forced GC, debugger, optimizer, scoring, automatic retry or dependency release.",
            "Private package/base pins omit generated metadata/bytecode and relocated non-site wheel data.",
            "Original failed admission and all historical runtime/source records remain unchanged."])
    if prepared:
        report["helper_environment"] = prepared_identity
    started = time.monotonic()
    try:
        probe_code = ("import gc,json,sys,numpy;print(json.dumps(dict(version=sys.version,"
            "executable=sys.executable,numpy=numpy.__version__,numpy_path=numpy.__file__,"
            "gc_enabled=gc.isenabled(),gc_thresholds=list(gc.get_threshold()))))")
        code, timeout = execute([str(python), "-I", "-B", "-c", probe_code], launcher, env,
                               output / "runtime.stdout", output / "runtime.stderr", 60)
        if code or timeout:
            raise RuntimeError("Private runtime probe failed; refinement not attempted")
        observed = json.loads((output / "runtime.stdout").read_bytes())
        if (not observed["version"].startswith("3.10.13") or observed["numpy"] != "2.2.6"
                or not observed["gc_enabled"] or observed["gc_thresholds"] != [700, 10, 10]
                or not Path(observed["numpy_path"]).is_relative_to(prefix)):
            raise ValueError("Private runtime identity or GC defaults differ")
        report["runtime_probe"] = observed
        report["refinement_attempts"] = 1
        code, timeout = execute(command, launcher, env, output / "stdout", output / "stderr", 360)
        report.update(returncode=code, timed_out=timeout)
        if code or timeout:
            raise RuntimeError("Private refinement failed; no retry")
        result_path = output / "refinement_repeat.json"
        result = json.loads(result_path.read_bytes())
        report["partition"] = validate_result(result, reference, output / "refinement_repeat.txt")
        names = (root / NATIVE / "payload/gene_names.txt").read_text().splitlines()
        report["coverage"] = coverage(output / "refinement_repeat.txt", names, reference["groups"])
        report["child_report"] = record(result_path)
        report["result"] = result
        check(records)
        report["status"] = "private_runtime_refinement_control_completed_not_admitted"
    except BaseException as error:
        report.update(status="private_runtime_control_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["wall_seconds_descriptive_only"] = time.monotonic() - started
        report["logs"] = [record(path) for path in sorted(output.glob("*stdout")) + sorted(output.glob("*stderr"))]
        with (output / "report.json").open("x") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return report


def run(root, output, protocol_sha, prepared_path=None, prepared_sha=None):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    try:
        if (prepared_path is None) != (prepared_sha is None):
            raise ValueError("Require both helper receipt and its SHA256")
        return _run(root, output, protocol_sha, prepared_path, prepared_sha)
    except BaseException as error:
        output.mkdir(parents=True, exist_ok=True)
        if not (output / "report.json").exists():
            with (output / "report.json").open("x") as handle:
                json.dump(dict(status="private_runtime_preflight_failed", refinement_attempts=0,
                    error_type=type(error).__name__, error=str(error), seed_admitted=False,
                    accuracy_evaluated=False, publication_ready=False), handle, indent=2, sort_keys=True)
                handle.write("\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--prepared-environment", type=Path)
    parser.add_argument("--prepared-sha256")
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute(), args.protocol_sha256,
        args.prepared_environment, args.prepared_sha256)
