"""Execute the prespecified nine dual-bracket CPU controls on the DGX."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_dual_cpu_brackets import read_point, compare
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.run_native_pressure_controls import ORDER, validate_trial as validate_pressure_trial
from benchmark_tools.run_native_pressure_controls import summarize as summarize_pressure

PROTOCOL_SHA = "0ded70961fb6b4fd556a899498009b1ca3b511434824089d1c81862cc0729617"


def validate_trial(row, job):
    expected = compare(*row["points"], job, enforce_gap=False)
    if expected != row["brackets"]:
        raise ValueError("Dual-bracket replay differs")
    pressure_row = {**row, "points": [p["native_pressure"] for p in row["points"]],
                    "pressure": expected["native_pressure"]}
    witness = validate_pressure_trial(pressure_row, job)
    begin = row["points"][0]["host"][1]["finished_monotonic_ns"]
    end = row["points"][1]["host"][0]["started_monotonic_ns"]
    for work in (row["native"], row["competitor"]):
        if work is not None and not begin < work["started_ns"] < work["finished_ns"] <= end - 200000000:
            raise ValueError("Control work not enclosed with required post-work delay")
    return witness


def summarize(rows):
    pressure_rows = [{**r, "pressure": r["brackets"]["native_pressure"]} if r["status"] == "validated"
                     else r for r in rows]
    pressure = summarize_pressure(pressure_rows)
    checks = []
    for row in rows:
        result = None
        if row["status"] == "validated":
            screen = row["brackets"]["narrow"]
            result = (screen["reasons"] == ["excess_unassigned_cpu"] and not screen["screen_passed"]
                      if row["mode"] == "contended" else screen["screen_passed"] and not screen["reasons"])
        checks.append(dict(block=row["block"], mode=row["mode"], narrow_expectation_met=result))
    passed = pressure["all_response_checks_passed"] and all(c["narrow_expectation_met"] is True for c in checks)
    return dict(pressure=pressure, cpu_checks=checks, all_control_checks_passed=passed,
                scientific_timings_admitted=False, controlled_workload_verified=False)


def trial(directory, block, mode, job):
    directory.mkdir()
    original = set(os.sched_getaffinity(0))
    worker = Path(__file__).with_name("run_native_pressure_controls.py")
    argv = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
            sys.executable, "-B", str(worker), "--worker", str(directory)]
    competitor = None
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            candidates = sorted(original - {ready["cpu"]})
            if not candidates or ready["cpu"] not in original:
                raise ValueError("Cannot separate observer/native allocation affinities")
            observer = candidates[0]
            os.sched_setaffinity(0, {observer})
            before = read_point(ready["pid"], ready["membership"], job, directory / "failed_before.json")
            save(directory / "before.json", before)
            save(directory / "go.json", dict(mode=mode))
            if mode == "contended":
                competitor = subprocess.Popen([sys.executable, "-B", str(worker), "--burn-on", str(ready["cpu"])],
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            native = wait_file(directory / "done.json")
            load = None
            if competitor is not None:
                stdout, stderr = competitor.communicate(timeout=15)
                save(directory / "competitor.json", dict(stdout=stdout, stderr=stderr, exit_code=competitor.returncode))
                if competitor.returncode:
                    raise ValueError("Competitor exited nonzero")
                load = json.loads(stdout)
            time.sleep(.2)
            after = read_point(ready["pid"], ready["membership"], job, directory / "failed_after.json")
            save(directory / "after.json", after)
        finally:
            try:
                if competitor is not None and competitor.poll() is None:
                    competitor.kill()
                    competitor.wait()
                for name, value in (("go.json", dict(mode="quiet")), ("release.json", dict(release=True))):
                    if not (directory / name).exists():
                        save(directory / name, value)
                code = process.wait(timeout=60)
            finally:
                os.sched_setaffinity(0, original)
    row = dict(block=block, mode=mode, ready=ready, observer_allowed=sorted(original), observer_cpu=observer,
               native=native, competitor=load, points=[before, after], worker_exit_code=code, argv=argv,
               brackets=compare(before, after, job, enforce_gap=False))
    row["validation"] = validate_trial(row, job)
    row["status"] = "validated"
    return row


def run(output):
    if (platform.node() != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "256"):
        raise ValueError("Require DGX two-CPU/256MiB control allocation")
    base = Path(__file__).resolve().parent
    protocol = base / "results/DUAL_BRACKET_CONTROL_PROTOCOL_20260919.md"
    sources = {str(p.relative_to(base)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [*sorted(base.glob("*.py")), protocol]}
    if sources[str(protocol.relative_to(base))] != PROTOCOL_SHA:
        raise ValueError("Frozen protocol changed")
    output.mkdir(exist_ok=False)
    job, rows = int(os.environ["SLURM_JOB_ID"]), []
    for block, modes in enumerate(ORDER):
        for mode in modes:
            directory = output / f"trial_{len(rows):02d}"
            try:
                row = trial(directory, block, mode, job)
            except Exception as error:
                row = dict(block=block, mode=mode, status="failed", error_type=type(error).__name__, error=str(error))
                directory.mkdir(exist_ok=True)
            save(directory / "trial.json", row)
            rows.append(row)
    for name, digest in sources.items():
        if hashlib.sha256((base / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Recipe source changed during controls")
    result = dict(status="dual_bracket_controls_completed", job_id=job, sources=sources, trials=rows,
                  summary=summarize(rows), python=sys.version, executable=sys.executable,
                  kernel=platform.release(), host=platform.node(), publication_ready=False,
                  limitations=["Single-core CPU controls only; not full-node native-tool calibration.",
                               "All failed trials retained; no timing admission, correction, or selective repeats."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output.resolve())
    raise SystemExit(0 if result["summary"]["all_control_checks_passed"] else 1)
