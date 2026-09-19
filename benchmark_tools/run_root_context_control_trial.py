"""One bounded root-context trial; all work and service cleanup remain witnessed."""

import json
from pathlib import Path
import subprocess
import sys
import threading
import time

from benchmark_tools.measure_native_root_context import measure
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_full_node_control_trial import wait
from benchmark_tools.root_context_control_design import MODES, unit_name, check_service_scope, service_command
from benchmark_tools.validate_full_node_control import validate_witnesses, integer


def read(path):
    return json.loads(path.read_text())


def coordinate(directory, mode, job, index, manager, stopped, outcome):
    process = None
    unit = unit_name(job, index)
    scope = None
    try:
        ready = wait(directory / "workload_ready.json", stopped)
        if mode == "user-contended":
            source = Path(__file__).with_name("full_node_control_workload.py")
            command = service_command(sys.executable, source, directory, unit, min(ready["cpus"]))
            outcome.update(command=command, unit=unit, manager=manager)
            with (directory / "service.log").open("x") as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                competitor = wait(directory / "competitor_ready.json", stopped)
                scope = check_service_scope(competitor["membership"], manager, unit)
                outcome["service_scope"] = scope
                save(directory / "workload_go.json", {"go": True})
                deadline = time.monotonic() + 50
                while process.poll() is None:
                    if stopped.wait(.02):
                        raise RuntimeError("Root-context service cancelled")
                    if time.monotonic() >= deadline:
                        raise TimeoutError("Owned user service exceeded deadline")
                outcome["service_exit_code"] = process.returncode
                if process.returncode != 0:
                    raise RuntimeError("Owned user service failed")
            observations, deadline = [], time.monotonic() + 15
            while True:
                exists = (Path("/sys/fs/cgroup") / scope.lstrip("/")).exists()
                observations.append(dict(monotonic_ns=time.monotonic_ns(), exists=exists))
                if not exists or time.monotonic() >= deadline:
                    break
                time.sleep(.1)
            outcome["removal"] = observations
            if exists:
                raise RuntimeError("Owned service scope did not disappear")
        else:
            save(directory / "workload_go.json", {"go": True})
        outcome["status"] = "completed"
    except Exception as error:
        outcome.update(status="failed", error_type=type(error).__name__, error=str(error))
    finally:
        if process is not None and process.poll() is None:
            try:
                result = subprocess.run(["systemctl", "--user", "stop", unit],
                                        capture_output=True, text=True, timeout=15)
                outcome["owned_service_stop"] = dict(returncode=result.returncode,
                    stdout=result.stdout, stderr=result.stderr)
                if result.returncode != 0:
                    outcome["cleanup_error"] = "Owned service stop failed"
            except Exception as error:
                outcome["cleanup_error"] = str(error)
            try:
                process.wait(timeout=50)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
                outcome["cleanup_error"] = "Owned service wait exceeded runtime allowance"
            outcome["service_exit_code"] = process.returncode


def validate_trial(directory, mode, job, index, manager, measured, context):
    if mode not in MODES:
        raise ValueError("Unknown root-context condition")
    native = measured["native"]
    if (native["exit_code"] != 0 or native["timed_out"] is not False
            or not all(integer(native[k]) for k in ("started_ns", "finished_ns"))):
        raise ValueError("Native control command failed")
    if mode == "idle":
        if native["finished_ns"] - native["started_ns"] < 20_000_000_000:
            raise ValueError("Idle command shorter than frozen duration")
        validation = dict(common_started_ns=native["started_ns"], common_finished_ns=native["finished_ns"],
                          scientific_timings_admitted=False)
        if any((directory / "workload").iterdir()):
            raise ValueError("Idle condition contains unexpected workload evidence")
    else:
        work = directory / "workload"
        controller = read(directory / "controller.json")
        if controller.get("status") != "completed" or "cleanup_error" in controller:
            raise ValueError("Workload coordinator failed")
        kwargs = {}
        if mode == "user-contended":
            before, after = read(work / "competitor_ready.json"), read(work / "competitor_done.json")
            scope = check_service_scope(after["membership"], manager, unit_name(job, index))
            if (controller.get("service_exit_code") != 0 or controller.get("service_scope") != scope
                    or not controller.get("removal") or controller["removal"][-1]["exists"] is not False):
                raise ValueError("Owned-service exit or removal not validated")
            kwargs = dict(competitor_ready=before, competitor=after,
                          expected_competitor_membership=after["membership"])
        elif any(work.glob("competitor*")) or (work / "service.log").exists():
            raise ValueError("Unexpected user-service evidence")
        validation = validate_witnesses("contended" if mode == "user-contended" else mode,
            read(work / "workload_ready.json"), read(work / "workload_done.json"), native,
            measured["points"][0]["native_membership"],
            measured["points"][0]["host"][0]["raw"]["cgroup_membership"], **kwargs)
    start, finish = validation["common_started_ns"], validation["common_finished_ns"]
    points = measured["points"]
    if not (points[0]["root_context"]["host_before"]["started_ns"] <= start
            < finish <= points[-1]["root_context"]["host_after"]["finished_ns"]):
        raise ValueError("Context observations do not enclose common workload")
    common = [i for i, (a, b) in enumerate(zip(points, points[1:]))
              if start <= a["host"][0]["started_monotonic_ns"]
              and b["root_context"]["host_after"]["finished_ns"] <= finish]
    if not common:
        raise ValueError("No jointly enclosed common-work intervals")
    user_cpu = context["context"]["observation_window"]["scope_cpu_usec"]["/user.slice"] / 1e6
    return dict(status="root_context_workload_validated", mode=mode, index=index, job_id=job,
        validation=validation, common_intervals=common, enclosing_user_cpu_s=user_cpu,
        positive_control_response=user_cpu >= 5 if mode == "user-contended" else None,
        scientific_timings_admitted=False, environmental_validity_established=False)


def trial(directory, mode, job, index, manager):
    if mode not in MODES:
        raise ValueError("Unknown root-context condition")
    directory.mkdir(exist_ok=False)
    work = directory / "workload"
    work.mkdir()
    source = Path(__file__).with_name("full_node_control_workload.py")
    command = ["/usr/bin/sleep", "20"] if mode == "idle" else [sys.executable, "-B", str(source),
        "--directory", str(work), "--worker", "churn" if mode == "churn" else "steady"]
    stopped, outcome = threading.Event(), {}
    controller = None
    if mode != "idle":
        controller = threading.Thread(target=coordinate,
            args=(work, mode, job, index, manager, stopped, outcome))
        controller.start()
    try:
        measured, context = measure(command, directory / "measurement", job)
        if controller is not None:
            controller.join(timeout=65)
    finally:
        stopped.set()
        if controller is not None:
            controller.join()
        save(directory / "controller.json", outcome)
    result = validate_trial(directory, mode, job, index, manager, measured, context)
    save(directory / "trial.json", result)
    return result
