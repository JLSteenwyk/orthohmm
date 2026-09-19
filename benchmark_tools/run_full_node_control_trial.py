"""One full-node diagnostic trial using the unchanged dual-bracket collector."""

import json
from pathlib import Path
import subprocess
import sys
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_dual_bracket_step import measure, evaluate
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.validate_full_node_control import validate_witnesses, common_intervals
from benchmark_tools.full_node_control_workload import require_start


def read(path):
    return json.loads(path.read_text())


def wait(path, stopped, seconds=30):
    deadline = time.monotonic() + seconds
    while not path.exists():
        if stopped.wait(.01):
            raise RuntimeError("Control start cancelled")
        if time.monotonic() >= deadline:
            raise TimeoutError(f"Missing control witness: {path.name}")
    return read(path)


def coordinate(directory, mode, stopped, outcome):
    process = None
    source = Path(__file__).with_name("full_node_control_workload.py")
    try:
        ready = wait(directory / "workload_ready.json", stopped)
        if mode == "contended":
            command = [sys.executable, "-B", str(source), "--directory", str(directory),
                       "--competitor", str(min(ready["cpus"]))]
            with (directory / "competitor.log").open("x") as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                wait(directory / "competitor_ready.json", stopped)
                save(directory / "workload_go.json", {"go": True})
                deadline = time.monotonic() + 30
                while process.poll() is None:
                    if stopped.wait(.02):
                        raise RuntimeError("Competitor cancelled")
                    if time.monotonic() >= deadline:
                        raise TimeoutError("Competitor exceeded deadline")
                outcome["competitor_exit_code"] = process.returncode
                if process.returncode != 0:
                    raise RuntimeError("Competitor failed")
        else:
            save(directory / "workload_go.json", {"go": True})
        outcome["status"] = "completed"
    except Exception as error:
        outcome.update(status="failed", error_type=type(error).__name__, error=str(error))
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
            outcome["competitor_exit_code"] = process.returncode


def trial(directory, mode, job):
    if mode not in {"steady", "churn", "contended"}:
        raise ValueError("Unknown control condition")
    directory.mkdir(exist_ok=False)
    work_dir = directory / "workload"
    work_dir.mkdir()
    measurement_dir = directory / "measurement"
    source = Path(__file__).with_name("full_node_control_workload.py")
    command = [sys.executable, "-B", str(source), "--directory", str(work_dir),
               "--worker", "churn" if mode == "churn" else "steady"]
    stopped, outcome = threading.Event(), {}
    controller = threading.Thread(target=coordinate, args=(work_dir, mode, stopped, outcome))
    controller.start()
    try:
        report = measure(command, measurement_dir, job, 20, 96 * 1024**3, 60, 1.)
        controller.join(timeout=35)
    finally:
        # Even a failed collector must not leave its competitor running.
        stopped.set()
        controller.join()
        save(directory / "controller.json", outcome)
    if outcome.get("status") != "completed":
        raise RuntimeError(f"Control controller failed: {outcome}")
    require_start(read(work_dir / "workload_go.json"))
    if evaluate(report["points"], report["native"], job) != report["screening"]:
        raise ValueError("Collector screening replay differs")
    competitor_ready = competitor = None
    if mode == "contended":
        competitor_ready = read(work_dir / "competitor_ready.json")
        competitor = read(work_dir / "competitor_done.json")
        if outcome.get("competitor_exit_code") != 0:
            raise ValueError("Competitor exit not validated")
    elif any(work_dir.glob("competitor*")):
        raise ValueError("Unexpected competitor evidence")
    validation = validate_witnesses(mode, read(work_dir / "workload_ready.json"),
        read(work_dir / "workload_done.json"), report["native"],
        report["points"][0]["native_membership"],
        report["points"][0]["host"][0]["raw"]["cgroup_membership"],
        competitor_ready, competitor)
    common = common_intervals(report["points"], validation)
    if not common:
        raise ValueError("No complete common-work observation intervals")
    screens = report["screening"]["narrow_intervals"]
    detected = any("excess_unassigned_cpu" in screens[i]["reasons"] for i in common)
    result = dict(mode=mode, job_id=job, status="workload_validated", validation=validation,
        common_intervals=common, common_narrow_flagged=[i for i in common if not screens[i]["screen_passed"]],
        positive_control_detected=detected if mode == "contended" else None,
        scientific_timings_admitted=False)
    save(directory / "trial.json", result)
    return result
