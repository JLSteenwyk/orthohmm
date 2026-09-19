"""Replay one archived full-node control without admitting scientific timing."""

from pathlib import Path

from benchmark_tools.replay_dual_native_measurement import replay as replay_measurement
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_full_node_control_trial import read
from benchmark_tools.full_node_control_workload import require_start
from benchmark_tools.validate_full_node_control import validate_witnesses, common_intervals


def inventory(directory):
    paths = sorted(directory.rglob("*"))
    if any(path.is_symlink() for path in paths):
        raise ValueError("Symlinks are not permitted in control evidence")
    return [record(path) for path in paths if path.is_file()]


def replay(directory, remote_directory, recipe_directory, python, mode, job):
    directory, remote_directory, recipe_directory = map(Path, (directory, remote_directory, recipe_directory))
    if mode not in {"steady", "churn", "contended"}:
        raise ValueError("Unknown control condition")
    if not remote_directory.is_absolute() or not recipe_directory.is_absolute() or not Path(python).is_absolute():
        raise ValueError("Expected runtime and deployed paths must be absolute")
    evidence = inventory(directory)
    work = directory / "workload"
    command = [python, "-B", str(recipe_directory / "benchmark_tools/full_node_control_workload.py"),
               "--directory", str(remote_directory / "workload"),
               "--worker", "churn" if mode == "churn" else "steady"]
    measured = replay_measurement(directory / "measurement", job, command, expected_timeout_s=60)
    report = measured["measured"]
    controller = read(directory / "controller.json")
    if controller.get("status") != "completed":
        raise ValueError("Controller did not complete")
    require_start(read(work / "workload_go.json"))
    ready, done = read(work / "workload_ready.json"), read(work / "workload_done.json")
    expected_names = {f"{kind}_{cpu}.json" for cpu in done["cpus"] for kind in ("ready", "done")}
    actual_names = {path.name for pattern in ("ready_*.json", "done_*.json") for path in work.glob(pattern)}
    if expected_names != actual_names or list(work.glob("failed_*.json")):
        raise ValueError("Worker witness inventory differs or contains failures")
    for before, after in zip(ready["workers"], done["workers"]):
        cpu = after["cpu"]
        if before != read(work / f"ready_{cpu}.json") or after != read(work / f"done_{cpu}.json"):
            raise ValueError("Embedded worker witnesses disagree with raw files")
    competitor_ready = competitor = None
    if mode == "contended":
        if type(controller.get("competitor_exit_code")) is not int or controller["competitor_exit_code"] != 0:
            raise ValueError("Competitor exit status invalid")
        competitor_ready, competitor = read(work / "competitor_ready.json"), read(work / "competitor_done.json")
    elif list(work.glob("competitor*")) or "competitor_exit_code" in controller:
        raise ValueError("Unexpected competitor")
    validation = validate_witnesses(mode, ready, done, report["native"],
        report["points"][0]["native_membership"],
        report["points"][0]["host"][0]["raw"]["cgroup_membership"], competitor_ready, competitor)
    common = common_intervals(report["points"], validation)
    if not common:
        raise ValueError("No enclosed common-work interval")
    screens = report["screening"]["narrow_intervals"]
    expected = dict(mode=mode, job_id=job, status="workload_validated", validation=validation,
        common_intervals=common, common_narrow_flagged=[i for i in common if not screens[i]["screen_passed"]],
        positive_control_detected=any("excess_unassigned_cpu" in screens[i]["reasons"] for i in common)
        if mode == "contended" else None, scientific_timings_admitted=False)
    if read(directory / "trial.json") != expected:
        raise ValueError("Trial summary does not reproduce")
    for item in evidence:
        check(item)
    if inventory(directory) != evidence:
        raise ValueError("Control evidence inventory changed during replay")
    return dict(status="full_node_control_replayed", trial=expected, evidence=evidence,
                measurement=measured, scientific_timings_admitted=False,
                limitations=["Scheduler, pinned recipe and runtime checks require separate validation.",
                             "CPU flags do not identify foreign activity or admit scientific timings."])
