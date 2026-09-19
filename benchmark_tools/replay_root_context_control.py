"""Replay raw control witnesses and both observation streams, without timing admission."""

from pathlib import Path

from benchmark_tools.full_node_control_workload import require_start
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.replay_full_node_control import inventory
from benchmark_tools.replay_native_root_context import replay as replay_measurement
from benchmark_tools.root_context_control_design import ORDER, service_command, unit_name
from benchmark_tools.run_root_context_control_trial import read, validate_trial
from benchmark_tools.validate_full_node_control import integer


def replay(directory, remote_directory, recipe_directory, python, mode, job, index, manager):
    directory, remote_directory, recipe_directory = map(Path, (directory, remote_directory, recipe_directory))
    unit = unit_name(job, index)
    if mode != [value for block in ORDER for value in block][index]:
        raise ValueError("Condition differs from frozen index")
    if not all(path.is_absolute() for path in (remote_directory, recipe_directory, Path(python))):
        raise ValueError("Require absolute expected deployment paths")
    evidence = inventory(directory)
    work = directory / "workload"
    source = recipe_directory / "benchmark_tools/full_node_control_workload.py"
    command = ["/usr/bin/sleep", "20"] if mode == "idle" else [python, "-B", str(source),
        "--directory", str(remote_directory / "workload"), "--worker", "churn" if mode == "churn" else "steady"]
    measured = replay_measurement(directory / "measurement", job, command, expected_timeout_s=60)
    report = measured["lineage"]["measured"]
    controller = read(directory / "controller.json")
    if mode == "idle":
        if controller != {}:
            raise ValueError("Idle controller evidence differs")
    else:
        require_start(read(work / "workload_go.json"))
        ready, done = read(work / "workload_ready.json"), read(work / "workload_done.json")
        names = {f"{kind}_{cpu}.json" for cpu in done["cpus"] for kind in ("ready", "done")}
        names.update({"workload_ready.json", "workload_done.json", "workload_go.json"})
        if mode == "user-contended":
            names.update({"competitor_ready.json", "competitor_done.json", "service.log"})
        if {path.name for path in work.iterdir()} != names:
            raise ValueError("Raw workload inventory differs")
        for before, after in zip(ready["workers"], done["workers"]):
            if (before != read(work / f"ready_{after['cpu']}.json")
                    or after != read(work / f"done_{after['cpu']}.json")):
                raise ValueError("Raw and embedded worker witnesses disagree")
        if mode == "user-contended":
            expected = service_command(python, source, remote_directory / "workload", unit, min(ready["cpus"]))
            if (controller.get("command") != expected or controller.get("unit") != unit
                    or controller.get("manager") != manager
                    or type(controller.get("service_exit_code")) is not int
                    or "owned_service_stop" in controller or "cleanup_error" in controller):
                raise ValueError("Owned-service command, identity or exit evidence differs")
            removal = controller.get("removal", [])
            finish = read(work / "competitor_done.json")["finished_ns"]
            times = [row["monotonic_ns"] for row in removal]
            if (not times or not all(integer(t) for t in times)
                    or times != sorted(set(times)) or times[0] < finish
                    or [row["exists"] for row in removal] != [True] * (len(removal)-1) + [False]
                    or any(type(row["exists"]) is not bool for row in removal)):
                raise ValueError("Owned-service removal chronology differs")
        elif controller != {"status": "completed"}:
            raise ValueError("Unexpected native-only controller evidence")
    expected = validate_trial(directory, mode, job, index, manager, report, dict(context=measured["context"]))
    if read(directory / "trial.json") != expected:
        raise ValueError("Trial summary does not reproduce")
    for item in evidence:
        check(item)
    if inventory(directory) != evidence:
        raise ValueError("Control evidence changed during replay")
    return dict(status="root_context_control_replayed", trial=expected, measurement=measured, evidence=evidence,
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Deployment source, runtime, scheduler and complete panel inventory require separate audit.",
            "User-slice response does not identify all CPU sources or establish scientific timing validity."])
