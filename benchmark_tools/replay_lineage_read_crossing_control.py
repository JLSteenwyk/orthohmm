"""Replay all three retained read-crossing controls; never admit scientific timings."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_cgroup_lineage import compare
from benchmark_tools.run_lineage_lifecycle_control import assess as assess_lifecycle
from benchmark_tools.run_lineage_read_crossing_control import assess

COMMIT = "6599c6e"
REMOTE = Path("/home/jlsteenwyk/projects/orthohmm-publication/lineage_read_crossing_source_v1")
SOURCES = {"probe_cgroup_frontier.py", "probe_cgroup_lineage.py", "probe_dgx_cpu_hierarchy.py",
           "probe_dgx_step_separation.py", "probe_host_counters.py", "run_lineage_lifecycle_control.py",
           "run_lineage_read_crossing_control.py", "results/LINEAGE_READ_CROSSING_PROTOCOL_20260919.md"}


def replay_trial(directory, index, identity):
    def read(name):
        return json.loads((directory / name).read_text())

    job, manager, target = (identity[k] for k in ("job", "manager", "target"))
    cpu = min(identity["affinity"])
    unit = f"orthohmm-lineage-{job}-{index}.service"
    command = read("lifecycle/command.json")
    expected = ["systemd-run", "--user", "--quiet", "--wait", "--pipe", "--collect", "--unit=" + unit,
        "--working-directory=" + str(REMOTE), "--property=RuntimeMaxSec=10s",
        "--property=MemoryMax=256M", "--property=TasksMax=8", "taskset", "-c", str(cpu),
        "/usr/bin/python3", "-B", "-m", "benchmark_tools.run_lineage_lifecycle_control", "--burn"]
    if command != dict(command=expected, cpu=cpu, manager=manager, target=target):
        raise ValueError("Owned-service command or resource limits differ")
    service = read("lifecycle/service.json")
    if service["returncode"] != 0:
        raise ValueError("Owned service failed")
    load = json.loads(service["stdout"])
    removal = read("lifecycle/removal.json")
    observations = removal["observations"]
    if (not observations or observations[-1]["exists"] is not False
            or removal["scope"] != load["cgroup"].strip().removeprefix("0::")
            or any(b["monotonic_ns"] <= a["monotonic_ns"] for a, b in zip(observations, observations[1:]))):
        raise ValueError("Invalid service disappearance evidence")
    first, last = read("lifecycle/before.json"), read("lifecycle/after.json")
    for key, scope in (("manager", manager), ("target", target)):
        if first[key]["target"] != scope or last[key]["target"] != scope:
            raise ValueError("Nested observation target differs")
    manager_delta = compare(first["manager"], last["manager"])
    outside_delta = compare(first["target"], last["target"])
    nested = dict(status="lifecycle_control_evaluated", index=index, load=load,
        manager=manager_delta, outside=outside_delta,
        result=assess_lifecycle(manager, target, load, unit, True, manager_delta, outside_delta, cpu))
    if nested != read("lifecycle_result.json"):
        raise ValueError("Nested lifecycle result does not reproduce")
    event = read("event.json")
    if not (event["started_ns"] <= first["manager"]["rows"][0]["started_ns"]
            < observations[-1]["monotonic_ns"] <= last["manager"]["rows"][0]["started_ns"]
            <= last["target"]["rows"][-1]["finished_ns"] <= event["finished_ns"]):
        raise ValueError("Nested service evidence lies outside callback")
    points = [read(name + ".json") for name in ("before", "crossing", "after")]
    if any(p["target"] != target for p in points):
        raise ValueError("Outer observation target differs")
    result = assess(*points, event, nested)
    if result != read("result.json"):
        raise ValueError("Read-crossing result does not reproduce")
    return dict(index=index, result=result, service_cpu_s=load["cpu_seconds"],
                manager_cpu_s=nested["result"]["manager_cpu_s"])


def replay(directory, scheduler, repo):
    evidence = [record(p) for p in sorted(directory.rglob("*")) if p.is_file()]
    if any(p.is_symlink() for p in directory.rglob("*")) or len(evidence) != 35:
        raise ValueError("Require complete nonsymlink 35-file control archive")
    evidence.append(record(scheduler))
    identity = json.loads((directory / "identity.json").read_text())
    summary = json.loads((directory / "report.json").read_text())
    fields = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", scheduler.read_text()))
    expected = dict(JobId=str(identity["job"]), JobState="COMPLETED", ExitCode="0:0", Requeue="0",
        Restarts="0", NodeList="spark-7ff0", NumCPUs="20", OverSubscribe="NO", MinMemoryNode="4G",
        TimeLimit="00:10:00", Command=str(repo / "benchmark_tools/run_dgx_lineage_read_crossing.sh"))
    if any(fields.get(k) != v for k, v in expected.items()) or fields.get("CPUs/Task") != "2":
        raise ValueError("Terminal scheduler identity or allocation differs")
    if (identity["uname"][1] != "spark-7ff0" or identity["affinity"] != [0, 1]
            or f"/job_{identity['job']}/step_0/" not in identity["target"]):
        raise ValueError("Observer identity or affinity differs")
    source_rows = []
    sources = {str(Path(p).relative_to(REMOTE / "benchmark_tools")): sha
               for p, sha in identity["sources"].items()}
    if set(sources) != SOURCES or summary["sources_unchanged"] is not True:
        raise ValueError("Changed or incomplete deployed sources")
    for name, sha in sources.items():
        code = subprocess.check_output(["git", "show", COMMIT + ":benchmark_tools/" + name], cwd=repo)
        if hashlib.sha256(code).hexdigest() != sha:
            raise ValueError("Deployed source differs from frozen commit")
        if hashlib.sha256((repo / "benchmark_tools" / name).read_bytes()).hexdigest() != sha:
            raise ValueError("Local replay dependency differs from frozen source")
        source_rows.append(dict(commit=COMMIT, path="benchmark_tools/" + name, sha256=sha))
    rows = [replay_trial(directory / f"trial_{i}", i, identity) for i in range(3)]
    if (summary["trials"] != [r["result"] for r in rows]
            or summary["all_controls_met"] != all(r["result"]["control_met"] for r in rows)):
        raise ValueError("Complete control summary does not reproduce")
    for item in evidence:
        check(item)
    return dict(status="three_read_crossing_controls_replayed", trials=rows,
        all_controls_met=summary["all_controls_met"], evidence=evidence, source=record(__file__),
        deployed_sources=source_rows, scientific_timings_admitted=False,
        environmental_validity_established=False,
        limitations=["Checks retained observations and source/command identity, not unobserved process activity.",
            "Fixed finite-service engineering control; no specificity, overhead or native timing admission."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("directory", "scheduler", "repo", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = replay(args.directory.resolve(), args.scheduler.resolve(), args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
