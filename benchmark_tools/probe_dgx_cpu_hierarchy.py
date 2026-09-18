"""Read-only Slurm hierarchy controls; no scientific timing admission."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_host_counters import snapshot, summarize, parse_group
from benchmark_tools.probe_dgx_step_separation import save, wait_file, scope_parts, validate_scopes


def usage(raw):
    values = {}
    for line in raw.splitlines():
        key, value = line.split()
        if key in values or int(value) < 0:
            raise ValueError("Invalid CPU counter")
        values[key] = int(value)
    return values["usage_usec"]


def read_cpu(scope):
    started = time.monotonic_ns()
    raw = (Path("/sys/fs/cgroup") / scope.lstrip("/") / "cpu.stat").read_text()
    result = dict(scope=scope, started_ns=started, finished_ns=time.monotonic_ns(), raw=raw)
    usage(raw)
    return result


def read_point(pid, native_membership, job):
    if Path(f"/proc/{pid}/cgroup").read_text() != native_membership:
        raise ValueError("Native process changed scope")
    before = snapshot()
    native = dict(raw=dict(cgroup_membership=native_membership))
    validate_scopes(before, native, job)
    native_scope = Path(*scope_parts(native, job))
    parent = native_scope.parent
    directory = Path("/sys/fs/cgroup") / str(parent).lstrip("/")
    names = sorted(p.name for p in directory.iterdir() if p.is_dir())
    if any(not name.startswith("step_") for name in names):
        raise ValueError("Unexpected immediate job child")
    parent_before = read_cpu(str(parent))
    children = [read_cpu(str(parent / name)) for name in names]
    parent_after = read_cpu(str(parent))
    after = snapshot()
    if names != sorted(p.name for p in directory.iterdir() if p.is_dir()):
        raise ValueError("Job hierarchy changed during read")
    if Path(f"/proc/{pid}/cgroup").read_text() != native_membership:
        raise ValueError("Native process disappeared or changed scope")
    result = dict(host=[before, after], parent=[parent_before, parent_after], children=children,
                  native_membership=native_membership, ticks=os.sysconf("SC_CLK_TCK"))
    validate_point(result, job)
    return result


def validate_point(point, job):
    before, after = point["host"]
    if before["errors"] or after["errors"]:
        raise ValueError("Host counter errors")
    summarize(before, after, point["ticks"])
    native = dict(raw=dict(cgroup_membership=point["native_membership"]))
    validate_scopes(before, native, job)
    validate_scopes(after, native, job)
    path = Path(*scope_parts(native, job))
    a, b = point["parent"]
    if a["scope"] != str(path.parent) or b["scope"] != a["scope"]:
        raise ValueError("Wrong job parent")
    children = point["children"]
    names = [c["scope"] for c in children]
    if (len(names) != len(set(names)) or names != sorted(names)
            or str(path) not in names or str(path.parent / "step_batch") not in names
            or any(Path(n).parent != path.parent or not Path(n).name.startswith("step_") for n in names)):
        raise ValueError("Invalid disjoint immediate step inventory")
    previous = before["finished_monotonic_ns"]
    for row in [a, *children, b]:
        start, end = row["started_ns"], row["finished_ns"]
        if type(start) is not int or type(end) is not int or not previous <= start <= end:
            raise ValueError("Counter reads are not enclosed and ordered")
        usage(row["raw"])
        previous = end
    if previous > after["started_monotonic_ns"] or usage(b["raw"]) < usage(a["raw"]):
        raise ValueError("Invalid parent counter bracket")


def compare(left, right, job):
    for p in (left, right):
        validate_point(p, job)
    if (left["native_membership"] != right["native_membership"] or left["ticks"] != right["ticks"]
            or [c["scope"] for c in left["children"]] != [c["scope"] for c in right["children"]]):
        raise ValueError("Hierarchy changed between observations")
    host = summarize(left["host"][0], right["host"][1], left["ticks"])
    if left["host"][1]["finished_monotonic_ns"] >= right["host"][0]["started_monotonic_ns"]:
        raise ValueError("Observation points overlap")
    def delta(a, b):
        if a["scope"] != b["scope"] or usage(b["raw"]) < usage(a["raw"]):
            raise ValueError("Scope changed or CPU decreased")
        return (usage(b["raw"])-usage(a["raw"])) / 1e6
    steps = {Path(a["scope"]).name: delta(a, b) for a, b in zip(left["children"], right["children"])}
    outer = delta(left["parent"][0], right["parent"][1])
    inner = delta(left["parent"][1], right["parent"][0])
    return dict(host_busy_cpu_s=host["accounted_host_busy_cpu_s"], step_cpu_s=steps,
                job_outer_cpu_s=outer, job_inner_cpu_s=inner,
                job_outer_minus_step_sum_cpu_s=outer-sum(steps.values()),
                host_minus_job_outer_cpu_s=host["accounted_host_busy_cpu_s"]-outer,
                scientific_timings_admitted=False, controlled_workload_verified=False)


def run(output):
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require DGX two-CPU engineering allocation")
    job = int(os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    helper = source.with_name("probe_dgx_step_separation.py")
    paths = [source, helper, source.with_name("probe_host_counters.py")]
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    trials = []
    for mode, repeats in (("quiet", 0), ("completed_burst", 1), ("sustained_batch", 4)):
        directory = output / mode
        directory.mkdir()
        command = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                   sys.executable, "-I", "-B", str(helper), "--worker", str(directory)]
        with (directory / "step.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                ready = wait_file(directory / "ready.json")
                membership = wait_file(directory / "native_before.json")["raw"]["cgroup_membership"]
                before = read_point(ready["pid"], membership, job)
                loads = []
                if not repeats:
                    time.sleep(1)
                for _ in range(repeats):
                    completed = subprocess.run([sys.executable, "-I", "-B", str(helper), "--burn"],
                        check=True, capture_output=True, text=True, timeout=15)
                    load = json.loads(completed.stdout)
                    if (parse_group(load["cgroup"]) != parse_group(before["host"][0]["raw"]["cgroup_membership"])
                            or not .75 <= load["cpu_seconds"] < .85):
                        raise ValueError("Unexpected batch load identity or duration")
                    loads.append(load)
                after = read_point(ready["pid"], membership, job)
                result = compare(before, after, job)
                native_name = Path(*scope_parts(dict(raw=dict(cgroup_membership=membership)), job)).name
                met = (None if not repeats else result["step_cpu_s"]["step_batch"] >= .5 * repeats
                       and result["step_cpu_s"][native_name] < .25)
                save(directory / "trial.json", dict(mode=mode, command=command, points=[before, after],
                     loads=loads, result=result, load_localization_expectation_met=met))
                trials.append(json.loads((directory / "trial.json").read_text()))
            finally:
                if not (directory / "release.json").exists():
                    save(directory / "release.json", {"release": True})
                if process.wait(timeout=45) != 0:
                    raise RuntimeError("Native sleeping step failed")
    if hashes != {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:
        raise ValueError("Probe sources changed")
    report = dict(status="hierarchical_cpu_controls_evaluated", job_id=job, sources=hashes, trials=trials,
        scientific_timings_admitted=False, controlled_workload_verified=False, publication_ready=False,
        limitations=["Three fixed-order controls with a sleeping native step, not native workload calibration.",
            "Parent and child counters have different read windows and accounting delays; residuals are not foreign-load bounds.",
            "Completed loads are in the same batch allocation, not unrelated user jobs.",
            "No timing threshold changed, no native time correction, no scientific panel admission."])
    save(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output.absolute())
