"""Read a disjoint cgroup frontier around one scope; never admit timing."""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_dgx_cpu_hierarchy import usage


def validate_scope(scope):
    path = PurePosixPath(scope)
    if not scope.startswith("/") or "//" in scope or str(path) != scope or ".." in path.parts or scope == "/":
        raise ValueError("Require normalized non-root absolute cgroup scope")
    return path


def inventory(root, target):
    path = validate_scope(target)
    ancestors = [PurePosixPath("/"), *list(reversed(path.parents))[1:]]
    scopes = {target}
    direct = {}
    for ancestor in ancestors:
        directory = root / str(ancestor).lstrip("/")
        if directory.is_symlink():
            raise ValueError("Symlink in cgroup hierarchy")
        direct[str(ancestor)] = len((directory / "cgroup.procs").read_text().splitlines())
        for child in directory.iterdir():
            if child.is_symlink():
                raise ValueError("Symlink in cgroup hierarchy")
            if child.is_dir():
                scope = ancestor / child.name
                if scope != path and scope not in path.parents:
                    scopes.add(str(scope))
    identities = {}
    for scope in sorted(scopes):
        directory = root / scope.lstrip("/")
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError("Missing or symlinked cgroup")
        stat = directory.stat()
        identities[scope] = [stat.st_dev, stat.st_ino]
    return dict(identities=identities, ancestor_direct_process_counts=direct)


def read_counter(root, scope):
    start = time.monotonic_ns()
    raw = (root / scope.lstrip("/") / "cpu.stat").read_text()
    end = time.monotonic_ns()
    usage(raw)
    return dict(scope=scope, started_ns=start, finished_ns=end, raw=raw)


def snapshot(root, target):
    before = inventory(root, target)
    root_before = read_counter(root, "/")
    rows = [read_counter(root, scope) for scope in before["identities"]]
    root_after = read_counter(root, "/")
    after = inventory(root, target)
    if before["identities"] != after["identities"]:
        raise ValueError("Cgroup frontier changed during sampling")
    return dict(target=target, boot_id=Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
                root=[root_before, root_after], rows=rows,
                inventory_before=before, inventory_after=after)


def validate(point):
    target = validate_scope(point["target"])
    ids = point["inventory_before"]["identities"]
    if ids != point["inventory_after"]["identities"]:
        raise ValueError("Cgroup identity changed")
    names = [row["scope"] for row in point["rows"]]
    if names != sorted(ids) or str(target) not in names or len(names) != len(set(names)):
        raise ValueError("Invalid frontier inventory")
    paths = [validate_scope(name) for name in names]
    if any(a in b.parents for a in paths for b in paths):
        raise ValueError("Overlapping frontier scopes")
    if any(row["scope"] != "/" for row in point["root"]) or len(point["root"]) != 2:
        raise ValueError("Invalid root bracket")
    previous = -1
    for row in [point["root"][0], *point["rows"], point["root"][1]]:
        if not previous <= row["started_ns"] <= row["finished_ns"]:
            raise ValueError("Counter read order invalid")
        usage(row["raw"])
        previous = row["finished_ns"]


def compare(left, right):
    for point in (left, right):
        validate(point)
    if (left["boot_id"] != right["boot_id"] or left["target"] != right["target"]
            or left["inventory_before"]["identities"] != right["inventory_before"]["identities"]):
        raise ValueError("Boot, target or frontier identity changed")
    if left["root"][1]["finished_ns"] >= right["root"][0]["started_ns"]:
        raise ValueError("Observation windows overlap")
    def delta(a, b):
        value = usage(b["raw"]) - usage(a["raw"])
        if value < 0:
            raise ValueError("Counter decreased")
        return value / 1e6
    scopes = {a["scope"]: delta(a, b) for a, b in zip(left["rows"], right["rows"])}
    root_delta = delta(left["root"][0], right["root"][1])
    return dict(scope_cpu_s=scopes, target_cpu_s=scopes[left["target"]],
                outside_target_frontier_cpu_s=sum(v for k, v in scopes.items() if k != left["target"]),
                root_outer_cpu_s=root_delta, root_minus_frontier_cpu_s=root_delta-sum(scopes.values()),
                scientific_timings_admitted=False, controlled_workload_verified=False,
                limitations=["Non-atomic counters; differences are not bounds or timing corrections.",
                    "Ancestor-direct tasks and accounting discrepancies remain in the signed root residual.",
                    "Stable endpoints cannot exclude a transient cgroup born and removed between samples.",
                    "Scope CPU activity does not establish causal interference with a native command."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root = Path("/sys/fs/cgroup")
    sources = [Path(__file__).resolve(), *[Path(__file__).with_name(name) for name in (
        "probe_dgx_cpu_hierarchy.py", "probe_host_counters.py", "probe_dgx_step_separation.py")]]
    hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    left = snapshot(root, args.target)
    time.sleep(2)
    right = snapshot(root, args.target)
    if hashes != {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}:
        raise ValueError("Probe source changed during observation")
    result = dict(host=os.uname().nodename, sources=hashes, points=[left, right], result=compare(left, right))
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
