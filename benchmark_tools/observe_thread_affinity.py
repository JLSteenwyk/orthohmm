"""Read-only, non-atomic cgroup thread-affinity evidence for timing collectors."""

import os
from pathlib import Path
import time


def identity(path):
    raw = path.read_text()
    # comm can contain whitespace and parentheses; fields after its final ')' are fixed.
    pid = int(raw.split("(", 1)[0].strip())
    fields = raw[raw.rindex(")") + 1:].split()
    return pid, int(fields[19])  # field 22: starttime


def inventory(scope):
    tasks = set()
    errors = []
    def failed(error):
        errors.append(dict(path=str(error.filename), error=str(error)))
    for root, _, _ in os.walk(scope, onerror=failed):
        path = Path(root) / "cgroup.threads"
        try:
            values = [int(v) for v in path.read_text().split()]
            if any(v <= 0 for v in values):
                raise ValueError("Nonpositive thread ID")
            tasks.update(values)
        except (OSError, ValueError) as error:
            errors.append(dict(path=str(path), error=str(error)))
    if not Path(scope).is_dir():
        errors.append(dict(path=str(scope), error="Missing scope"))
    return sorted(tasks), errors


def observe(scope, allowed_cpus, *, cgroup_root=Path("/sys/fs/cgroup"),
            proc_root=Path("/proc"), get_affinity=os.sched_getaffinity):
    """Sample every enumerated thread, including threads in descendant cgroups.

    A subset affinity is valid (workers may pin themselves). Races, inaccessible
    threads and membership changes are gaps, never evidence of compliance.
    """
    allowed = list(allowed_cpus)
    if (not allowed or any(type(v) is not int or v < 0 for v in allowed)
            or len(set(allowed)) != len(allowed)):
        raise ValueError("Require distinct nonnegative CPU IDs")
    cgroup_root = Path(cgroup_root).resolve(strict=True)
    scope = Path(scope).resolve(strict=True)
    relative = scope.relative_to(cgroup_root)
    if relative == Path("."):
        raise ValueError("Refuse whole-host scope")
    expected = "/" + relative.as_posix()
    start = time.monotonic_ns()
    tids, errors = inventory(scope)
    threads = []
    for tid in tids:
        path = Path(proc_root) / str(tid)
        try:
            before = identity(path / "stat")
            membership_before = (path / "cgroup").read_text()
            affinity = sorted(get_affinity(tid))
            membership_after = (path / "cgroup").read_text()
            after = identity(path / "stat")
            memberships = [line[3:] for line in membership_before.splitlines()
                           if line.startswith("0::")]
            if before != after or before[0] != tid:
                raise ValueError("Thread identity changed")
            if (membership_before != membership_after or len(memberships) != 1
                    or not (memberships[0] == expected or memberships[0].startswith(expected + "/"))):
                raise ValueError("Thread scope changed or mismatched")
            if not affinity:
                raise ValueError("Empty affinity")
            threads.append(dict(tid=tid, start_ticks=before[1], cgroup=membership_before,
                                affinity=affinity, outside_cpus=sorted(set(affinity) - set(allowed))))
        except (OSError, ValueError, IndexError) as error:
            errors.append(dict(tid=tid, error=str(error)))
    final_tids, final_errors = inventory(scope)
    errors.extend(final_errors)
    if tids != final_tids:
        errors.append(dict(error="Thread inventory changed", before=tids, after=final_tids))
    if not tids:
        errors.append(dict(error="Empty thread inventory"))
    violations = [row["tid"] for row in threads if row["outside_cpus"]]
    return dict(scope=str(scope), allowed_cpus=sorted(allowed), started_ns=start,
                finished_ns=time.monotonic_ns(), initial_tids=tids, final_tids=final_tids,
                threads=threads, errors=errors, violating_tids=violations,
                status="violation" if violations else "incomplete" if errors else "observed_within_affinity",
                full_run_affinity_verified=False,
                limitations=["Non-atomic periodic evidence misses short-lived threads and between-sample changes.",
                             "Affinity is not a CPU quota; this observation does not establish host isolation."])
