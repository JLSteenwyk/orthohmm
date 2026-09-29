"""One owned Slurm diagnostic of parent exit versus live descendant detection."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.observe_thread_affinity import observe, identity


def run(output):
    if output.exists():
        raise FileExistsError(output)
    job = int(os.environ["SLURM_JOB_ID"])
    membership = Path("/proc/self/cgroup").read_text()
    rows = [line[3:] for line in membership.splitlines() if line.startswith("0::")]
    if len(rows) != 1 or f"job_{job}" not in Path(rows[0]).parts:
        raise ValueError("Wrong job membership")
    scope = Path("/sys/fs/cgroup") / rows[0].lstrip("/")
    while scope.name != "user":
        if scope == scope.parent or scope.name.startswith("step_"):
            raise ValueError("Missing native user scope")
        scope = scope.parent
    allowed = sorted(os.sched_getaffinity(0))
    before = observe(scope, allowed)
    subprocess.run(["/bin/true"], check=True)
    normal_finished = time.monotonic_ns()
    normal = completion_evidence(before, observe(scope, allowed), os.getpid(), normal_finished)
    if normal["errors"]:
        raise ValueError("Initial clean control failed")
    command = [sys.executable, "-I", "-B", "-c",
        "import subprocess; p=subprocess.Popen(['/bin/sleep','20'], start_new_session=True, "
        "stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); print(p.pid, flush=True)"]
    parent = subprocess.run(command, check=True, capture_output=True, text=True, start_new_session=True)
    finished = time.monotonic_ns()
    pid = int(parent.stdout.strip())
    fd = os.pidfd_open(pid)
    try:
        child_identity = identity(Path(f"/proc/{pid}/stat"))
        child_membership = Path(f"/proc/{pid}/cgroup").read_text()
        if child_membership != membership or child_identity[0] != pid:
            raise ValueError("Owned child identity or scope changed")
        escaped_session = os.getsid(pid) == pid
        after = observe(scope, allowed)
        descendant = completion_evidence(before, after, os.getpid(), finished)
    finally:
        signal.pidfd_send_signal(fd, signal.SIGTERM)
        os.close(fd)
    result = dict(job_id=job, host=os.uname().nodename, membership=membership,
        command=command, parent_exit_code=parent.returncode, child_pid=pid,
        child_identity=child_identity, child_membership=child_membership,
        child_started_new_session=escaped_session, normal=normal, descendant=descendant,
        negative_control_passed=bool(descendant["errors"]) and pid in after["final_tids"] and escaped_session,
        scientific_timings_admitted=False, complete_job_accounting_verified=False,
        sources=[dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for path in (Path(__file__).resolve(), Path(__file__).with_name("native_completion.py").resolve(),
                              Path(__file__).with_name("observe_thread_affinity.py").resolve())])
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not result["negative_control_passed"]:
        raise ValueError("Live descendant negative control failed")
    print(json.dumps(dict(job_id=job, clean_control=True, descendant_rejected=True)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output)
