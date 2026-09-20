"""Hold one bounded submitting session for the frozen five-hour overhead panel."""

import argparse
from pathlib import Path
import re
import shlex
import subprocess
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.submit_root_context_session import decoded, job_id
from benchmark_tools.verify_root_context_overhead_provenance import RECIPE_ROOT, SUBMISSION_SCRIPT

LOCAL_TIMEOUT = 18150
REMOTE_TIMEOUT = 18120
KILL_GRACE = 10
QUEUE_COMMAND = ["squeue", "-h", "-p", "spark", "-o", "%i %T %j"]


def command(recipe_sha):
    if not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Require pinned recipe SHA-256")
    remote = ["timeout", "--signal=TERM", f"--kill-after={KILL_GRACE}s", f"{REMOTE_TIMEOUT}s",
        "sbatch", "--wait", "--parsable", "--chdir=" + str(RECIPE_ROOT), str(SUBMISSION_SCRIPT), recipe_sha]
    return ["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
            "jlsteenwyk@10.10.10.2", shlex.join(remote)]


def run(output, recipe_sha):
    launched = command(recipe_sha)
    source = record(__file__)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="queue_check_pending", recipe_sha256=recipe_sha, source=source,
        scientific_timings_admitted=False, scheduler_terminal_verified=False)
    try:
        started = time.time_ns()
        queue = subprocess.run(QUEUE_COMMAND, capture_output=True, text=True, timeout=15)
        save(output / "queue.json", dict(command=QUEUE_COMMAND, started_unix_ns=started,
            finished_unix_ns=time.time_ns(), returncode=queue.returncode, stdout=queue.stdout, stderr=queue.stderr))
        if queue.returncode != 0 or queue.stdout.strip():
            result["status"] = "not_submitted_queue_not_verified_empty"
        else:
            save(output / "launch.json", dict(command=launched, started_unix_ns=time.time_ns(),
                recipe_sha256=recipe_sha, local_timeout_s=LOCAL_TIMEOUT, remote_timeout_s=REMOTE_TIMEOUT,
                remote_kill_grace_s=KILL_GRACE, scientific_timings_admitted=False))
            waited = subprocess.run(launched, capture_output=True, text=True, timeout=LOCAL_TIMEOUT)
            result.update(status="wait_returned", returncode=waited.returncode, stdout=waited.stdout,
                stderr=waited.stderr, job_id=job_id(waited.stdout), finished_unix_ns=time.time_ns())
    except subprocess.TimeoutExpired as error:
        result.update(status="observation_timeout", stdout=decoded(error.stdout), stderr=decoded(error.stderr),
            job_id=job_id(decoded(error.stdout)), finished_unix_ns=time.time_ns(),
            error="Observation timeout is not terminal evidence; inspect the same job, never automatically resubmit.")
    except OSError as error:
        result.update(status="observation_error", error_type=type(error).__name__, error=str(error),
                      finished_unix_ns=time.time_ns())
    finally:
        save(output / "result.json", result)
    check(source)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    args = parser.parse_args()
    result = run(args.output.resolve(), args.recipe_sha)
    print(result)
    raise SystemExit(0 if result["status"] == "wait_returned" and result["returncode"] == 0 else 1)
