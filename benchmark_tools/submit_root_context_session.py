"""Submit one fresh panel while a bounded SSH session waits for its termination."""

import argparse
from pathlib import Path
import re
import shlex
import subprocess
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")
RECIPE = ROOT / "root_context_controls_recipe_v2"


def command(recipe_sha):
    if not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Require pinned recipe SHA-256")
    remote = ["timeout", "--signal=TERM", "--kill-after=10s", "1020s", "sbatch", "--wait", "--parsable",
        "--chdir=" + str(RECIPE), str(RECIPE / "benchmark_tools/run_dgx_root_context_session.sh"), recipe_sha]
    return ["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "jlsteenwyk@10.10.10.2", shlex.join(remote)]


def decoded(value):
    return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""


def job_id(stdout):
    match = re.fullmatch(r"([1-9][0-9]*)(?:;[A-Za-z0-9_.-]+)?\s*", stdout)
    return int(match[1]) if match else None


def run(output, recipe_sha):
    launched = command(recipe_sha)
    source = record(__file__)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="queue_check_pending", recipe_sha256=recipe_sha, source=source,
        scientific_timings_admitted=False, scheduler_terminal_verified=False)
    try:
        queue_command = ["squeue", "-h", "-p", "spark", "-o", "%i %T %j"]
        start = time.time_ns()
        queue = subprocess.run(queue_command, capture_output=True, text=True, timeout=15)
        save(output / "queue.json", dict(command=queue_command, started_unix_ns=start,
            finished_unix_ns=time.time_ns(), returncode=queue.returncode, stdout=queue.stdout, stderr=queue.stderr))
        if queue.returncode != 0 or queue.stdout.strip():
            result["status"] = "not_submitted_queue_not_verified_empty"
        else:
            save(output / "launch.json", dict(command=launched, started_unix_ns=time.time_ns(),
                recipe_sha256=recipe_sha, local_timeout_s=1050, remote_timeout_s=1020,
                remote_kill_grace_s=10, scientific_timings_admitted=False))
            completed = subprocess.run(launched, capture_output=True, text=True, timeout=1050)
            result.update(status="wait_returned", returncode=completed.returncode,
                stdout=completed.stdout, stderr=completed.stderr, job_id=job_id(completed.stdout),
                finished_unix_ns=time.time_ns())
    except subprocess.TimeoutExpired as error:
        result.update(status="observation_timeout", stdout=decoded(error.stdout), stderr=decoded(error.stderr),
            job_id=job_id(decoded(error.stdout)), finished_unix_ns=time.time_ns(),
            error="Timeout does not establish job termination; inspect Slurm, do not resubmit.")
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
