"""Execute one pinned 32-CPU integrated OrthoBench reproducibility job."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess


def record(path):
    path = Path(path).absolute()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def validate(path, digest):
    if record(path)["sha256"] != digest:
        raise ValueError("Changed integrated run plan")
    plan = json.loads(path.read_text())
    if (plan["dataset"] != "orthobench" or plan["attempts"] != 1
            or plan["checkpoint_reuse"] is not False or plan["controlled_timing"] is not False
            or plan["resources"] != dict(cpus=32, memory_gib=128, hours=24, node="bizon")
            or plan["expected_genes"] != 251378 or plan["expected_refogs"] != 70):
        raise ValueError("Changed prespecified full-data scope")
    if record(__file__) != plan["launcher"]:
        raise ValueError("Changed launcher")
    for row in plan["checked_records"]:
        if record(row["path"]) != row:
            raise ValueError("Changed frozen dependency: " + row["path"])
    return plan


def run(path, digest):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32":
        raise ValueError("Require a 32-CPU Slurm allocation")
    plan = validate(path, digest)
    directory = path.parent
    if (directory / "run").exists() or (directory / "execution_started.json").exists():
        raise FileExistsError("Existing attempt; no retry or resume")
    save(directory / "execution_started.json", dict(plan=record(path), job_id=os.environ["SLURM_JOB_ID"]))
    command = ["/usr/bin/time", "-v", "-o", str(directory / "time.txt"), *plan["command"]]
    result = dict(status="running", plan=record(path), job_id=os.environ["SLURM_JOB_ID"], command=command)
    child = None
    try:
        with (directory / "workflow.log").open("x") as log:
            child = subprocess.Popen(command, cwd=directory, env=plan["environment"],
                                     stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            result["returncode"] = child.wait(timeout=23 * 3600)
        if result["returncode"]:
            raise RuntimeError("Integrated workflow failed")
        validate(path, digest)
        complete = directory / "run/complete.json"
        if json.loads(complete.read_text())["status"] != "integrated_install_inference_readback_scoring_complete":
            raise ValueError("Missing integrated completion")
        result.update(status="integrated_complete_pending_independent_admission", complete=record(complete))
    except BaseException as error:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
        result.update(status="failed", error_type=type(error).__name__, error=str(error), retry=False)
        raise
    finally:
        result["logs"] = [record(directory / n) for n in ("workflow.log", "time.txt") if (directory / n).exists()]
        save(directory / "execution.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    run(args.plan.resolve(), args.sha256)
