"""Freeze and run a paired candidate-only QfO hit-order experiment."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys


ARMS = ("retained_order", "canonical_order")
INSTALL_SHA = "f1636e6e3222c135d2f174bb3f2c247f74e646d8653d087f068960f4926f2559"


def record(path):
    path = Path(path).absolute()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def validate_plan(path, sha):
    if record(path)["sha256"] != sha:
        raise ValueError("Changed paired plan")
    plan = json.loads(Path(path).read_text())
    if (plan["arms"] != list(ARMS) or plan["attempts_per_arm"] != 1
            or plan["profile"] != "satellite_v2" or plan["accuracy_evaluated"] is not False):
        raise ValueError("Wrong experimental scope")
    for item in plan["checked_records"]:
        if record(item["path"]) != item:
            raise ValueError("Changed pinned record: " + item["path"])
    return plan


def prepare(repo, directory):
    from benchmark_tools.audit_qfo_candidate_order import ADMISSION_SHA
    from benchmark_tools.audit_recovery_install import audit
    if directory.exists():
        raise FileExistsError(directory)
    admission_path = repo / "benchmark_tools/results/qfo_corrected_candidate_admission_21759.json"
    if record(admission_path)["sha256"] != ADMISSION_SHA:
        raise ValueError("Changed candidate admission")
    admission = json.loads(admission_path.read_text())
    prepared_record = admission["prepared_manifest"]
    if record(prepared_record["path"]) != prepared_record:
        raise ValueError("Changed candidate preparation")
    prepared = json.loads(Path(prepared_record["path"]).read_text())
    arm = prepared["candidate_arms"]["p1_c1"]
    checkpoint = prepared["numeric_checkpoint"]["manifest"]
    if record(checkpoint["path"]) != checkpoint:
        raise ValueError("Changed checkpoint")
    data = json.loads(Path(checkpoint["path"]).read_text())
    if data["genes"] != 984137 or data["hits"] != 90687327 or data["complete"] is not True:
        raise ValueError("Wrong checkpoint universe")
    installed = audit(repo, repo / "benchmarks/work/publication_recovery_install_20260926")
    installed_bytes = (json.dumps(installed, indent=2, sort_keys=True) + "\n").encode()
    if hashlib.sha256(installed_bytes).hexdigest() != INSTALL_SHA:
        raise ValueError("Recovery package audit differs")
    checked = [record(admission_path), prepared_record, checkpoint,
        arm["seed_partition"], arm["candidate_partition"], arm["membership_constraints"],
        record(__file__), record(repo / "benchmark_tools/candidate_hit_order_policy.py")]
    for name, identity in data["files"].items():
        checked.append(dict(path=str(Path(checkpoint["path"]).parent / name), **identity))
    site = Path(installed["runtime"]["module"]).parent.parent
    for row in installed["wheels"]:
        for member in row["matched"]:
            checked.append(dict(path=str(site / member["member"]),
                bytes=member["bytes"], sha256=member["sha256"]))
    python = repo / "benchmarks/work/publication_recovery_install_20260926/venv/bin/python"
    checked.append(record(python))
    plan = dict(repo=str(repo), directory=str(directory), python=str(python),
        checked_records=checked, arms=list(ARMS), attempts_per_arm=1, profile="satellite_v2",
        checkpoint=checkpoint, seed=arm["seed_partition"], historical=arm["candidate_partition"],
        historical_constraints=arm["membership_constraints"], expected_genes=data["genes"], expected_species=78,
        environment=dict(HOME=os.environ["HOME"], LANG="C.UTF-8", PATH="/usr/bin:/bin",
            PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"),
        accuracy_evaluated=False, publication_ready=False,
        resource_request=dict(cpus=2, memory_gib=128, hours=4),
        endpoints=["Complete candidate partitions: retained vs canonical and each vs historical",
                   "Ordered semantic source/target membership constraints in the same contrasts"],
        limitations=["Candidate-only cached-input experiment; upstream runtime equivalence untested",
                     "No phylogeny, accuracy scoring, parameter tuning or native retry",
                     "Shared-host resources are descriptive"])
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed preparation record")
    directory.mkdir(parents=True)
    save(directory / "install_preparation.json", installed)
    save(directory / "plan.json", plan)
    return record(directory / "plan.json")


def worker(path, sha, arm):
    # Import installed scientific code before exposing the harness namespace.
    import orthohmm
    from orthohmm import orthohmm as pipeline
    from orthohmm.accuracy import load_accuracy_checkpoint
    import numpy as np
    plan = validate_plan(path, sha)
    if arm not in ARMS or sys.executable != plan["python"]:
        raise ValueError("Wrong arm/interpreter")
    if not Path(orthohmm.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise ValueError("Scientific code must be installed")
    if any(os.environ.get(k) != v for k, v in plan["environment"].items()):
        raise ValueError("Changed native environment")
    target = Path(plan["directory"]) / arm
    target.mkdir(exist_ok=False)
    save(target / "started.json", dict(plan=record(path), arm=arm, executable=sys.executable,
        source=record(__file__), pipeline=record(pipeline.__file__), attempts=1))
    try:
        sys.path.insert(0, plan["repo"])
        from benchmark_tools.candidate_hit_order_policy import canonical_hit_order_v1
        names, species, q, t, s = load_accuracy_checkpoint(Path(plan["checkpoint"]["path"]).parent)
        if len(names) != plan["expected_genes"]:
            raise ValueError("Wrong gene universe")
        hits = (q, t, s)
        if arm == "canonical_order":
            hits = canonical_hit_order_v1(names, *hits)
        working = target / "orthohmm_working_res"
        working.mkdir()
        shutil.copy2(plan["seed"]["path"], working / "orthohmm_edges_clustered.txt")
        details = pipeline._expand_phylogeny_candidates(str(target), names, species, hits,
                                                        profile=plan["profile"])
        details.pop("_membership_constraints", None)
        if len(np.unique(species)) != plan["expected_species"]:
            raise ValueError("Wrong species universe")
        validate_plan(path, sha)
        save(target / "complete.json", dict(status="candidate_arm_complete_pending_readback",
            plan=record(path), arm=arm, details=details, outputs=[record(p) for p in sorted(working.iterdir())
            if p.is_file()], accuracy_evaluated=False))
    except BaseException as error:
        save(target / "failure.json", dict(type=type(error).__name__, error=str(error), retry=False))
        raise


def run(path, sha):
    plan = validate_plan(path, sha)
    directory = Path(plan["directory"])
    if os.environ.get("SLURM_CPUS_PER_TASK") != "2" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled two-CPU execution")
    save(directory / "execution_started.json", dict(plan=record(path), job_id=os.environ["SLURM_JOB_ID"]))
    results = []
    try:
        for arm in ARMS:
            command = ["/usr/bin/time", "-v", "-o", str(directory / (arm + ".time.txt")),
                plan["python"], "-I", str(Path(__file__).resolve()), "worker", "--plan", str(path),
                "--sha256", sha, "--arm", arm]
            with (directory / (arm + ".log")).open("x") as stream:
                process = subprocess.Popen(command, cwd=directory, env=plan["environment"],
                    stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=6600)
                except BaseException:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
            results.append(dict(arm=arm, command=command, returncode=code,
                log=record(directory / (arm + ".log")), time=record(directory / (arm + ".time.txt"))))
            if code:
                raise RuntimeError("Failed arm; no retry")
        validate_plan(path, sha)
        save(directory / "execution.json", dict(status="paired_candidates_complete_pending_readback",
            plan=record(path), job_id=os.environ["SLURM_JOB_ID"], arms=results, accuracy_evaluated=False))
    except BaseException as error:
        save(directory / "execution_failure.json", dict(error=str(error), type=type(error).__name__, arms=results))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "worker", "run"))
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--arm", choices=ARMS)
    args = parser.parse_args()
    if args.mode == "prepare":
        print(json.dumps(prepare(args.repo.resolve(), args.directory.resolve()), indent=2))
    elif args.mode == "worker":
        worker(args.plan.resolve(), args.sha256, args.arm)
    else:
        run(args.plan.resolve(), args.sha256)
