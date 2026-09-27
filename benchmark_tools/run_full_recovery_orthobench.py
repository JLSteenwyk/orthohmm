"""Pin and execute one full canonical recovery OrthoBench reproduction."""

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess

from benchmark_tools.audit_recovery_install import audit as install_audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.run_publication_pipeline import save

PARENT_SHA = "5fd8dc70c337951706fe9246c2bf08d1191196939da066c5b74bba0bd3dd2b12"
INSTALL_SHA = "f1636e6e3222c135d2f174bb3f2c247f74e646d8653d087f068960f4926f2559"


def prepare(repo, directory):
    if directory.exists():
        raise FileExistsError(directory)
    parent = repo / "benchmarks/work/publication_installed_orthobench_20260926/plan.json"
    if record(parent)["sha256"] != PARENT_SHA:
        raise ValueError("Changed baseline plan")
    old = json.loads(parent.read_text())
    input_dir = parent.parent / "input"
    inputs = [record(p) for p in sorted(input_dir.iterdir()) if p.is_file()]
    expected = {r["path"]: r for r in old["checked_records"]}
    if (len(inputs) != 12 or any(expected.get(r["path"]) != r for r in inputs)
            or len(fasta_ids([Path(r["path"]) for r in inputs])) != 251378):
        raise ValueError("Changed OrthoBench input universe")
    installation = repo / "benchmarks/work/publication_recovery_install_20260926"
    installed = install_audit(repo, installation)
    mafft_path = repo / "benchmark_tools/results/publication_mafft_build_20260926.json"
    if record(mafft_path)["sha256"] != "4c4e92a29c1dc4b27e4b47f4010e9ea9648c17beb02b5a00995a59535394ac76":
        raise ValueError("Changed MAFFT evidence")
    mafft = json.loads(mafft_path.read_text())
    fasttree = old["command"][old["command"].index("--tree_builder") + 1]
    checked = [record(parent), record(mafft_path), *inputs, mafft["launcher"],
               *mafft["built_helpers"], expected[fasttree], old["baseline_partition"]]
    checked.extend(record(p) for p in sorted((repo / "benchmark_tools").glob("*.py")))
    checked.extend(installed["checked_records"])
    site = Path(installed["runtime"]["module"]).parent.parent
    for row in installed["wheels"]:
        for member in row["matched"]:
            item = record(site / member["member"])
            if (item["bytes"], item["sha256"]) != (member["bytes"], member["sha256"]):
                raise ValueError("Installed wheel member changed")
            checked.append(item)
    for item in checked:
        check(item)
    directory.mkdir(parents=True)
    save(directory / "install_preparation.json", installed)
    if record(directory / "install_preparation.json")["sha256"] != INSTALL_SHA:
        raise ValueError("Recovery installation differs from validated state")
    python = installation / "venv/bin/python"
    command = [str(python), "-I", str(repo / "benchmark_tools/run_publication_pipeline.py"),
               "--input", str(input_dir), "--output", str(directory / "native"),
               "--cpu", "32", "--aligner", mafft["launcher"]["path"], "--tree-builder", fasttree]
    plan = dict(status="prepared_not_run", repo=str(repo), directory=str(directory),
        source=record(__file__), checked_records=checked, command=command,
        installation=str(installation), installation_audit_sha256=INSTALL_SHA,
        environment=dict(HOME=os.environ["HOME"], LANG="C.UTF-8", **old["environment_overrides"]),
        resource_request=dict(cpus=32, memory_gib=128, hours=24), attempts=1,
        expected_genes=251378, expected_species=12, checkpoint_reuse=False,
        baseline_partition=old["baseline_partition"], baseline_scores=old["baseline_scores"],
        scientific_scores_admitted=False, publication_ready=False,
        endpoint="Exact historical root partition and all 70 frozen RefOG score records; report any mismatch.",
        timing_scope="Shared-host reproduction, not controlled comparative timing.")
    save(directory / "plan.json", plan)
    return record(directory / "plan.json")


def validate_plan(path, sha):
    if record(path)["sha256"] != sha:
        raise ValueError("Changed or unpinned plan")
    plan = json.loads(path.read_text())
    if (plan["attempts"] != 1 or plan["checkpoint_reuse"] is not False
            or plan["scientific_scores_admitted"] is not False
            or plan["expected_genes"] != 251378 or plan["expected_species"] != 12
            or plan["resource_request"] != dict(cpus=32, memory_gib=128, hours=24)):
        raise ValueError("Wrong full reproduction scope")
    if plan["source"] != record(__file__):
        raise ValueError("Changed executor")
    for item in plan["checked_records"]:
        check(item)
    return plan


def run(path, sha):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32":
        raise ValueError("Require the specified 32-CPU Slurm allocation")
    plan = validate_plan(path, sha)
    directory = Path(plan["directory"])
    if any((directory / n).exists() for n in ("native", "execution_started.json", "execution.json")):
        raise FileExistsError("Existing attempt; no automatic retry/resume")
    save(directory / "execution_started.json", dict(plan=record(path), job_id=os.environ["SLURM_JOB_ID"]))
    command = ["/usr/bin/time", "-v", "-o", str(directory / "time.txt"), *plan["command"]]
    result = dict(status="running", plan=record(path), job_id=os.environ["SLURM_JOB_ID"], command=command)
    child = None
    try:
        before = install_audit(Path(plan["repo"]), Path(plan["installation"]))
        save(directory / "install_before.json", before)
        if record(directory / "install_before.json")["sha256"] != plan["installation_audit_sha256"]:
            raise ValueError("Installation changed before execution")
        with (directory / "native.log").open("x") as log:
            child = subprocess.Popen(command, cwd=directory, env=plan["environment"],
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            result["returncode"] = child.wait(timeout=23 * 3600)
        if result["returncode"]:
            raise RuntimeError("Native pipeline failed")
        validate_plan(path, sha)
        completed = json.loads((directory / "native/complete.json").read_text())
        if completed["status"] != "native_complete_pending_scientific_readback":
            raise ValueError("Missing native completion")
        after = install_audit(Path(plan["repo"]), Path(plan["installation"]))
        save(directory / "install_after.json", after)
        if record(directory / "install_after.json")["sha256"] != plan["installation_audit_sha256"]:
            raise ValueError("Installation changed during execution")
        result.update(status="native_complete_pending_independent_readback", native_complete=record(directory / "native/complete.json"))
    except BaseException as error:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
        result.update(status="failed", error_type=type(error).__name__, error=str(error), retry=False)
        raise
    finally:
        result["logs"] = [record(directory / n) for n in ("native.log", "time.txt") if (directory / n).exists()]
        save(directory / "execution.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--prepare", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--plan-sha256")
    args = parser.parse_args()
    if args.prepare and args.repo and not args.run:
        print(json.dumps(prepare(args.repo.resolve(), args.prepare.resolve())))
    elif args.run and args.plan_sha256 and not args.prepare:
        run(args.run.resolve(), args.plan_sha256)
    else:
        parser.error("Require --repo/--prepare or --run/--plan-sha256")
