"""Bounded two-distribution replay through pre-phylogeny candidate grouping."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import venv


def driver(repo, directory, checkpoint, inputs, cpu):
    # Load the installed scientific package before exposing the audit executor.
    import orthohmm
    scientific_path = Path(orthohmm.__file__).resolve().parent
    sys.path.insert(0, str(repo))
    import numpy as np
    import importlib.metadata as metadata
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    from benchmark_tools.probe_installed_ob_clustering import write_json
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from orthohmm.accuracy import load_accuracy_checkpoint
    import orthohmm.externals as externals
    import igraph._igraph
    import leidenalg._c_leiden
    if Path(externals.__file__).resolve().parent != scientific_path:
        raise ValueError("Scientific imports escaped installed package")
    runtime = dict(executable=sys.executable, python=sys.version,
        versions={n:metadata.version(n) for n in ("numpy","igraph","leidenalg")},
        sources=[record(p) for p in sorted(scientific_path.rglob("*.py"))],
        native=[record(igraph._igraph.__file__),record(leidenalg._c_leiden.__file__)])
    child_code = ("import json,importlib.metadata as m; import orthohmm.externals,leidenalg._c_leiden; "
                  "print(json.dumps([m.version('leidenalg'),orthohmm.externals.__file__,leidenalg._c_leiden.__file__]))")
    child = json.loads(subprocess.check_output([sys.executable,"-c",child_code],text=True,timeout=60))
    if child != [runtime["versions"]["leidenalg"],externals.__file__,leidenalg._c_leiden.__file__]:
        raise ValueError("Clustering child import resolution differs from driver")
    runtime["plain_child_import_probe"] = child
    write_json(directory / "runtime.json",runtime)
    original = replay.execute_leiden
    snapshots = []
    def observed(*args, **kwargs):
        index = len(snapshots)
        if index >= 4:
            raise ValueError("Unexpected extra clustering call")
        edges = kwargs["edges"]
        target = directory / f"clustering_{index}"
        target.mkdir()
        np.savez(target / "graph.npz",sources=edges.sources,targets=edges.targets,weights=edges.weights)
        write_json(target / "before.json",dict(graph=record(target / "graph.npz"),
                   resolution=args[0],include_isolates=kwargs["include_isolates"],seed=kwargs["seed"]))
        original(*args, **kwargs)
        out = Path(args[1]) / "orthohmm_working_res/orthohmm_edges_clustered.txt"
        shutil.copyfile(out,target / "partition.txt")
        snapshots.append(dict(index=index,graph=record(target / "graph.npz"),
                              partition=record(target / "partition.txt")))
    replay.execute_leiden = observed
    inference = directory / "replay"
    rc = replay.main(["--accuracy-checkpoint",str(checkpoint),"--checkpoint-sha256",
        record(checkpoint / "manifest.json")["sha256"],"--fasta-directory",str(inputs),
        "--output-directory",str(inference),"--json",str(directory / "replay.json"),
        "--cpu",str(cpu),"--matrix","BLOSUM62","--cpm-resolution","0.1",
        "--leiden-seed","4","--profile-iterations","1","--profile-min-species","1"])
    if rc != 0 or len(snapshots) != 4:
        raise ValueError("Replay did not complete exactly four clustering stages")
    names,species,q,t,s = load_accuracy_checkpoint(checkpoint)
    shutil.copyfile(inference / "orthogroups_profiles_refined.txt",
                    inference / "orthohmm_working_res/orthohmm_edges_clustered.txt")
    candidate = _expand_phylogeny_candidates(str(inference),names,species,(q,t,s),profile="satellite_v2")
    candidate.pop("_membership_constraints",None)
    write_json(directory / "stage_report.json",dict(status="prephylogeny_dependency_replay_complete",
        runtime=runtime,snapshots=snapshots,candidates=candidate,
        candidate_partition=record(inference / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"),
        accuracy_evaluated=False,phylogeny_run=False))


def prepare(repo, directory):
    sys.path.insert(0,str(repo))
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    if directory.exists():
        raise FileExistsError(directory)
    path = repo / "benchmark_tools/results/ob_leiden_overlay_probe_20260926.json"
    if record(path)["sha256"] != "3cc5b1795d1d0b361a02e4d18fc98e0a1bc5307ab6c2c086fe28e7bc0e0c4e4a":
        raise ValueError("Changed distribution evidence")
    overlay_report = json.loads(path.read_text())
    # Locate distribution root without assuming which file sorts first.
    overlay = Path(overlay_report["plan"]["path"]).parent / "distribution"
    run = repo / "benchmarks/work/publication_installed_orthobench_20260926"
    plan_path = run / "plan.json"
    if record(plan_path)["sha256"] != "5fd8dc70c337951706fe9246c2bf08d1191196939da066c5b74bba0bd3dd2b12":
        raise ValueError("Changed installed input/source plan")
    installed_plan = json.loads(plan_path.read_text())
    checkpoint = run / "inference/orthohmm_working_res/high_sensitivity_checkpoint"
    records = [record(path),record(plan_path),*installed_plan["checked_records"],
               *[r["copy"] for r in overlay_report["copied"]],
               *[record(p) for p in sorted(checkpoint.iterdir())],record(__file__),
               record(repo / "benchmark_tools/replay_high_sensitivity.py")]
    for item in records:
        check(item)
    site = repo / "benchmarks/work/publication_frozen_overlay_20260926/venv_clean/lib/python3.10/site-packages"
    directory.mkdir(parents=True)
    arms = []
    for label in ("leiden012","leiden011"):
        arm = directory / label
        environment = directory / (label + "_venv")
        venv.EnvBuilder(with_pip=False).create(environment)
        paths = ([str(overlay)] if label == "leiden011" else []) + [str(site)]
        pth = environment / "lib/python3.10/site-packages/frozen_diagnostic.pth"
        pth.write_text("\n".join(paths) + "\n")
        arm.mkdir()
        arms.append(dict(label=label,directory=str(arm),python=str(environment / "bin/python"),pth=record(pth)))
    write_json(directory / "plan.json",dict(arms=arms,checked_records=records,repo=str(repo),
        checkpoint=str(checkpoint),inputs=str(run / "input"),cpu=32,timeout_seconds=21600,
        attempts_per_arm=1,source=record(__file__),scoring=False,phylogeny=False))


def load_plan(plan_path, expected_sha):
    raw = plan_path.read_bytes()
    if not expected_sha or hashlib.sha256(raw).hexdigest() != expected_sha:
        raise ValueError("Changed or unpinned execution plan")
    return json.loads(raw)


def execute(plan_path, expected_sha):
    plan = load_plan(plan_path, expected_sha)
    sys.path.insert(0,plan["repo"])
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    for item in plan["checked_records"]:
        check(item)
    env = {"PATH":"/usr/bin:/bin","HOME":os.environ["HOME"],"LANG":"C.UTF-8",
           "OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1","MKL_NUM_THREADS":"1","PYTHONHASHSEED":"0"}
    for arm in plan["arms"]:
        directory = Path(arm["directory"])
        if (directory / "native.log").exists():
            raise FileExistsError("Existing arm attempt; no implicit resume")
        check(arm["pth"])
        command = [arm["python"],"-I",plan["source"]["path"],"--driver",str(directory),
                   "--repo",plan["repo"],"--checkpoint",plan["checkpoint"],"--inputs",plan["inputs"],
                   "--cpu",str(plan["cpu"])]
        start,timeout = time.monotonic(),False
        with (directory / "native.log").open("x") as log:
            p = subprocess.Popen(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                rc = p.wait(timeout=plan["timeout_seconds"])
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGKILL)
                rc,timeout = p.wait(),True
        write_json(directory / "execution.json",dict(command=command,env=env,returncode=rc,
            timed_out=timeout,attempts=1,job_id=os.environ.get("SLURM_JOB_ID"),
            wall_seconds=time.monotonic()-start,log=record(directory / "native.log")))
        if rc or timeout:
            raise RuntimeError("Arm failed; no retry or next-arm execution")
        runtime = json.loads((directory / "runtime.json").read_text())
        expected = "0.11.0" if arm["label"] == "leiden011" else "0.12.0"
        if runtime["versions"]["leidenalg"] != expected:
            raise ValueError("Wrong Leiden distribution")
    for item in plan["checked_records"]:
        check(item)
    write_json(plan_path.parent / "execution.json",dict(status="native_complete_pending_readback",plan=record(plan_path)))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path)
    p.add_argument("--prepare",type=Path)
    p.add_argument("--execute",type=Path)
    p.add_argument("--plan-sha256")
    p.add_argument("--driver",type=Path)
    p.add_argument("--checkpoint",type=Path)
    p.add_argument("--inputs",type=Path)
    p.add_argument("--cpu",type=int,default=32)
    a = p.parse_args()
    if a.prepare:
        prepare(a.repo.resolve(),a.prepare.resolve())
    elif a.execute:
        execute(a.execute.resolve(),a.plan_sha256)
    elif a.driver:
        driver(a.repo.resolve(),a.driver.resolve(),a.checkpoint.resolve(),a.inputs.resolve(),a.cpu)
    else:
        p.error("Require --prepare, --execute or --driver")
