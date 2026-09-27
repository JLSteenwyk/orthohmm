"""Isolate the retained leidenalg 0.11 distribution in the clean runtime."""

import argparse
import base64
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_installed_ob_clustering import write_json
from benchmark_tools.compare_installed_ob_search import partition
from benchmark_tools.audit_installed_orthobench import compare_partitions

PRIOR_SHA = "0cd3e0f86ed1033af1cb65f4759fce38052d1d2a996b1ce5c7e7c3f9f6559fff"


def allowed_distribution_path(path):
    path = Path(path)
    if path.is_absolute() or ".." in path.parts or not path.parts or path.parts[0] not in {
            "leidenalg", "leidenalg.libs", "leidenalg-0.11.0.dist-info"}:
        raise ValueError("Unexpected distribution path")
    return "__pycache__" not in path.parts and path.suffix != ".pyc"


def stage_distribution(destination):
    distribution = importlib.metadata.distribution("leidenalg")
    if distribution.version != "0.11.0":
        raise ValueError("Expected retained 0.11.0 distribution")
    destination.mkdir()
    copied = []
    for entry in distribution.files:
        if not allowed_distribution_path(entry):
            continue
        source = Path(distribution.locate_file(entry))
        if source.is_symlink() or not source.is_file():
            raise ValueError("Invalid source distribution file")
        data = source.read_bytes()
        if entry.hash:
            if entry.hash.mode != "sha256":
                raise ValueError("Unsupported RECORD digest")
            digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip("=")
            if digest != entry.hash.value:
                raise ValueError("Source differs from installed RECORD")
        target = destination / entry
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
        copied.append(dict(source=record(source), copy=record(target)))
    return copied


def run(repo, directory):
    if directory.exists():
        raise FileExistsError(directory)
    prior_path = repo / "benchmark_tools/results/ob_clustering_runtime_probe_20260926.json"
    if record(prior_path)["sha256"] != PRIOR_SHA:
        raise ValueError("Changed runtime contrast")
    prior = json.loads(prior_path.read_text())
    for item in prior["checked_records"]:
        check(item)
    for item in prior["installed_runtime"]["files"]:
        check(item)
    check(prior["arms"][0]["prediction"])
    installed_item = next(r for r in prior["checked_records"]
                          if Path(r["path"]).name == "installed_ob_clustering_probe_20260926.json")
    installed = json.loads(Path(installed_item["path"]).read_text())
    python = installed["runtime"]["executable"]
    env = installed["settings"]["env"]
    source_payload = Path(installed["execution_details"][0]["command"][-1])
    names = (source_payload / "gene_names.txt").read_text().splitlines()
    expected = partition(Path(prior["arms"][0]["prediction"]["path"]), set(names))
    baseline = partition(Path(installed["arms"][0]["prediction"]["path"]), set(names))
    directory.mkdir(parents=True)
    overlay = directory / "distribution"
    copied = stage_distribution(overlay)
    old_extension = next(r for r in prior["runtime"]["files"] if "_c_leiden" in r["path"])
    copied_extension = next(r for r in copied if "_c_leiden" in r["source"]["path"])
    if copied_extension["source"] != old_extension:
        raise ValueError("Overlay extension differs from prior runtime")
    write_json(directory / "plan.json", dict(source=record(__file__),prior=record(prior_path),
        copied=copied, python=python, env=env, arms=["overlay", "overlay_repeat"],
        attempts_per_arm=1, timeout_seconds=1800, cpm_resolution=0.1,seed=4,include_isolates=True))
    # Capture imports in the same child process that invokes the frozen worker.
    code = (
        "import sys,json,importlib.metadata as m; from pathlib import Path; "
        "overlay=sys.argv.pop(1); evidence=sys.argv.pop(1); sys.path.insert(0,overlay); "
        "import numpy,igraph._igraph,leidenalg._c_leiden; "
        "from orthohmm import externals,helpers,leiden_worker; "
        "r=dict(python=sys.version,executable=sys.executable,"
        "versions={n:m.version(n) for n in ('numpy','igraph','leidenalg')},"
        "files=[x.__file__ for x in (numpy,externals,helpers,leiden_worker,igraph._igraph,leidenalg._c_leiden)],"
        "maps=Path('/proc/self/maps').read_text()); "
        "Path(evidence).write_text(json.dumps(r,indent=2)); leiden_worker.main()"
    )
    arms, groups = [], []
    for label in ("overlay", "overlay_repeat"):
        arm = directory / label
        payload = arm / "payload"
        payload.mkdir(parents=True)
        (arm / "orthohmm_working_res").mkdir()
        for item in installed["execution_details"][0]["inputs"]:
            check(item)
            source = Path(item["path"])
            if source.name != "metadata.json":
                (payload / source.name).write_bytes(source.read_bytes())
        write_json(payload / "metadata.json", dict(cpm_resolution=0.1,seed=4,
                   include_isolates=True,output_directory=str(arm)))
        inputs = [record(p) for p in sorted(payload.iterdir())]
        command = [python,"-I","-B","-c",code,str(overlay),str(arm / "runtime.json"),str(payload)]
        started, timeout = time.monotonic(), False
        with (arm / "native.log").open("x") as stream:
            process = subprocess.Popen(command,cwd=arm,env=env,stdout=stream,
                                       stderr=subprocess.STDOUT,start_new_session=True)
            try:
                rc = process.wait(timeout=1800)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL)
                rc, timeout = process.wait(), True
        execution = dict(command=command,returncode=rc,timed_out=timeout,attempts=1,
                         wall_seconds=time.monotonic()-started,inputs=inputs,log=record(arm / "native.log"))
        write_json(arm / "execution.json",execution)
        if rc or timeout:
            raise RuntimeError(f"{label} failed; no retry")
        runtime = json.loads((arm / "runtime.json").read_text())
        runtime["files"] = [record(p) for p in runtime["files"]]
        unchanged = {Path(r["path"]).name:r for r in installed["runtime"]["files"]
                     if "_c_leiden" not in r["path"]}
        for item in runtime["files"]:
            if Path(item["path"]).name in unchanged and item != unchanged[Path(item["path"]).name]:
                raise ValueError("Unexpected scientific/igraph module substitution")
        if runtime["versions"] != {"numpy":"2.2.6","igraph":"1.0.0","leidenalg":"0.11.0"}:
            raise ValueError("Unexpected overlay versions")
        runtime["mapping_record"] = record(arm / "runtime.json")
        runtime["mapped_leiden_libraries"] = sorted({line.split()[-1] for line in runtime.pop("maps").splitlines()
                                                    if "leidenalg" in line and "/" in line})
        out = arm / "orthohmm_working_res/orthohmm_edges_clustered.txt"
        group = partition(out,set(names))
        groups.append(group)
        arms.append(dict(arm=label,runtime=runtime,execution=execution,prediction=record(out),
            versus_installed=compare_partitions(baseline,group),
            versus_current_anaconda=compare_partitions(expected,group)))
        for item in inputs:
            check(item)
    for item in copied:
        check(item["source"])
        check(item["copy"])
    for item in prior["installed_runtime"]["files"]:
        check(item)
    write_json(directory / "report.json", dict(status="private_leiden_distribution_contrast_complete",
        source=record(__file__),plan=record(directory / "plan.json"),prior=record(prior_path),
        copied=copied,arms=arms,repeat=compare_partitions(*groups),historical_scores_replaced=False,
        limitations=["Private distribution overlay, not a replacement installation or recommended dependency downgrade.",
            "Varies the complete leidenalg artifact including bundled libigraph/libleidenalg, not Python version text alone.",
            "Conditional on this graph, seed, indexing and clean runtime; two planned calls without retries.",
            "No downstream F1 attribution, historical environment attestation or general determinism claim."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    run(args.repo.resolve(),args.output.resolve())
