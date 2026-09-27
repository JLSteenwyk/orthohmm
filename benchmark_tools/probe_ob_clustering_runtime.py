"""Compare current replay-launcher runtime to retained installed clustering."""

import argparse
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

PRIOR_SHA = "2fd69ad44dfe5ff5d1bb6f4d82b31523d6cce06d2a2e8173ad036ef68f527f13"


def verify_sources(observed, expected):
    left = {Path(r["path"]).name: r["sha256"] for r in observed}
    right = {Path(r["path"]).name: r["sha256"] for r in expected}
    names = {"externals.py", "helpers.py", "leiden_worker.py"}
    if not names <= left.keys() or any(left[n] != right[n] for n in names):
        raise ValueError("Scientific clustering sources differ")


def run(repo, directory):
    if directory.exists():
        raise FileExistsError(directory)
    prior_path = repo / "benchmark_tools/results/installed_ob_clustering_probe_20260926.json"
    if record(prior_path)["sha256"] != PRIOR_SHA:
        raise ValueError("Changed installed clustering receipt")
    prior = json.loads(prior_path.read_text())
    records = [record(prior_path), prior["plan"], prior["arms"][0]["prediction"],
               *prior["execution_details"][0]["inputs"]]
    for item in records:
        check(item)
    source_payload = Path(prior["execution_details"][0]["command"][-1])
    names = (source_payload / "gene_names.txt").read_text().splitlines()
    original = partition(Path(prior["arms"][0]["prediction"]["path"]), set(names))
    core = repo / "benchmarks/work/publication_method_native_v2"
    python = "/home/bizon/anaconda3/bin/python"
    env = prior["settings"]["env"]
    prefix = "import sys; sys.path.insert(0,sys.argv.pop(1)); "
    metadata = prefix + (
        "import json,importlib.metadata as m; import igraph._igraph,leidenalg._c_leiden; "
        "from orthohmm import externals,helpers,leiden_worker; "
        "print(json.dumps(dict(python=sys.version,executable=sys.executable,"
        "versions={n:m.version(n) for n in ('numpy','igraph','leidenalg')},"
        "files=[x.__file__ for x in (externals,helpers,leiden_worker,igraph._igraph,leidenalg._c_leiden)])))"
    )
    runtime = json.loads(subprocess.check_output([python,"-I","-c",metadata,str(core)],
                                                env=env,text=True,timeout=60))
    runtime["files"] = [record(p) for p in runtime["files"]]
    verify_sources(runtime["files"], prior["runtime"]["files"])
    records.extend(runtime["files"])
    directory.mkdir(parents=True)
    write_json(directory / "plan.json", dict(arms=["current_launcher", "current_launcher_repeat"],
        attempts_per_arm=1, timeout_seconds=1800, source=record(__file__), checked_records=records,
        runtime=runtime, installed_runtime=prior["runtime"], env=env,
        cpm_resolution=0.1, seed=4, include_isolates=True))
    results, partitions = [], []
    for label in ("current_launcher", "current_launcher_repeat"):
        arm = directory / label
        payload = arm / "payload"
        payload.mkdir(parents=True)
        (arm / "orthohmm_working_res").mkdir()
        for item in prior["execution_details"][0]["inputs"]:
            path = Path(item["path"])
            if path.name != "metadata.json":
                (payload / path.name).write_bytes(path.read_bytes())
        write_json(payload / "metadata.json", dict(cpm_resolution=0.1,seed=4,
                   include_isolates=True,output_directory=str(arm)))
        inputs = [record(p) for p in sorted(payload.iterdir())]
        code = prefix + "import runpy; runpy.run_module('orthohmm.leiden_worker',run_name='__main__')"
        command = [python,"-I","-c",code,str(core),str(payload)]
        start, timeout = time.monotonic(), False
        with (arm / "native.log").open("x") as stream:
            process = subprocess.Popen(command,cwd=arm,env=env,stdout=stream,
                                       stderr=subprocess.STDOUT,start_new_session=True)
            try:
                exitcode = process.wait(timeout=1800)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL)
                exitcode, timeout = process.wait(), True
        execution = dict(command=command,returncode=exitcode,timed_out=timeout,attempts=1,
                         inputs=inputs,log=record(arm / "native.log"),wall_seconds=time.monotonic()-start)
        write_json(arm / "execution.json", execution)
        if exitcode or timeout:
            raise RuntimeError(f"{label} failed; no retry")
        for item in inputs:
            check(item)
        out = arm / "orthohmm_working_res/orthohmm_edges_clustered.txt"
        groups = partition(out,set(names))
        partitions.append(groups)
        results.append(dict(arm=label,execution=execution,prediction=record(out),
                            versus_installed=compare_partitions(original,groups)))
    for item in records:
        check(item)
    write_json(directory / "report.json", dict(status="current_launcher_runtime_contrast_complete",
        source=record(__file__),plan=record(directory / "plan.json"),runtime=runtime,
        installed_runtime=prior["runtime"],checked_records=records,arms=results,
        repeat=compare_partitions(*partitions),historical_runtime_identity_established=False,
        accuracy_evaluated=False,historical_scores_replaced=False,
        limitations=["Current Anaconda environment at the recorded historical launcher path, not attested historical dependency identity.",
            "Same graph, weights, ordered gene/edge indices, seed, CPM setting and scientific worker sources.",
            "Whole-runtime comparison, not an isolated causal test of the leidenalg version.",
            "Initial clustering only; no final-score or whole-pipeline causal attribution.",
            "Two planned calls and one retained installed result; shared-host times are descriptive only."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    run(args.repo.resolve(),args.output.resolve())
