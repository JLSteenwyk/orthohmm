"""Isolated initial-Leiden comparison of retained historical/fresh weights."""

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.audit_installed_orthobench import compare_partitions, PLAN_SHA
from benchmark_tools.compare_installed_ob_search import partition

GRAPH_SHA = "e930906950cc65a1b1141aa49b8e6b25d8b9d3e7ac15319505170bceb955c27f"
ARMS = (("historical", "historical_order_historical_scores"),
        ("fresh", "historical_order_fresh_scores"),
        ("historical_repeat", "historical_order_historical_scores"))


def validate_graph(names, sources, targets, weights):
    if not names or len(set(names)) != len(names) or any(not n or len(n.split()) != 1 for n in names):
        raise ValueError("Invalid gene names")
    if (any(a.ndim != 1 for a in (sources, targets, weights))
            or not len(sources) == len(targets) == len(weights)
            or sources.dtype.kind not in "iu" or targets.dtype.kind not in "iu"):
        raise ValueError("Invalid graph shape or index dtype")
    if len(sources) and (sources.min() < 0 or targets.max() >= len(names)
                        or np.any(sources >= targets)):
        raise ValueError("Invalid canonical graph endpoints")
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Invalid graph weights")
    keys = sources.astype(np.int64) * len(names) + targets
    if np.any(np.diff(keys) <= 0):
        raise ValueError("Graph endpoints must be unique and sorted")


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def run(repo, directory):
    if directory.exists():
        raise FileExistsError(directory)
    graph_path = repo / "benchmark_tools/results/installed_ob_graph_probe_20260926.json"
    plan_path = repo / "benchmark_tools/results/installed_orthobench_plan_20260926.json"
    if record(graph_path)["sha256"] != GRAPH_SHA or record(plan_path)["sha256"] != PLAN_SHA:
        raise ValueError("Changed graph or installed plan")
    graphs, plan = json.loads(graph_path.read_text()), json.loads(plan_path.read_text())
    records = [record(graph_path), record(plan_path), graphs["gene_names"],
               *graphs["artifacts"].values(), *plan["checked_records"]]
    for item in records:
        check(item)
    python = repo / "benchmarks/work/publication_frozen_overlay_20260926/venv_clean/bin/python"
    env = {"PATH": "/usr/bin:/bin", "HOME": os.environ["HOME"], "LANG": "C.UTF-8",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "PYTHONHASHSEED": "0"}
    metadata_code = (
        "import json,sys,importlib.metadata as m; import numpy,igraph,leidenalg; "
        "import igraph._igraph,leidenalg._c_leiden; "
        "from orthohmm import externals,helpers,leiden_worker; "
        "print(json.dumps(dict(python=sys.version, executable=sys.executable, "
        "versions={n:m.version(n) for n in ('numpy','igraph','leidenalg')}, "
        "files=[x.__file__ for x in (externals,helpers,leiden_worker,igraph._igraph,leidenalg._c_leiden)])))"
    )
    runtime = json.loads(subprocess.check_output([str(python), "-I", "-c", metadata_code],
                                                env=env, text=True, timeout=60))
    runtime["files"] = [record(p) for p in runtime["files"]]
    records.extend(runtime["files"])
    names = Path(graphs["gene_names"]["path"]).read_text().splitlines()
    directory.mkdir(parents=True)
    write_json(directory / "plan.json", dict(arms=ARMS, runtime=runtime, env=env,
        cpm_resolution=0.1, seed=4, include_isolates=True, attempts_per_arm=1,
        per_arm_timeout_seconds=1800, checked_records=records, source=record(__file__)))
    completed, groups = [], {}
    endpoints = None
    for label, graph_key in ARMS:
        arm = directory / label
        payload = arm / "payload"
        payload.mkdir(parents=True)
        working = arm / "orthohmm_working_res"
        working.mkdir()
        with np.load(graphs["artifacts"][graph_key]["path"], allow_pickle=False) as archive:
            arrays = [archive[k] for k in ("sources", "targets", "weights")]
        validate_graph(names, *arrays)
        if endpoints is not None and not all(np.array_equal(a, b) for a, b in zip(endpoints, arrays[:2])):
            raise ValueError("Graph endpoints differ across arms")
        endpoints = arrays[:2]
        for key, array in zip(("sources", "targets", "weights"), arrays):
            np.save(payload / f"{key}.npy", array, allow_pickle=False)
        (payload / "gene_names.txt").write_text("\n".join(names) + "\n")
        write_json(payload / "metadata.json", dict(cpm_resolution=0.1, seed=4,
                    include_isolates=True, output_directory=str(arm)))
        inputs = [record(p) for p in sorted(payload.iterdir())]
        command = [str(python), "-I", "-m", "orthohmm.leiden_worker", str(payload)]
        started = time.monotonic()
        timed_out = False
        with (arm / "native.log").open("x") as log:
            process = subprocess.Popen(command, cwd=arm, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=1800)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                code, timed_out = process.wait(), True
        execution = dict(command=command, returncode=code, timed_out=timed_out,
                         wall_seconds=time.monotonic()-started, inputs=inputs,
                         log=record(arm / "native.log"), attempts=1)
        write_json(arm / "execution.json", execution)
        if code or timed_out:
            raise RuntimeError(f"{label} failed; outputs retained, no retry")
        for item in inputs:
            check(item)
        out = working / "orthohmm_edges_clustered.txt"
        groups[label] = partition(out, set(names))
        completed.append(dict(arm=label, execution=record(arm / "execution.json"),
                              prediction=record(out), groups=len(groups[label])))
    for item in records:
        check(item)
    result = dict(status="installed_initial_clustering_weight_probe_complete", source=record(__file__),
        plan=record(directory / "plan.json"), runtime=runtime, arms=completed,
        historical_vs_fresh=compare_partitions(groups["historical"], groups["fresh"]),
        historical_repeat=compare_partitions(groups["historical"], groups["historical_repeat"]),
        accuracy_evaluated=False, historical_scores_replaced=False,
        limitations=["Fixed canonical gene indexing, seed 4, CPM 0.1 and installed runtime; not full historical execution.",
            "Three planned initial clustering calls, no search/profile/refinement/candidate/phylogeny runs.",
            "One repeat checks this case, not universal determinism; no score tuning or seed selection.",
            "Wall times are shared-host descriptive observations, not controlled efficiency evidence."])
    write_json(directory / "report.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve())
