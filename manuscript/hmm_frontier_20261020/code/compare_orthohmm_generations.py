"""Fresh small-dataset generation comparison, with an immutable execution plan."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import signal
from statistics import median
import subprocess
import sys
import tarfile
import time

import psutil


METHODS = ("original", "sensitive", "phylogeny")
LABELS = {"original": "OrthoHMM1 (v0.1.0)", "sensitive": "OrthoHMM2 sensitive",
          "phylogeny": "OrthoHMM2 + phylogeny"}
ORIGINAL = "7d1540c624ab19a9d6083d9bcbdfe55c7fcd26c2"
KERNELS = ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so")


def now():
    return datetime.now(timezone.utc).isoformat()


def record(path):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=str(path.resolve()), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temporary, path)


def sequences(directory):
    result = []
    for path in sorted(Path(directory).glob("*.faa")):
        ids = [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]
        if not ids or len(ids) != len(set(ids)) or max(map(len, ids)) > 50 or len(path.name) > 50:
            raise ValueError("FASTA must satisfy original implementation identifier limits")
        result.append((path, ids))
    all_ids = [name for _, ids in result for name in ids]
    if len(all_ids) != len(set(all_ids)):
        raise ValueError("Protein identifiers must be globally unique")
    return result


def partition(path, expected):
    groups = [line.split()[1:] for line in Path(path).read_text().splitlines() if line.strip()]
    values = [name for group in groups for name in group]
    if not groups or any(not group for group in groups) or len(values) != len(set(values)) or set(values) != set(expected):
        raise ValueError("Completed partition must cover every input protein exactly once")
    return dict(orthogroups=len(groups), proteins=len(values), sha256=record(path)["sha256"])


def search_complete(directory, files):
    for left in files:
        for right in files:
            path = Path(directory) / f"{left}_2_{right}.phmmerout.txt"
            if not path.is_file() or path.read_text().rstrip().splitlines()[-1] != "# [ok]":
                raise ValueError("Original phmmer search is incomplete")


def schedule(repeats=3):
    if repeats != 3:
        raise ValueError("This balanced technical-repeat design requires three repeats")
    return [dict(index=3*r+j, repeat=r+1, method=METHODS[(r+j) % 3])
            for r in range(3) for j in range(3)]


def snapshot(repo, revision, destination):
    destination.mkdir()
    data = subprocess.check_output(["git", "archive", revision, "orthohmm"], cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        for member in archive.getmembers():
            target = destination / member.name
            if not target.resolve().is_relative_to(destination.resolve()) or member.issym() or member.islnk():
                raise ValueError("Unsafe source archive member")
        archive.extractall(destination)


def prepare(repo, output, python, *, mcl=None, modern_revision="HEAD"):
    repo, output = repo.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name in ("input", "sources", "runs"):
        (output / name).mkdir()
    selected = sorted(sequences(repo / "bacterial_scaling/subsets/n005"), key=lambda v:(len(v[1]), v[0].name))[:3]
    for path, _ in selected:
        shutil.copyfile(path, output / "input" / path.name)
    current = subprocess.check_output(["git", "rev-parse", modern_revision], cwd=repo, text=True).strip()
    snapshot(repo, ORIGINAL, output / "sources/original")
    snapshot(repo, current, output / "sources/current")
    for name in KERNELS:
        shutil.copyfile(repo / "orthohmm/search/csrc" / name, output / "sources/current/orthohmm/search/csrc" / name)
    shutil.copyfile(repo / "benchmark_tools/time_original_orthohmm.py", output / "sources/observe_original.py")
    executables = {name: shutil.which(name) for name in ("phmmer", "mcl", "mafft", "FastTree", "taskset")}
    if mcl is not None:
        executables["mcl"] = str(mcl.resolve())
    if not all(executables.values()):
        raise ValueError("Required external executable is unavailable")
    mcl_help = subprocess.run([executables["mcl"], "-h"], capture_output=True, text=True, check=False)
    if "--abc" not in mcl_help.stdout + mcl_help.stderr:
        raise ValueError("The original CLI requires an MCL build supporting --abc; pass --mcl explicitly")
    cpus = list(range(64, 80))
    if not set(cpus).issubset(os.sched_getaffinity(0)):
        raise ValueError("Planned sixteen-core CPU set is unavailable")
    environment = {name:"1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")}
    environment.update(PYTHONUNBUFFERED="1", PYTHONHASHSEED="0", CUDA_VISIBLE_DEVICES="", LC_ALL="C")
    plan = dict(schema_version=1, created_at_utc=now(), original_revision=ORIGINAL, modern_revision=current,
        cpu_set=cpus, cpu_budget=16, original_cli_cpu=4, modern_cli_cpu=16, modern_threads_per_worker=4,
        python=record(python), executables={k:record(v) for k,v in executables.items()},
        mcl_version=subprocess.check_output([executables["mcl"], "--version"], text=True),
        environment=environment, platform=platform.platform(),
        lscpu=subprocess.check_output(["lscpu", "-p=CPU,CORE,SOCKET,NODE"], text=True),
        source_selection="Three smallest complete proteomes by protein count from the existing n005 set, chosen before timing.",
        input=[dict(record(path), proteins=len(ids)) for path,ids in sequences(output / "input")],
        source_input=[record(path) for path,_ in selected],
        sources=[record(p) for p in sorted((output / "sources").rglob("*")) if p.is_file()],
        harness=record(Path(__file__)), runs=schedule(), per_run_timeout_seconds=1200,
        scope="Sequential fresh end-to-end orthogroup inference to --stop infer; sequence FASTA exports excluded. Original uses MCL, modern uses high-sensitivity Leiden with and without inferred phylogeny. Different methods/outputs, not an equal-accuracy comparison.",
        limitations=["Shared Threadripper host; no process isolation, exclusive CPUs or contention correction. OS affinity limits eligible CPUs only.",
            "Original --cpu 4 creates up to four phmmer processes each with four compute threads. Modern --cpu 16 and four threads per worker use the same maximum sixteen eligible physical cores. Original graph/MCL phases retain its four-worker setting.",
            "Three technical repeats on one small bacterial dataset do not establish general scalability or accuracy. Native process startup/import costs are included, and OS caches are not flushed.",
            "This frozen v0.1.0 tag is the explicitly identified original implementation, not a claimed version 1.0 release."])
    save(output / "plan.json", plan)
    return plan


def command(output, run, plan):
    directory = output / "runs" / f'{run["index"]:02d}_{run["method"]}'
    metrics = directory / "metrics.json"
    common = [str(output / "input"), "-o", str(directory / "output"), "-x", "BLOSUM62", "-e", "0.0001", "--stop", "infer"]
    if run["method"] == "original":
        source = output / "sources/original"
        args = [plan["python"]["path"], str(output / "sources/observe_original.py"), "--source", str(source),
                "--metrics", str(metrics), "--", *common, "-c", "4", "-p", plan["executables"]["phmmer"]["path"],
                "-m", plan["executables"]["mcl"]["path"], "-i", "1.5"]
    else:
        source = output / "sources/current"
        args = [plan["python"]["path"], "-m", "orthohmm", *common, "-c", "16", "--threads_per_worker", "4",
                "--accuracy_profile", "high_sensitivity", "--search_result_order", "canonical", "--clustering", "leiden",
                "--cpm_resolution", "0.1", "--metrics_json", str(metrics)]
        if run["method"] == "phylogeny":
            args += ["--phylogeny", "reconcile", "--species_tree_mode", "infer", "--phylogeny_candidates", "satellite_v2",
                     "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
                     "--species_tree_rooting", "min_variance", "--aligner", plan["executables"]["mafft"]["path"],
                     "--tree_builder", plan["executables"]["FastTree"]["path"]]
    return directory, source, [plan["executables"]["taskset"]["path"], "-c", ",".join(map(str,plan["cpu_set"])), *args]


def host():
    return dict(cpu_times=psutil.cpu_times()._asdict(), load_average=list(os.getloadavg()),
                available_memory_bytes=psutil.virtual_memory().available)


def run_panel(output):
    output = output.resolve()
    plan = json.loads((output / "plan.json").read_text())
    for ref in [plan["python"], *plan["executables"].values(), *plan["sources"], *plan["input"]]:
        if record(ref["path"]) != {k:ref[k] for k in ("path", "bytes", "sha256")}:
            raise ValueError("Planned input or runtime bytes changed")
    expected = [name for _,ids in sequences(output / "input") for name in ids]
    results = []
    for run in plan["runs"]:
        directory, source, args = command(output, run, plan)
        directory.mkdir(exist_ok=False)
        (directory / "output").mkdir()
        environment = os.environ.copy()
        environment.update(plan["environment"], PYTHONPATH=str(source))
        row = dict(run, command=args, cwd=str(source), started_at_utc=now(), host_before=host(), status="running")
        save(directory / "attempt.json", row)
        print(json.dumps({"started":run}), flush=True)
        started = time.perf_counter()
        with (directory / "native.log").open("x") as stream:
            child = subprocess.Popen(args, cwd=source, env=environment, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                row["returncode"] = child.wait(timeout=plan["per_run_timeout_seconds"])
                row["elapsed_s"] = time.perf_counter() - started
                row["status"] = "completed" if row["returncode"] == 0 else "failed"
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                row.update(status="timed_out", elapsed_s=time.perf_counter() - started)
        row.update(finished_at_utc=now(), host_after=host())
        try:
            if row["status"] != "completed":
                raise ValueError("Native execution did not complete")
            metrics = json.loads((directory / "metrics.json").read_text())
            if metrics["status"] not in {"complete", "completed"}:
                raise ValueError("Native metrics are incomplete")
            row["search_s"] = float(metrics["stages"]["search"]["wall_s"])
            row["partition"] = partition(directory / "output/orthohmm_orthogroups.txt", expected)
            if run["method"] == "original":
                search_complete(directory / "output/orthohmm_working_res", [Path(r["path"]).name for r in plan["input"]])
            row["status"] = "validated"
        except (ValueError, KeyError, FileNotFoundError, json.JSONDecodeError) as error:
            row.update(validation_error=str(error), status=row["status"] if row["status"] != "completed" else "invalid_outputs")
        row["artifacts"] = [record(p) for p in (directory / "metrics.json", directory / "native.log") if p.exists()]
        save(directory / "attempt.json", row)
        results.append(row)
        save(output / "results.json", dict(plan=record(output / "plan.json"), runs=results, status="running" if len(results)<9 else "finished"))
        print(json.dumps({"finished":{k:row[k] for k in ("index", "method", "status", "elapsed_s")}}), flush=True)


def summarize(data):
    runs = data["runs"]
    if len(runs) != 9 or [{k:r[k] for k in ("index", "repeat", "method")} for r in runs] != schedule():
        raise ValueError("Require all nine prespecified attempts in order")
    rows = []
    for method in METHODS:
        values = [r for r in runs if r["method"] == method]
        if len(values) != 3 or any(r["status"] != "validated" for r in values):
            raise ValueError("Incomplete cells cannot be summarized as successful speed comparisons")
        row = dict(method=method, label=LABELS[method])
        for key in ("elapsed_s", "search_s"):
            scores = [float(r[key]) for r in values]
            if any(not math.isfinite(v) or v <= 0 for v in scores):
                raise ValueError("Runtime must be finite and positive")
            row[key] = dict(median=median(scores), minimum=min(scores), maximum=max(scores), values=scores)
        rows.append(row)
    for row in rows:
        row["original_to_modern_elapsed_ratio"] = rows[0]["elapsed_s"]["median"] / row["elapsed_s"]["median"]
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "summarize"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--mcl", type=Path)
    parser.add_argument("--modern-revision", default="HEAD")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.repo, args.output, args.python, mcl=args.mcl, modern_revision=args.modern_revision)
    elif args.action == "run":
        run_panel(args.output)
    else:
        summary = summarize(json.loads((args.output / "results.json").read_text()))
        save(args.output / "summary.json", summary)
        print(json.dumps(summary, indent=2))
