"""Execute one frozen recovery query batch with durable, non-admitting status."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_replay_panel import PANEL_SHA
from benchmark_tools.run_qfo_corrected_blast import verify, environment
from benchmark_tools.run_simulation_methods import read_frozen

BATCHES_SHA = "5e59e7deedb71feb6453057f0b379bda0e06a3556651b13eb7275b26155f45e0"
PROTOCOL_SHA = "4f120db4cb0733fc8d128461465bab877d91d31f668c03b45a14a05be3bb9410"


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def save_status(path, report):
    temporary = path.with_suffix(".tmp")
    with temporary.open("w") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    sync_directory(path.parent)


def command_for(original, query, output):
    command = list(original)
    for flag in ("-i", "-o", "-d", "-a"):
        if command.count(flag) != 1 or command.index(flag) == len(command)-1:
            raise ValueError("Ambiguous legacy search command")
    if command[command.index("-a")+1] != "180":
        raise ValueError("Changed frozen thread count")
    command[command.index("-i")+1] = str(query)
    command[command.index("-o")+1] = str(output)
    return command


def prepare(root, index):
    if type(index) is not int or index not in range(20):
        raise ValueError("Unknown recovery batch")
    results = root / "benchmark_tools/results"
    batches_path = results / "qfo_blast_recovery_batches_20260923.json"
    batches = read_frozen(batches_path, BATCHES_SHA)
    panel_path = results / "qfo_blast_replay_panel_20260923.json"
    panel = read_frozen(panel_path, PANEL_SHA)
    protocol = record(results / "QFO_BLAST_RECOVERY_PROTOCOL_20260923.md")
    if protocol["sha256"] != PROTOCOL_SHA or batches["total_queries"] != 98913 or len(batches["batches"]) != 20:
        raise ValueError("Changed recovery plan")
    plan_path, runtime_path = (Path(panel[k]["path"]) for k in ("plan", "runtime"))
    plan = verify(plan_path, runtime_path)
    batch = batches["batches"][index]
    if batch["index"] != index or batch["queries"] != (5000 if index < 19 else 3913):
        raise ValueError("Wrong batch identity")
    directory = root / "benchmarks/results/qfo_blast_recovery_v1" / f"batch_{index:02d}"
    checked = [record(batches_path), record(panel_path), protocol, record(__file__),
        *batches["inputs"], batch["input"], panel["plan"], panel["runtime"], *panel["database"],
        record("/usr/bin/time"), record(Path(__file__).with_name("run_qfo_corrected_blast.py"))]
    for item in checked:
        check(item)
    return dict(index=index, batch=batch, directory=str(directory),
        command=command_for(plan["search_commands"]["blast"], batch["input"]["path"], directory / "hits.blast.partial"),
        cwd=str(Path(plan["output_root"]) / "work"), environment=environment(),
        checked_records=checked, plan_path=str(plan_path), runtime_path=str(runtime_path))


def execute(preflight, job_id, array_job_id):
    directory = Path(preflight["directory"])
    directory.parent.mkdir(parents=True, exist_ok=True)
    directory.mkdir(exist_ok=False)
    sync_directory(directory.parent)
    status = directory / "status.json"
    report = dict(status="starting", preflight=preflight, job_id=job_id, array_job_id=array_job_id,
        index=preflight["index"], node=os.uname().nodename, started_epoch=time.time(),
        search_admitted=False, reuse_authorized=False, publication_ready=False)
    save_status(directory / "preflight.json", report)
    save_status(status, report)
    try:
        log, timing = directory / "blast.log", directory / "blast.time.txt"
        with log.open("xb") as stream:
            result = subprocess.run(["/usr/bin/time", "-v", "-o", str(timing), *preflight["command"]],
                cwd=preflight["cwd"], env=preflight["environment"], stdout=stream, stderr=subprocess.STDOUT)
        report.update(exit_code=result.returncode, finished_epoch=time.time())
        partial = directory / "hits.blast.partial"
        for path in (log, timing, partial):
            if path.is_file():
                with path.open("rb") as stream:
                    os.fsync(stream.fileno())
        sync_directory(directory)
        report["outputs"] = [record(p) for p in (log, timing, partial) if p.is_file()]
        if result.returncode:
            raise RuntimeError(f"Native recovery batch exited {result.returncode}")
        if not all(p.is_file() for p in (log, timing, partial)):
            raise ValueError("Missing completed native artifacts")
        for item in preflight["checked_records"]:
            check(item)
        verify(Path(preflight["plan_path"]), Path(preflight["runtime_path"]))
        # Native bytes reach stable storage before the completed status marker.
        final = directory / "hits.blast"
        if final.exists() or final.is_symlink():
            raise FileExistsError(final)
        partial.rename(final)
        sync_directory(directory)
        report.update(status="native_batch_completed_pending_admission", outputs=[record(p) for p in (log, timing, final)])
        save_status(status, report)
    except BaseException as error:
        report.update(status="failed_or_interrupted", error_type=type(error).__name__, error=str(error))
        save_status(status, report)
        raise
    return report


def run(root, index):
    if (os.environ.get("SLURM_CPUS_PER_TASK") != "180" or not os.environ.get("SLURM_JOB_ID")
            or not os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled 180-CPU matching array task")
    return execute(prepare(root.resolve(), index), os.environ["SLURM_JOB_ID"], os.environ["SLURM_ARRAY_JOB_ID"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(20))
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only:
        prepare(args.root.resolve(), args.index)
    else:
        run(args.root, args.index)
