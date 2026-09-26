"""Execute one frozen search-sensitivity dataset once, retaining failed attempts."""

import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.search_sensitivity_hmm import input_inventory


PLAN_SHA = "b3b19927821f46f0e79fcc444987a51b88f699e7665bd364fa6639a4239323bc"
INSTALL_SHA = "467036ee88d5e6f15d322aea4c94397591ef7ba521ac9f8ddde4b7cb783b97d2"


def diamond_commands(binary, queries, target, destination):
    database = destination / (target.stem + ".dmnd")
    hits = destination / (target.stem + ".tsv")
    return [
        [str(binary), "makedb", "--in", str(target), "--db", str(database), "--threads", "4"],
        [str(binary), "blastp", "--query", str(queries), "--db", str(database),
         "--out", str(hits), "--outfmt", "6", "qseqid", "sseqid", "qlen", "slen",
         "score", "bitscore", "evalue", "--very-sensitive", "--matrix", "BLOSUM62",
         "--gapopen", "11", "--gapextend", "1", "--comp-based-stats", "1",
         "--masking", "1", "--max-target-seqs", "0", "--max-hsps", "1",
         "--evalue", "1.0", "--threads", "4"],
    ]


def execute(command, directory, name, env, deadline, stages):
    stage = dict(name=name, command=command, started_unix=time.time())
    stages.append(stage)
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("Cell deadline exceeded")
    timed = ["/usr/bin/time", "-v", "-o", str(directory / (name + ".time.txt")), *command]
    with (directory / (name + ".log")).open("x") as log:
        process = subprocess.Popen(timed, cwd=directory, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=remaining)
            stage["returncode"] = code
            if code:
                raise RuntimeError(f"{name} exited {code}")
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            stage["returncode"] = process.returncode
            stage["ended_unix"] = time.time()


def verify_inputs(row):
    observed, ids = input_inventory(Path(row["input"]))
    if sorted(observed, key=lambda r: r["path"]) != sorted(row["inputs"], key=lambda r: r["path"]):
        raise ValueError("Input inventory differs from frozen plan")
    if len(ids) != row["species"] or sum(map(len, ids.values())) != row["genes"]:
        raise ValueError("Input gene/species count differs from plan")


def run(plan_path, index, output):
    plan_record = record(plan_path)
    if plan_record["sha256"] != PLAN_SHA:
        raise ValueError("Unrecognized frozen plan")
    plan = json.loads(plan_path.read_text())
    if len(plan["datasets"]) != 70 or not 0 <= index < 70:
        raise ValueError("Invalid dataset index/panel size")
    row = plan["datasets"][index]
    output.mkdir(parents=True, exist_ok=False)
    receipt = dict(status="failed", dataset_index=index, condition=row["condition"], seed=row["seed"],
                   split=row["split"], plan=plan_record, stages=[], attempt=1,
                   slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                   slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"),
                   slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                   timing_comparability="shared_host_descriptive_only")
    deadline = time.monotonic() + 3300
    try:
        checked = [plan_record, plan["diamond"]]
        install_record = next(r for r in plan["checked_evidence"] if r["sha256"] == INSTALL_SHA)
        check(install_record)
        installed = json.loads(Path(install_record["path"]).read_text())
        checked.extend([install_record, *installed["checked_records"]])
        staging = next(Path(r["path"]) for r in installed["checked_records"] if Path(r["path"]).name == "staging.json")
        python = staging.parent / "venv_clean/bin/python"
        helper = Path(__file__).with_name("search_sensitivity_hmm.py").resolve()
        checked.extend([record(__file__), record(helper), record(python), record("/usr/bin/time")])
        receipt["checked_records"] = checked
        receipt["python_command_path"] = str(python)
        for item in checked:
            check(item)
        verify_inputs(row)
        private = output / "input"
        private.mkdir()
        for item in row["inputs"]:
            shutil.copyfile(item["path"], private / Path(item["path"]).name)
        private_records, _ = input_inventory(private)
        expected = {(Path(r["path"]).name, r["bytes"], r["sha256"]) for r in row["inputs"]}
        if {(Path(r["path"]).name, r["bytes"], r["sha256"]) for r in private_records} != expected:
            raise ValueError("Private input copy mismatch")
        receipt["private_inputs"] = private_records
        env = dict(PATH="/usr/bin:/bin", HOME=str(output), LC_ALL="C", OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
        receipt["environment"] = env
        execute([str(python), "-I", str(helper), "--input", str(private), "--output", str(output / "hmm")],
                output, "hmm", env, deadline, receipt["stages"])
        queries = output / "queries.fasta"
        with queries.open("xb") as stream:
            for path in sorted(private.iterdir()):
                with path.open("rb") as source:
                    shutil.copyfileobj(source, stream)
                stream.write(b"\n")
        diamond = output / "diamond"
        diamond.mkdir()
        diamond_outputs = []
        for target in sorted(private.iterdir()):
            for number, command in enumerate(diamond_commands(plan["diamond"]["path"], queries, target, diamond)):
                execute(command, output, f"diamond_{target.stem}_{number}", env, deadline, receipt["stages"])
            diamond_outputs.extend([record(diamond / (target.stem + ".dmnd")),
                                    record(diamond / (target.stem + ".tsv"))])
        verify_inputs(row)
        if input_inventory(private)[0] != private_records:
            raise ValueError("Private input changed during search")
        for item in checked:
            check(item)
        receipt["outputs"] = [record(output / "hmm/hits.tsv"), record(output / "hmm/receipt.json"),
                              record(queries), *diamond_outputs]
        receipt["status"] = "native_completed_pending_independent_readback"
    except BaseException as error:
        receipt["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        receipt["logs"] = [record(p) for p in sorted(output.iterdir())
                           if p.is_file() and (p.name.endswith(".log") or p.name.endswith(".time.txt"))]
        with (output / "execution.json").open("x") as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write("\n")


def interrupted(signum, frame):
    raise RuntimeError(f"Received signal {signum}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    run(args.plan.resolve(), args.index, args.output.absolute())
