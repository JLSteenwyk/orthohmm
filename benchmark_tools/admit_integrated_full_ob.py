"""Independently admit only completed full integrated OrthoBench job 22337."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_frozen_overlay_install import local_install_wheels
from benchmark_tools.audit_installed_orthobench import read_root_hogs
from benchmark_tools.audit_recovery_advisories import verify_inventory
from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.readback_canonical_ob_phylogeny import comparison, SCORER_SHA
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.run_integrated_full_job import validate as validate_plan
from benchmark_tools.run_integrated_publication_workflow import install_commands, validate_data
from benchmark_tools.run_publication_pipeline import save
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.validate_reader_upgrade import verify_lock
from benchmark_tools.verify_ygob_validation import require_completed_job

JOB = 22337
PLAN_SHA = "5cab1cd592e6f6d69276c84d0ab5e1b9dc91e1e12a44ace757db87d08df2d899"


def read(path):
    return json.loads(path.read_text())


def expected_stages(plan, directory):
    command = plan["command"]
    value = lambda key: command[command.index(key) + 1]
    root = directory / "run"
    if Path(value("--output")) != root or value("--cpu") != "32":
        raise ValueError("Wrong integrated output or CPU scope")
    assets = Path(value("--assets"))
    environment = dict(HOME=str(root / "home"), PATH="/usr/bin:/bin", LANG="C.UTF-8",
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    stages = []
    for name, wheels, lock in (("inference", assets / "wheels", assets / "benchmark_tools/results/publication_recovery_requirements_20260926.txt"),
                              ("reader", Path(value("--reader-wheels")), Path(value("--reader-lock")))):
        commands = install_commands(Path(value("--base-python")), Path(value("--installer-python")), wheels, lock, root / name)
        stages.extend((f"{name}_install_{i}", c, environment) for i, c in enumerate(commands))
    native = [str(root / "inference/bin/python"), "-I", "-B", str(assets / "benchmark_tools/run_publication_pipeline.py"),
              "--input", str(root / "input"), "--output", str(root / "native"), "--cpu", "32",
              "--aligner", str(assets / "mafft/bin/mafft"), "--tree-builder", str(assets / "FastTree")]
    stages.append(("native", native, dict(environment, MAFFT_BINARIES=str(assets / "mafft/libexec/mafft"))))
    stages.append(("readback_score", [str(root / "reader/bin/python"), "-I", "-B", command[4],
                   "--score-worker", value("--readers"), "--output", str(root)], environment))
    return stages


def verify_stages(root, complete, stages):
    if len(complete["outcomes"]) != 8 or len(stages) != 8:
        raise ValueError("Expected all eight completed stages")
    checked = []
    for (name, command, environment), outcome in zip(stages, complete["outcomes"]):
        start, finish = root / (name + "_started.json"), root / (name + "_finished.json")
        if (root / (name + "_failed.json")).exists():
            raise ValueError("Failed stage cannot be admitted")
        if read(start) != dict(command=command, environment=environment):
            raise ValueError("Stage command or environment differs")
        expected = dict(command=command, returncode=0, log=record(root / (name + ".log")))
        if outcome != expected or read(finish) != expected:
            raise ValueError("Stage completion differs")
        checked.extend([record(start), record(finish), expected["log"]])
    return checked


def verify_private_inputs(original, private, root):
    validate_data(original)
    validate_data(private)
    if {k: v for k, v in private.items() if k not in ("fasta", "references", "uncertain")} != {
            k: v for k, v in original.items() if k not in ("fasta", "references", "uncertain")}:
        raise ValueError("Private dataset metadata differs")
    for role in ("fasta", "references", "uncertain"):
        directory = root / ("input" if role == "fasta" else "scoring_inputs/" + role)
        expected = [dict(row, path=str(directory / Path(row["path"]).name)) for row in original[role]]
        if private[role] != expected or set(directory.iterdir()) != {Path(r["path"]) for r in expected}:
            raise ValueError("Private input inventory or byte pins differ")


def admission(directory, job):
    if job != JOB:
        raise ValueError("Admission is bound to the prespecified job")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    if scheduler["AllocCPUS"] != "32" or scheduler["NodeList"] != "bizon" or scheduler["ReqMem"] not in {"128G", "128Gn"}:
        raise ValueError("Scheduler resource allocation differs")
    plan_path = directory / "plan.json"
    plan = validate_plan(plan_path, PLAN_SHA)
    plan_record = record(plan_path)
    submitted, started, executed = [read(directory / name) for name in ("submission.json", "execution_started.json", "execution.json")]
    if (submitted["job_id"] != str(job) or submitted["plan"] != plan_record
            or started != dict(plan=plan_record, job_id=str(job))
            or executed["status"] != "integrated_complete_pending_independent_admission"
            or executed["returncode"] != 0 or executed["job_id"] != str(job) or executed["plan"] != plan_record
            or executed["command"] != ["/usr/bin/time", "-v", "-o", str(directory / "time.txt"), *plan["command"]]):
        raise ValueError("Scheduler submission and execution binding differs")
    check(submitted["protocol"])
    protocol = read(Path(submitted["protocol"]["path"]))
    if protocol["plan"] != plan_record:
        raise ValueError("Presubmission protocol plan differs")
    root = directory / "run"
    if (root / "failure.json").exists() or (root / "native/failure.json").exists():
        raise ValueError("Failed workflow cannot be admitted")
    complete = read(root / "complete.json")
    if (executed["complete"] != record(root / "complete.json")
            or complete["status"] != "integrated_install_inference_readback_scoring_complete"
            or complete["dataset"] != "orthobench" or complete["native_checkpoint_reuse"] is not False
            or complete["controlled_timing"] is not False or complete["started"] != record(root / "started.json")
            or complete["score"] != record(root / "score.json")):
        raise ValueError("Integrated completion binding differs")
    source = Path(install_commands.__code__.co_filename)
    if record(source)["sha256"] != record(Path(plan["command"][4]))["sha256"]:
        raise ValueError("Controller helper differs from the frozen controller")
    checked = verify_stages(root, complete, expected_stages(plan, directory))
    workflow_started = read(root / "started.json")
    pinned = {r["path"]: r for r in plan["checked_records"]}
    if (workflow_started["source"] != record(Path(plan["command"][4]))
            or workflow_started["cpu"] != 32 or workflow_started["dataset"] != "orthobench"
            or workflow_started["attempts"] != 1 or workflow_started["native_checkpoint_reuse"] is not False
            or not workflow_started["inputs"]
            or any(pinned.get(r["path"]) != r for r in workflow_started["inputs"])):
        raise ValueError("Integrated start scope or source/input binding differs")
    native_command = expected_stages(plan, directory)[6][1]
    native_started = read(root / "native/started.json")
    if (native_started["command"] != native_command[3:]
            or native_started["executable"] != native_command[0]
            or native_started["attempts"] != 1 or native_started["checkpoint_reuse"] is not False):
        raise ValueError("Native command, interpreter or attempt scope differs")
    original = read(directory / "data.json")
    private = read(root / "data.json")
    verify_private_inputs(original, private, root)
    if executed["logs"] != [record(directory / n) for n in ("workflow.log", "time.txt")]:
        raise ValueError("Missing launcher resource or execution logs")
    for row in executed["logs"]:
        check(row)
    checked.extend(record(directory / n) for n in ("submission.json", "execution_started.json", "execution.json"))
    checked.extend([record(root / "complete.json"), record(root / "score.json"), record(root / "data.json")])
    return plan, dict(accounting=accounting, scheduler=scheduler, checked_records=checked)


def audit(directory, job, output):
    if output.exists():
        raise FileExistsError(output)
    plan, admitted = admission(directory, job)
    root = directory / "run"
    output.mkdir(parents=True)
    save(output / "native_admission.json", admitted)
    package_audits = {}
    for name, wheels, lock in (("inference", directory / "assets/wheels", directory / "assets/benchmark_tools/results/publication_recovery_requirements_20260926.txt"),
                              ("reader", directory / "reader_wheels", directory / "reader_requirements.txt")):
        report = read(root / (name + "_install.json"))
        rows = local_install_wheels(report, wheels)
        verify_lock(lock.read_text(), rows)
        probe = "import json,importlib.metadata as m;print(json.dumps([dict(name=d.metadata['Name'],version=d.version) for d in m.distributions()]))"
        observed = json.loads(subprocess.check_output([str(root / name / "bin/python"), "-I", "-c", probe], text=True, timeout=60))
        inventory = verify_inventory(report, observed)
        payloads = [dict(name=w["name"], **installed_payload(Path(w["wheel"]["path"]), root / name / "lib/python3.10/site-packages")) for w in rows]
        package_audits[name] = dict(inventory=inventory, packages=payloads)
    save(output / "package_audits.json", package_audits)
    readers = directory / "readers"
    code = ("import sys;from pathlib import Path;sys.path.insert(0," + repr(str(readers)) + ");"
            "from benchmark_tools.audit_publication_pipeline import audit;"
            "audit(Path(" + repr(str(root / "native")) + "),Path(" + repr(str(output / "scientific")) + "))")
    command = [str(root / "reader/bin/python"), "-I", "-B", "-c", code]
    environment = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with (output / "scientific.log").open("x") as log:
        subprocess.run(command, env=environment, cwd=output, stdout=log, stderr=subprocess.STDOUT,
                       timeout=3600, check=True)
    independent = read(output / "scientific/result.json")
    recorded = read(root / "score.json")
    if independent["summary"] != recorded["summary"] or recorded["scientific_readback"] != record(root / "readback/result.json"):
        raise ValueError("Independent scientific readback differs")
    data = read(root / "data.json")
    universe = fasta_ids([Path(r["path"]) for r in data["fasta"]])
    if len(universe) != 251378:
        raise ValueError("Wrong full input universe")
    paths = dict(historical=Path(plan["historical_partition"]["path"]),
                 current=root / "native/inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv")
    partitions = {k: read_root_hogs(p, universe) for k, p in paths.items()}
    if recorded["dataset"] != "orthobench" or recorded["genes"] != len(universe) or recorded["groups"] != len(partitions["current"]):
        raise ValueError("Recorded score dataset dimensions differ")
    refs, uncertain = [{Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines()) for r in data[role]}
                       for role in ("references", "uncertain")]
    scorer = record(Path(score_partition.__code__.co_filename))
    if scorer["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    scores = {k: score_partition(v, refs, uncertain) for k, v in partitions.items()}
    if scores["current"] != recorded["score"] or scores["historical"] != read(Path(plan["historical_readback"]["path"]))["scores"]["current"]:
        raise ValueError("Recomputed score differs from its own recorded result")
    _, after = admission(directory, job)
    if after != admitted:
        raise ValueError("Admission changed during independent verification")
    result = dict(status="full_integrated_orthobench_independently_admitted", job_id=job,
        plan=record(directory / "plan.json"), source=record(__file__), scorer=scorer,
        scores=scores, comparison=comparison(partitions["historical"], partitions["current"], scores["historical"], scores["current"]),
        partitions={k: record(p) for k, p in paths.items()}, summary=independent["summary"],
        evidence=[record(output / n) for n in ("native_admission.json", "package_audits.json", "scientific/result.json")],
        historical_scores_replaced=False, controlled_timing=False, publication_ready=False,
        limitations=["Shared-host full workflow reproduction, not independent biological validation or controlled timing.",
            "Post-run wheel payload audit excludes generated metadata/bytecode and relocated non-site data.",
            "Cross-host restoration, data rights and remaining scientific requirements remain separate."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(audit(args.directory.resolve(), args.job, args.output.resolve())["status"])
