"""Assemble fully terminal simulation panels with explicit native failures."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.run_simulation_methods import read_frozen, verify_inputs
from benchmark_tools.simulation_conditions import score_pairs
from benchmark_tools.simulation_method_outputs import load_predictions
from benchmark_tools.summarize_simulation_panel import METHODS, PANELS, summarize
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure, validate_orthofinder, validate_orthohmm

TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL", "DEADLINE", "PREEMPTED"}


def terminal_tasks(accounting, array_id, count):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    tasks = []
    for index in range(count):
        matches = [r for r in rows if r.get("JobID") == f"{array_id}_{index}"]
        if len(matches) != 1 or matches[0]["State"].split()[0] not in TERMINAL:
            raise ValueError(f"Task {array_id}_{index} not uniquely terminal; do not score a partial panel")
        row = matches[0]
        if row["State"] == "COMPLETED" and row["ExitCode"] != "0:0":
            raise ValueError("Contradictory scheduler completion evidence")
        tasks.append(row)
    return tasks


def input_universe(records, truth):
    owners = {}
    species = []
    for record in records:
        path = Path(record["absolute_path"])
        verify_file(path, record)
        name = path.stem
        if name in species:
            raise ValueError("Duplicate input species")
        species.append(name)
        for sequence in SeqIO.parse(path, "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty input sequence")
            owners[sequence.id] = name
    if len(owners) != truth["extant_genes"] or set(species) != set(truth["species"]):
        raise ValueError("Truth and sequence input universes differ")
    if len(truth["ortholog_pairs"]) != truth["ortholog_pair_count"]:
        raise ValueError("Truth pair-count metadata differs")
    return owners, species


def verify_execution_status(status, dataset, task, verified, manifest_hash, generation_hash, executor):
    if status["dataset"] != dataset["label"] or status["verified_inputs"] != verified:
        raise ValueError("Execution dataset or verified inputs changed")
    provenance = status["provenance"]
    if provenance["method_manifest_sha256"] != manifest_hash or provenance["generation_manifest_sha256"] != generation_hash:
        raise ValueError("Execution manifest provenance mismatch")
    if str(provenance["slurm_job_id"]) != task["JobIDRaw"] or str(provenance["slurm_array_task_id"]) != task["JobID"].split("_")[-1]:
        raise ValueError("Execution belongs to a different scheduler task")
    expected = [file_record(executor / "benchmark_tools" / name, executor / "benchmark_tools") for name in
                ("run_simulation_methods.py", "verify_simulation_histories.py", "run_simulation_generation.py", "benchmark_production.py")]
    if provenance["sources"] != expected:
        raise ValueError("Inference executor sources differ from recorded versions")
    if task["State"] == "COMPLETED" and status["status"] not in {"finished_pending_native_validation", "inapplicable"}:
        raise ValueError("Scheduler success without terminal executor evidence")


def admit_method(method, dataset, status, verified, manifest):
    parent = "orthofinder_full" if method == "orthofinder_sequence_only" else method
    record = status["methods"].get(parent)
    config = dataset["methods"][parent]
    if record is None or record["status"] == "running":
        return {"status": "failed", "failure_stage": "execution_interrupted",
                "reason": "Terminal scheduler task has no completed process record"}
    if record["argv"] != config["argv"]:
        raise ValueError("Recorded method command differs from frozen command")
    if record["status"] == "failed":
        return {"status": "failed", "failure_stage": "execution",
                "reason": record.get("error") or f"Method process exit code {record.get('exit_code')}",
                "exit_code": record.get("exit_code")}
    if record["status"] != "process_succeeded":
        raise ValueError("Unknown process record status")
    try:
        if parent == "orthofinder_full":
            evidence = validate_orthofinder(config, record, verified["inputs"])
        else:
            evidence = validate_orthohmm(parent, config, record, verified["inputs"], manifest)
    except NativeOutputFailure as error:
        return {"status": "failed", "failure_stage": "native_output", "reason": str(error), "exit_code": record["exit_code"]}
    return {"status": "admitted", "native_validation": evidence}


def dataset_records(dataset, task, manifest, manifest_hash, generation, generation_hash, panel, executor):
    verified = verify_inputs(dataset, generation, panel, generation_hash)
    evidence = Path(dataset["methods"]["orthohmm_high_sensitivity"]["output"]).parent / "execution"
    status_path = evidence / "status.json"
    common = {"condition": dataset["condition"], "seed": dataset["seed"], "scheduler": task}
    if not status_path.exists():
        if task["State"] == "COMPLETED":
            raise ValueError("Successful scheduler task has no execution evidence")
        return [{**common, "method": method, "status": "failed", "failure_stage": "execution_preflight_or_scheduler",
                 "reason": "Terminal unsuccessful scheduler task has no execution record"} for method in METHODS]
    status = json.loads(status_path.read_text())
    verify_execution_status(status, dataset, task, verified, manifest_hash, generation_hash, executor)
    common["execution_evidence"] = dict(file_record(status_path, evidence), absolute_path=str(status_path))
    if verified["status"] == "inapplicable":
        return [{**common, "method": method, "status": "inapplicable", "reason": verified["reason"]} for method in METHODS]
    truth_path = Path(dataset["truth"])
    truth = json.loads(truth_path.read_text())
    owners, species = input_universe(verified["inputs"], truth)
    common["truth_sha256"] = verified["truth"]["sha256"]
    result = []
    for method in METHODS:
        row = {**common, "method": method, **admit_method(method, dataset, status, verified, manifest)}
        parent = "orthofinder_full" if method == "orthofinder_sequence_only" else method
        row["resource_evidence"] = [dict(file_record(path, evidence), absolute_path=str(path)) for path in
                                    (evidence / f"{parent}.time.log", evidence / f"{parent}.log") if path.is_file()]
        if method == "orthofinder_sequence_only":
            row.update(parent_method=parent, independent_timing=False)
        if row["status"] == "admitted":
            pairs, artifacts = load_predictions(method, Path(dataset["methods"][method]["output"]), owners, species)
            inventoried = {Path(r["absolute_path"]).resolve() for r in status["methods"][parent]["outputs"]}
            if not {p.resolve() for p in artifacts}.issubset(inventoried):
                raise ValueError("Converted predictions not present in verified output inventory")
            row.update(status="complete", score=score_pairs(pairs, truth["ortholog_pairs"], owners),
                       prediction_artifacts=[dict(file_record(p, p.parent), absolute_path=str(p)) for p in artifacts])
        result.append(row)
    return result


def verify_scoring_dependencies(manifest):
    # Preserve frozen conversion/scoring semantics; admission gates have their own provenance.
    names = {"simulation_method_outputs.py", "simulation_conditions.py", "orthofinder_to_pairwise.py",
             "orthofinder_mcl_to_orthogroups.py", "report_ygob_validation.py", "score_ygob_groups.py"}
    records = {r["path"]: r for r in manifest["adapter_sources"]}
    for name in names:
        verify_file(Path(__file__).with_name(name), records[name])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--panel-variant", choices=sorted(PANELS), required=True)
    parser.add_argument("--array-job", type=int, required=True)
    parser.add_argument("--executor", type=Path, required=True)
    parser.add_argument("--executor-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to replace assembled scientific results")
    manifest = read_frozen(args.manifest, args.manifest_sha256)
    generation_hash = manifest["generation_manifest"]["sha256"]
    generation = read_frozen(Path(manifest["generation_manifest"]["absolute_path"]), generation_hash)
    if generation.get("panel_variant", "fixed_length_v1") != args.panel_variant:
        raise ValueError("Requested analysis variant differs from frozen generation manifest")
    datasets = manifest["datasets"]
    if len(datasets) != 70 or {d["seed"] for d in datasets} != set(PANELS[args.panel_variant][0]):
        raise ValueError("Wrong scientific panel dimensions or seeds")
    executor = args.executor.resolve()
    commit = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    if commit != args.executor_commit:
        raise ValueError("Wrong frozen inference executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--quiet", "HEAD", "--", "benchmark_tools"], check=True)
    accounting = subprocess.check_output(["sacct", "-j", str(args.array_job), "--parsable2",
                                         "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = terminal_tasks(accounting, args.array_job, 70)
    verify_scoring_dependencies(manifest)
    rows = []
    for dataset, task in zip(datasets, tasks):
        rows.extend(dataset_records(dataset, task, manifest, args.manifest_sha256, generation, generation_hash,
                                    args.panel.resolve(), executor))
    report = summarize(rows, args.panel_variant)
    report.update(method_manifest_sha256=args.manifest_sha256, generation_manifest_sha256=generation_hash,
                  inference_executor_commit=commit, accounting_raw=accounting,
                  command=[sys.executable, *sys.argv], python=sys.version,
                  sources=[file_record(Path(__file__).with_name(n), Path(__file__).parent) for n in
                           ("assemble_simulation_results.py", "validate_simulation_outputs.py", "summarize_simulation_panel.py",
                            "run_simulation_methods.py", "verify_simulation_histories.py", "verify_ygob_validation.py",
                            "run_simulation_generation.py", "benchmark_production.py")])
    report["limitations"].extend(["Original fixed-length and heterogeneous-length panels are not pooled.",
        "Resource logs were obtained on a shared machine, not controlled scaling hardware.",
        "A sequence-only checkpoint is admitted only when its full parent passes native completion and finite-weight checks."])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"Assembled {len(rows)} explicit method outcomes for {args.panel_variant}")


if __name__ == "__main__":
    main()
