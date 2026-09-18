"""Assess admitted QfO factorial pairs, reusing only identical baseline assessments."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.admit_qfo_recovered_assessment import ENV_SHA, PAIRS_SHA, EXECUTOR, validate_stage
from benchmark_tools.assemble_simulation_results import terminal_tasks

CONVERTER_COMMIT = "d786352dc57f3560bf245c86453044de32630190"
ADMITTED_SHA = "89b683d0fc7fe9964ce5b6182832bb8eb28758fad3f6bbb4d706926143805956"


def cell_label(index):
    if type(index) is not int or not 0 <= index < 8:
        raise ValueError("Require factorial cell index0-7")
    return f"p{index // 4}_c{(index // 2) % 2}_r{index % 2}"


def verify_conversion_identity(report, index, accounting):
    label = cell_label(index)
    array = 21675 if index % 2 else 21674
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    matches = [r for r in rows if r["JobID"] == f"{array}_{index}"]
    if len(matches) != 1 or matches[0]["State"] != "COMPLETED" or matches[0]["ExitCode"] != "0:0":
        raise ValueError("Conversion task is not terminal successful")
    scheduler = matches[0]
    semantics = "native phylogenetically inferred pairs" if index % 2 else "cross-species group-derived clique pairs"
    expected = {"status": "cell_pairs_prepared_unscored", "accuracy_evaluated": False, "index": index,
                "cell": label, "participant": f"ohmm_qfo_factorial_{label}", "semantics": semantics,
                "array_job_id": str(array), "array_task_id": str(index), "job_id": scheduler["JobIDRaw"]}
    if any(report.get(k) != v for k, v in expected.items()):
        raise ValueError("Conversion identity or prediction semantics changed")
    if not 0 < report["retained_pairs"] <= report["total_pairs"] or report["removed_mapping_pairs"] != report["total_pairs"] - report["retained_pairs"]:
        raise ValueError("Conversion counts inconsistent")
    return scheduler


def verify_reuse_binding(conversion, native_stage, index):
    if index not in (0, 4):
        raise ValueError("Only identical unexpanded baselines may reuse assessments")
    if native_stage["status"] != "admitted" or native_stage["index"] != (1 if index == 0 else 3):
        raise ValueError("Wrong admitted baseline stage")
    old = native_stage["conversion"]
    if conversion.get("reused_conversion", {}).get("stage") != old["stage"]:
        raise ValueError("Baseline conversion was not reused")
    for key in ("pairs", "filtered_pairs"):
        if conversion[key] != old[key]:
            raise ValueError("Baseline prediction file differs")
    for key in ("total_pairs", "retained_pairs", "removed_mapping_pairs"):
        if conversion[key] != old[key]:
            raise ValueError("Baseline prediction coverage differs")
    a, b = conversion["candidate_partition"], old["partition"]
    if (a["sha256"], a["bytes"]) != (b["sha256"], b["bytes"]):
        raise ValueError("Baseline partition differs")


def run(root, index, pairs_sha):
    if Path.cwd().resolve() != root:
        raise ValueError("Run from repository verification directory")
    label = cell_label(index)
    pairs_path = root / "benchmarks/results/qfo_factorial_pairs_v1" / label / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    array = 21675 if index % 2 else 21674
    accounting = subprocess.check_output(["sacct", "-j", str(array), "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = verify_conversion_identity(stage, index, accounting)
    converter = root / "benchmarks/work/publication_qfo_factorial_pairs_v1"
    if subprocess.check_output(["git", "-C", str(converter), "rev-parse", "HEAD"], text=True).strip() != CONVERTER_COMMIT:
        raise ValueError("Converter checkout changed")
    subprocess.run(["git", "-C", str(converter), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "qfo_benchmark/og_to_pairwise.py"], check=True)
    if stage["source"] != record(converter / "benchmark_tools/prepare_qfo_factorial_pairs.py"):
        raise ValueError("Wrong pair converter source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    records = [*environment_records(environment), *stage["checked_records"], stage["pairs"], stage["filtered_pairs"], record(pairs_path)]
    for item in records:
        check(item)
    output = root / "benchmarks/results/qfo_factorial_assessment_v1" / label
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    report = {"status": "running", "accuracy_admitted": False, "source": record(__file__), "cell": label,
              "index": index, "stage": stage, "environment_manifest": record(env_path), "pairs_manifest": record(pairs_path),
              "conversion_scheduler": scheduler, "job_id": os.environ.get("SLURM_JOB_ID"), "verified_records": records,
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                ("admit_qfo_recovered_assessment.py", "run_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")]}
    (output / "preflight.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    try:
        if index in (0, 4):
            admission_path = root / "benchmark_tools/results/qfo_recovered_assessment_admitted_20260917.json"
            admission = read_frozen(admission_path, ADMITTED_SHA)
            old_pairs_path = root / "benchmarks/results/qfo_recovered_stage_pairs_v1/results.json"
            old_pairs = read_frozen(old_pairs_path, PAIRS_SHA)
            if admission["environment_manifest"] != record(env_path) or old_pairs["mapping"] != stage["mapping"] or old_pairs["input_fastas"] != stage["input_fastas"]:
                raise ValueError("Baseline reference environment or input provenance differs")
            for item in [admission["source"], *admission["helper_sources"]]:
                check(item)
            old_executor = root / "benchmarks/work/publication_qfo_recovered_assessment_v1"
            if subprocess.check_output(["git", "-C", str(old_executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
                raise ValueError("Original assessment executor changed")
            subprocess.run(["git", "-C", str(old_executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
            old_accounting = subprocess.check_output(["sacct", "-j", "21548", "--parsable2",
                "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
            original_index = 1 if index == 0 else 3
            tasks = terminal_tasks(old_accounting, 21548, 4)
            native = validate_stage(root, original_index, tasks[original_index], environment, old_pairs, old_executor, env_path, old_pairs_path)
            if native != admission["records"][original_index]:
                raise ValueError("Original assessment differs on independent revalidation")
            verify_reuse_binding(stage, native, index)
            report.update(status="admitted_reused_assessment", accuracy_admitted=True,
                reused_admission=record(admission_path), original_participant=native["participant"], admitted_stage=native,
                original_accounting_recheck=old_accounting,
                scoring_rerun=False, limitation="Reused saved native FAS sample; not an independent scoring replicate or new biological validation.")
        else:
            work = root / "qfo_benchmark/w" / f"qfx_{index}"
            results = root / "qfo_benchmark/scoring" / f"factorial_v1_{index}"
            if work.exists() or results.exists():
                raise FileExistsError("Existing scoring namespace")
            command = command_for(root, stage, environment, work, results)
            env = {**os.environ, **environment["environment_overrides"],
                   "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
            report.update(command=command, cwd=str(output), work=str(work), results=str(results))
            with (output / "scoring.log").open("x") as log:
                process = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
            report.update(exit_code=process.returncode, log=record(output / "scoring.log"),
                outputs=[record(p) for p in sorted(results.rglob("*")) if p.is_file()],
                status="process_succeeded_pending_independent_admission" if process.returncode == 0 else "scoring_failed")
        for item in [*records, report["source"], *report["helper_sources"]]:
            check(item)
        read_frozen(env_path, ENV_SHA)
        read_frozen(pairs_path, pairs_sha)
    except BaseException as error:
        report.update(status="failed", accuracy_admitted=False, error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(8), required=True)
    parser.add_argument("--pairs-sha256", required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.index, args.pairs_sha256)
    raise SystemExit(0 if result["status"] in {"admitted_reused_assessment", "process_succeeded_pending_independent_admission"} else 1)
