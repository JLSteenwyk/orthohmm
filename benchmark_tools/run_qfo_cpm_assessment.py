"""Run frozen QfO endpoints for independently admitted CPM native pairs."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_qfo_cpm_phylogeny import verify_sources
from benchmark_tools.run_qfo_cpm_variant import ARMS
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen

CONVERSION_JOB = "21978"
CONVERTER_COMMIT = "ba4005be49a050c538306f2d7e64a0c4ac387eeb"
CONVERTER_SHA = "53124b59b10ee111151d2b4332f643d69306c9c2915736dbcca4a848b6e3a197"


def completed_conversion(accounting, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM assessment index")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{CONVERSION_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require successful terminal CPM pair conversion")
    return rows[0]


def validate_stage(stage, index, scheduler):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM assessment index")
    if (stage["status"] != "cpm_native_pairs_prepared_unscored"
            or type(stage["index"]) is not int or stage["index"] != index or stage["arm"] != ARMS[index]
            or stage["participant"] != "ohmm_qfo_parameter_" + ARMS[index]
            or stage["semantics"] != "native phylogenetically inferred pairs"
            or stage["job_id"] != scheduler["JobIDRaw"] or stage["array_task_id"] != str(index)
            or stage["accuracy_evaluated"] is not False or stage["publication_ready"] is not False):
        raise ValueError("Wrong CPM conversion identity or semantics")
    counts = [stage[k] for k in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")]
    if any(type(v) is not int for v in counts) or not 0 < counts[0] == counts[1] == counts[2] or counts[3] != 0:
        raise ValueError("Invalid CPM conversion or mapping counts")
    if any(stage["pairs"][k] != stage["filtered_pairs"][k] for k in ("bytes", "sha256")):
        raise ValueError("Reference filtering changed corrected predictions")
    if stage["native_admission_recheck"] not in stage["checked_records"]:
        raise ValueError("Fresh native admission missing")


def prepare(root, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM assessment index")
    accounting = subprocess.check_output(["sacct", "-j", CONVERSION_JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_conversion(accounting, index)
    path = root / "benchmarks/results/qfo_cpm_pairs_v1" / ARMS[index] / "results.json"
    pair_record = record(path)
    stage = read_frozen(path, pair_record["sha256"])
    validate_stage(stage, index, scheduler)
    converter = root / "benchmarks/work/publication_qfo_cpm_pairs_v1"
    if subprocess.check_output(["git", "-C", str(converter), "rev-parse", "HEAD"], text=True).strip() != CONVERTER_COMMIT:
        raise ValueError("Frozen CPM converter revision changed")
    subprocess.run(["git", "-C", str(converter), "diff", "--exit-code", "HEAD", "--",
                    "benchmark_tools", "orthohmm", "qfo_benchmark/og_to_pairwise.py"], check=True)
    source = record(converter / "benchmark_tools/prepare_qfo_cpm_pairs.py")
    if source["sha256"] != CONVERTER_SHA or source != stage["source"] or source not in stage["checked_records"]:
        raise ValueError("Wrong CPM converter source")
    verified = verify_sources(root, index)
    if stage["context"] != verified["context"] or stage["input_fastas"] != verified["manifest"]["input_fastas"]:
        raise ValueError("CPM conversion context differs from frozen arm")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if mappings != [stage["mapping"]]:
        raise ValueError("Conversion and assessment mapping differ")
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith("benchmark_tools.") and getattr(m, "__file__", None)]
    records = [pair_record, record(env_path), *environment_records(environment), *stage["checked_records"],
               *verified["checked_records"], stage["pairs"], stage["filtered_pairs"], stage["conversion_counts"], *helpers]
    for item in records:
        check(item)
    counts = json.loads(Path(stage["conversion_counts"]["path"]).read_text())
    if counts != {k: stage[k] for k in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")}:
        raise ValueError("Conversion count sidecar differs")
    output = root / "benchmarks/results/qfo_cpm_assessment_v1" / ARMS[index]
    work = root / "qfo_benchmark/w" / f"qcpv{index}"
    results = root / "qfo_benchmark/scoring" / f"cpm_v1_{index}"
    for namespace in (output, work, results):
        if namespace.exists() or namespace.is_symlink():
            raise FileExistsError(namespace)
    return {"status": "prepared_unrun", "index": index, "arm": ARMS[index], "context": verified["context"],
        "stage": stage, "source": record(__file__), "pairs_manifest": pair_record,
        "environment_manifest": record(env_path), "conversion_scheduler": scheduler,
        "conversion_accounting": accounting, "converter_commit": CONVERTER_COMMIT,
        "command": command_for(root, stage, environment, work, results), "cwd": str(output),
        "work": str(work), "results": str(results), "verified_records": records,
        "environment_overrides": environment["environment_overrides"], "accuracy_admitted": False}


def run(root, index):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "8"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled eight-CPU bizon CPM assessment array")
    report = prepare(root, index)
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"], array_task_id=str(index))
    with (output / "preflight.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    env = {**os.environ, **report["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            process = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        for item in [report["source"], *report["verified_records"]]:
            check(item)
        report["outputs"] = [record(p) for p in sorted(Path(report["results"]).rglob("*")) if p.is_file()]
        if process.returncode:
            raise RuntimeError(f"Native scoring failed: {process.returncode}")
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        log = output / "scoring.log"
        if log.exists():
            report["log"] = record(log)
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.index)
