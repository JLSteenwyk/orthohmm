"""Convert independently admitted CPM native pairs to the frozen QfO reference."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_parameter_pairs import convert
from benchmark_tools.run_qfo_cpm_phylogeny import verify_sources
from benchmark_tools.run_qfo_cpm_variant import ARMS
from benchmark_tools.run_simulation_methods import read_frozen

ADMISSION_JOB = "22070"
ADMITTER_COMMIT = "448c87d9330af1822598693c26c9b7dc546d3355"
ADMITTER_SHA = "395cb6ea4dd0c07ebfb55a4c526dbb9d9729c00b80ed24abdae8207daee881fa"


def completed_admission(accounting, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM pair conversion index")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{ADMISSION_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require successful terminal CPM native admission")
    return rows[0]


def validate_native(native, index, source, verified):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM pair conversion index")
    if (native["status"] != "cpm_native_pairs_verified_unscored"
            or native["arm"] != ARMS[index] or type(native["index"]) is not int or native["index"] != index
            or native["cell"]["label"] != "candidate_" + ARMS[index]
            or native["source"] != source or native["context"] != verified["context"]
            or native["candidate_admission"] != verified["admission_record"]
            or native["accuracy_evaluated"] is not False or native["scoring_admitted"] is not False
            or native["publication_ready"] is not False or type(native["native_pair_count"]) is not int
            or native["native_pair_count"] <= 0 or native["native_pairs"] not in native["checked_records"]):
        raise ValueError("Wrong or incomplete CPM native admission")


def prepare(root, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM pair conversion index")
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled two-CPU bizon CPM conversion array")
    accounting = subprocess.check_output(["sacct", "-j", ADMISSION_JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_admission(accounting, index)
    executor = root / "benchmarks/work/publication_qfo_cpm_native_admission_v2"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER_COMMIT:
        raise ValueError("CPM native admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    checker = executor / "benchmark_tools/admit_qfo_cpm_phylogeny.py"
    checker_record = record(checker)
    if checker_record["sha256"] != ADMITTER_SHA:
        raise ValueError("CPM native admission source changed")
    native_path = root / f"benchmarks/work/qfo_cpm_native_admission_{ADMISSION_JOB}_{index}.json"
    admission = record(native_path)
    native = read_frozen(native_path, admission["sha256"])
    verified = verify_sources(root, index)
    validate_native(native, index, checker_record, verified)
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    mapping = mappings[0]
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith(("benchmark_tools.", "qfo_benchmark.")) and getattr(m, "__file__", None)]
    checked = [record(__file__), admission, checker_record, *native["checked_records"], *native["helpers"],
               *verified["checked_records"], record(env_path), mapping, *helpers]
    for item in checked:
        check(item)
    output = root / "benchmarks/results/qfo_cpm_pairs_v1" / ARMS[index]
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=False)
    recheck = output / "native_admission_recheck.json"
    command = [sys.executable, "-B", str(checker), "--root", str(root), "--index", str(index), "--output", str(recheck)]
    report = {"status": "preparing", "arm": ARMS[index], "index": index, "context": verified["context"],
        "participant": "ohmm_qfo_parameter_" + ARMS[index], "job_id": os.environ["SLURM_JOB_ID"],
        "array_task_id": str(index), "admission_scheduler": scheduler, "admission_accounting": accounting,
        "native_admission": admission, "native_input": native["native_pairs"], "mapping": mapping,
        "input_fastas": verified["manifest"]["input_fastas"], "native_recheck_command": command,
        "semantics": "native phylogenetically inferred pairs", "source": record(__file__), "helpers": helpers,
        "checked_records": checked, "accuracy_evaluated": False, "publication_ready": False}
    try:
        with (output / "native_recheck.log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True, cwd=root)
        fresh = record(recheck)
        if read_frozen(recheck, fresh["sha256"]) != native:
            raise ValueError("Fresh CPM native admission differs from supplied report")
        checked.extend([fresh, record(output / "native_recheck.log")])
        metadata = native["native_group_integrity"]["native_manifest"]
        check(metadata)
        checked.append(metadata)
        owners, _ = gene_ownership(verified["manifest"], read_frozen(Path(metadata["path"]), metadata["sha256"]),
                                  Path(verified["arm"]["partition"]["path"]))
        counts = convert(native, owners, Path(mapping["path"]), output)
        for item in checked:
            check(item)
        if verify_sources(root, index) != verified:
            raise ValueError("CPM conversion inputs changed during execution")
        pairs, filtered = output / "pairs.tsv", output / "pairs.qfo.tsv"
        (output / "pairs.partial.tsv").rename(pairs)
        (output / "pairs.qfo.partial.tsv").rename(filtered)
        report.update(status="cpm_native_pairs_prepared_unscored", **counts,
            pairs=record(pairs), filtered_pairs=record(filtered), native_admission_recheck=fresh,
            conversion_counts=record(output / "conversion_counts.json"))
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index)
