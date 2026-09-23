"""Convert independently admitted parameter-variant native pairs for QfO."""

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
from benchmark_tools.prepare_qfo_factorial_pairs import write_native_pairs
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_qfo_parameter_phylogeny import VARIANTS, verify_sources
from benchmark_tools.run_simulation_methods import read_frozen

ADMITTER_COMMIT = "7f8ff221bb9dd0b2a3213496ba109bd776dfe13f"
ADMISSION_JOB = "22035"


def completed_admission(accounting, index):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown parameter variant")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{ADMISSION_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require successful terminal native admission")
    return rows[0]


def validate_native(native, index, source):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown parameter variant")
    if (native["status"] != "corrected_qfo_parameter_native_pairs_verified"
            or native["variant"] != VARIANTS[index] or native["index"] != index
            or native["cell"]["label"] != "candidate_" + VARIANTS[index]
            or native["source"] != source or native["accuracy_evaluated"] is not False
            or native["scoring_admitted"] is not False or native["publication_ready"] is not False
            or type(native["native_pair_count"]) is not int or native["native_pair_count"] <= 0
            or native["native_pairs"] not in native["checked_records"]):
        raise ValueError("Wrong or incomplete native parameter admission")


def convert(native, owners, mapping, output):
    pairs, filtered = output / "pairs.partial.tsv", output / "pairs.qfo.partial.tsv"
    count = write_native_pairs(Path(native["native_pairs"]["path"]), pairs, owners, native["native_pair_count"])
    total, retained = filter_pairs(pairs, filtered, load_mapping(mapping))
    counts = {"written_pairs": count, "total_pairs": total, "retained_pairs": retained,
              "removed_mapping_pairs": total - retained}
    with (output / "conversion_counts.json").open("x") as stream:
        json.dump(counts, stream, indent=2, sort_keys=True)
        stream.write("\n")
    if not 0 < count == total == retained:
        raise ValueError("Pair count mismatch or corrected reference mapping loss")
    return counts


def prepare(root, index):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown parameter variant")
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled two-CPU bizon conversion array")
    accounting = subprocess.check_output(["sacct", "-j", ADMISSION_JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_admission(accounting, index)
    executor = root / "benchmarks/work/publication_qfo_parameter_native_admission_v2"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER_COMMIT:
        raise ValueError("Native admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    checker = executor / "benchmark_tools/admit_qfo_parameter_phylogeny.py"
    native_path = root / f"benchmarks/work/qfo_parameter_native_admission_{ADMISSION_JOB}_{index}.json"
    admission_record = record(native_path)
    native = read_frozen(native_path, admission_record["sha256"])
    validate_native(native, index, record(checker))
    arm, manifest, _, _, _, _, inputs, _ = verify_sources(root, index)
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    mapping = mappings[0]
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith(("benchmark_tools.", "qfo_benchmark.")) and getattr(m, "__file__", None)]
    checked = [record(__file__), admission_record, record(checker), *native["checked_records"],
               *native["helpers"], *inputs, record(env_path), mapping, *helpers]
    for item in checked:
        check(item)
    output = root / "benchmarks/results/qfo_parameter_pairs_v1" / VARIANTS[index]
    output.mkdir(parents=True, exist_ok=False)
    recheck = output / "native_admission_recheck.json"
    command = [sys.executable, str(checker), "--root", str(root), "--index", str(index), "--output", str(recheck)]
    report = {"status": "preparing", "variant": VARIANTS[index], "index": index,
        "participant": "ohmm_qfo_parameter_" + VARIANTS[index], "job_id": os.environ["SLURM_JOB_ID"],
        "array_task_id": str(index), "admission_scheduler": scheduler, "admission_accounting": accounting,
        "native_admission": admission_record, "native_input": native["native_pairs"], "mapping": mapping,
        "input_fastas": manifest["input_fastas"], "native_recheck_command": command,
        "semantics": "native phylogenetically inferred pairs", "checked_records": checked,
        "accuracy_evaluated": False, "publication_ready": False}
    try:
        with (output / "native_recheck.log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True, cwd=root)
        if json.loads(recheck.read_text()) != native:
            raise ValueError("Fresh native admission differs from supplied report")
        checked.append(record(recheck))
        native_record = native["native_group_integrity"]["native_manifest"]
        check(native_record)
        owners, _ = gene_ownership(manifest, json.loads(Path(native_record["path"]).read_text()), Path(arm["partition"]["path"]))
        counts = convert(native, owners, Path(mapping["path"]), output)
        for item in checked:
            check(item)
        verify_sources(root, index)
        pairs, filtered = output / "pairs.tsv", output / "pairs.qfo.tsv"
        (output / "pairs.partial.tsv").rename(pairs)
        (output / "pairs.qfo.partial.tsv").rename(filtered)
        report.update(status="corrected_parameter_native_pairs_prepared_unscored", **counts,
            pairs=record(pairs), filtered_pairs=record(filtered), native_admission_recheck=record(recheck),
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
    parser.add_argument("--index", type=int, choices=range(4), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index)
