"""Convert corrected phylogenetic predictions after fresh independent admission."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_factorial_pairs import write_native_pairs
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs

ADMITTER_COMMIT = "5a9d18d8fa5cf26bd94c2e1e5f83089c3e7fa3e6"


def validate_native(native, index, candidate_record, prepared_record, source):
    if type(index) is not int or index not in (1, 3, 5, 7):
        raise ValueError("Only R-on cells use native phylogenetic pairs")
    if (native["status"] != "corrected_qfo_native_pair_output_verified"
            or native["cell"] != CELLS[index] or native["index"] != index // 2
            or native["accuracy_evaluated"] is not False or native["scoring_admitted"] is not False
            or native["publication_ready"] is not False):
        raise ValueError("Require unscored corrected native admission for this cell")
    if (native["candidate_admission"] != candidate_record or native["prepared"] != prepared_record
            or native["source"] != source):
        raise ValueError("Native candidate/preparation/admitter binding differs")
    if type(native["native_pair_count"]) is not int or native["native_pair_count"] <= 0:
        raise ValueError("Invalid admitted pair count")
    if native["native_pairs"] not in native["checked_records"]:
        raise ValueError("Native pair file is not in checked inventory")


def prepare(root, index, candidate_path, candidate_sha, candidate_job, native_path, native_sha):
    if type(index) is not int or index not in (1, 3, 5, 7):
        raise ValueError("Require R-on index 1, 3, 5 or 7")
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled 2-CPU native conversion")
    admission, manifest, cell, _, _, _ = verify_admission(root, candidate_path, candidate_sha, candidate_job, index // 2)
    native = read_frozen(native_path, native_sha)
    executor = root / "benchmarks/work/publication_qfo_corrected_reconcile_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER_COMMIT:
        raise ValueError("Independent admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    checker = executor / "benchmark_tools/admit_qfo_corrected_factorial_cell.py"
    candidate_record = record(candidate_path)
    validate_native(native, index, candidate_record, admission["prepared_manifest"], record(checker))
    if cell["label"] != CELLS[index] or cell["reconciliation"] is not True:
        raise ValueError("Wrong corrected reconciliation cell")
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    mapping = mappings[0]
    sources = [record(__file__), record(checker), *[record(Path(__file__).with_name(n)) for n in (
        "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py", "simulation_method_outputs.py",
        "admit_qfo_corrected_factorial_cell.py", "run_qfo_corrected_factorial_cell.py")],
        record(Path(__file__).resolve().parent.parent / "qfo_benchmark/og_to_pairwise.py")]
    checked = [record(native_path), candidate_record, admission["prepared_manifest"], record(environment_path),
               mapping, *native["checked_records"], *sources]
    for item in checked:
        check(item)
    output = root / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / cell["label"]
    output.mkdir(parents=True, exist_ok=False)
    recheck = output / "native_admission_recheck.json"
    command = [sys.executable, str(checker), "--root", str(root), "--index", str(index // 2),
        "--job", native["scheduler"]["JobIDRaw"], "--admission", str(candidate_path),
        "--admission-sha256", candidate_sha, "--admission-job", str(candidate_job), "--output", str(recheck)]
    report = {"status": "preparing", "index": index, "cell": cell["label"],
        "participant": "ohmm_qfo_corrected_factorial_" + cell["label"], "job_id": os.environ["SLURM_JOB_ID"],
        "candidate_admission": candidate_record, "native_admission": record(native_path),
        "prepared": admission["prepared_manifest"], "mapping": mapping, "input_fastas": manifest["input_fastas"],
        "native_input": native["native_pairs"], "native_recheck_command": command,
        "semantics": "native phylogenetically inferred pairs", "checked_records": checked,
        "accuracy_evaluated": False, "publication_ready": False}
    try:
        with (output / "native_recheck.log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True, cwd=root)
        if json.loads(recheck.read_text()) != native:
            raise ValueError("Fresh independent native admission differs from supplied report")
        checked.append(record(recheck))
        native_manifest_path = native["native_group_integrity"]["native_manifest"]
        check(native_manifest_path)
        owners, _ = gene_ownership(manifest, json.loads(Path(native_manifest_path["path"]).read_text()),
                                   Path(cell["candidate_partition"]))
        pairs, filtered = output / "pairs.partial.tsv", output / "pairs.qfo.partial.tsv"
        count = write_native_pairs(Path(native["native_pairs"]["path"]), pairs, owners, native["native_pair_count"])
        total, retained = filter_pairs(pairs, filtered, load_mapping(Path(mapping["path"])))
        if not 0 < count == total == retained:
            raise ValueError("Pair count mismatch or corrected reference mapping loss")
        for item in checked:
            check(item)
        verify_admission(root, candidate_path, candidate_sha, candidate_job, index // 2)
        final, mapped = output / "pairs.tsv", output / "pairs.qfo.tsv"
        pairs.rename(final)
        filtered.rename(mapped)
        report.update(status="corrected_factorial_native_pairs_prepared_unscored", pairs=record(final),
            filtered_pairs=record(mapped), total_pairs=total, retained_pairs=retained, removed_mapping_pairs=0,
            native_admission_recheck=record(recheck), recheck_log=record(output / "native_recheck.log"))
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
    for name in ("root", "candidate-admission", "native-admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("candidate-sha256", "candidate-job", "native-sha256"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--index", type=int, choices=(1, 3, 5, 7), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index, args.candidate_admission.resolve(), args.candidate_sha256,
            args.candidate_job, args.native_admission.resolve(), args.native_sha256)
