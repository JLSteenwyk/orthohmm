"""Convert admitted sequence-control groups to QfO cross-species clique pairs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.compare_qfo_search_coverage import frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_corrected_group_pairs import expected_pairs
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_qfo_sequence_search_control import verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from qfo_benchmark.og_to_pairwise import index_genes

ADMITTER = "f1e21b09c28f270dc3ef2243bdcad86f212b58a0"


def validate_admission(admission, label):
    if (label not in ("all_hits", "top100") or admission["variant"] != label
            or admission["status"] != "corrected_sequence_graph_admitted"
            or admission["accuracy_evaluated"] is not False or admission["publication_ready"] is not False
            or admission["numeric"]["status"] != "numeric_checkpoint_verified"
            or admission["numeric"]["summary"]["genes"] != 984137
            or admission["numeric"]["summary"]["species"] != 78
            or [stage["label"] for stage in admission["coverage"]] != ["multipass", "multipass_refined"]
            or admission["prediction"] != admission["coverage"][-1]["output"]
            or admission["prediction"] not in admission["checked_records"]):
        raise ValueError("Wrong sequence graph admission or refined prediction")
    return admission["prediction"]


def convert(partition, fasta, converter, valid_ids, directory):
    expected = expected_pairs(partition, index_genes(fasta))
    command = [sys.executable, str(converter), str(partition), str(fasta)]
    pairs, filtered = directory / "pairs.partial.tsv", directory / "pairs.qfo.partial.tsv"
    with pairs.open("x") as stream, (directory / "conversion.log").open("x") as log:
        subprocess.run(command, stdout=stream, stderr=log, check=True)
    total, retained = filter_pairs(pairs, filtered, valid_ids)
    if total != expected or retained != total:
        raise ValueError("Independent pair count differs or corrected-input mapping loss")
    return dict(command=command, expected_pairs=expected, total_pairs=total,
                retained_pairs=retained, removed_mapping_pairs=0)


def prepare(root, label, admission_path, admission_sha, admission_job):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled 2-CPU conversion")
    scheduler = completed(admission_job, 2, "192G")
    executor = frozen(root, "publication_qfo_sequence_graph_admission_v1", ADMITTER)
    admission = read_frozen(admission_path, admission_sha)
    prediction = validate_admission(admission, label)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_sequence_graph.py"):
        raise ValueError("Unexpected graph admission source")
    search_path = root / "benchmarks/work/qfo_sequence_search_control_v1/manifest.json"
    search = verify_plan(search_path)
    if record(search_path) not in admission["checked_records"]:
        raise ValueError("Search plan absent from graph provenance")
    fastas = search["inputs"]
    parents = {Path(item["path"]).parent for item in fastas}
    if len(parents) != 1:
        raise ValueError("Ambiguous input FASTA directory")
    fasta = next(iter(parents))
    if {str(path) for path in fasta.glob("*.fasta")} != {item["path"] for item in fastas}:
        raise ValueError("Input FASTA inventory differs")
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [item for item in environment["reference_files"] if Path(item["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require unique frozen reference mapping")
    mapping = mappings[0]
    converter = Path(__file__).resolve().parent.parent / "qfo_benchmark/og_to_pairwise.py"
    checked = [record(admission_path), admission["source"], *admission["checked_records"],
               record(search_path), *fastas, record(environment_path), mapping, record(__file__), record(converter),
               record(Path(__file__).with_name("qfo_filter_pairs.py")),
               record(Path(__file__).with_name("prepare_qfo_corrected_group_pairs.py"))]
    for item in checked:
        check(item)
    directory = root / "benchmarks/results/qfo_sequence_pairs_v1" / label
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(status="preparing", source=record(__file__), variant=label,
        participant="ohmm_qfo_corrected_sequence_" + label, job_id=os.environ["SLURM_JOB_ID"],
        graph_admission=record(admission_path), admission_scheduler=scheduler, prediction=prediction,
        mapping=mapping, input_fastas=fastas, checked_records=checked,
        semantics="cross-species group-derived clique pairs", accuracy_evaluated=False, publication_ready=False)
    try:
        report.update(convert(Path(prediction["path"]), fasta, converter, load_mapping(Path(mapping["path"])), directory))
        for item in checked:
            check(item)
        for name in ("pairs", "pairs.qfo"):
            (directory / (name + ".partial.tsv")).rename(directory / (name + ".tsv"))
        report.update(status="corrected_sequence_group_pairs_prepared_unscored", pairs=record(directory / "pairs.tsv"),
                      filtered_pairs=record(directory / "pairs.qfo.tsv"), conversion_log=record(directory / "conversion.log"))
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (directory / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("admission-sha256", "admission-job"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--variant", required=True, choices=("all_hits", "top100"))
    args = parser.parse_args()
    prepare(args.root.resolve(), args.variant, args.admission.resolve(), args.admission_sha256, args.admission_job)
