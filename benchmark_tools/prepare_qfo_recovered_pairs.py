"""Prepare all four admitted QfO stage predictions without running inference or scoring."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_SHA = "905c8accae2a0c98e73e028ec30e8a40d5a3d36d3223169afb6bedc01127758e"
STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined")


def record(path):
    path = Path(path).resolve()
    return {**file_record(path, path.parent), "path": str(path)}


def stage_outputs(admission):
    if (admission["status"] != "checked_full_replay_recovered_verified" or admission["accuracy_evaluated"] is not False
            or admission["run_version"] != "v2" or admission["recovery"]["inference_rerun"] is not False):
        raise ValueError("Require exact admitted recovered replay")
    stages = admission["native_replay"]["stages"]
    if [r["label"] for r in stages] != list(STAGES):
        raise ValueError("Require all four stages in frozen order")
    coverage = admission["coverage"]
    if [r["stage"] for r in coverage] != list(STAGES):
        raise ValueError("Incomplete admitted stage coverage")
    for stage, row in zip(stages, coverage):
        if stage["output"] != row["observed"] or row["partition_equal"] is not True:
            raise ValueError("Partition differs from admitted coverage")
    if len({r["output"]["path"] for r in stages}) != 4:
        raise ValueError("Stage paths are not unique")
    return stages


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    admission_path = root / "benchmark_tools/results/qfo_checked_full_replay_recovered_verified_20260917.json"
    admission = read_frozen(admission_path, ADMISSION_SHA)
    stages = stage_outputs(admission)
    accounting = subprocess.check_output(["sacct", "-j", "21480", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21480)
    check(admission["source"])
    check(admission["source_report"])
    for item in admission["provenance_checked"]:
        check(item)
    inputs = root / "qfo_benchmark/input"
    fastas = sorted(inputs.glob("*.fasta"))
    if len(fastas) != 78:
        raise ValueError("Unexpected QfO input species inventory")
    input_records = [record(p) for p in fastas]
    audit_link = admission["recovery"]["input_recheck"]
    check(audit_link)
    audited_inputs = json.loads(Path(audit_link["path"]).read_text())["input_fastas"]
    if input_records != sorted(audited_inputs, key=lambda r: r["path"]):
        raise ValueError("Conversion FASTAs differ from admitted replay inputs")
    converter = root / "qfo_benchmark/og_to_pairwise.py"
    mapping = root / "qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz"
    sources = [record(__file__), record(converter), record(Path(__file__).with_name("qfo_filter_pairs.py"))]
    mapping_record = record(mapping)
    valid_ids = load_mapping(mapping)
    output.mkdir(parents=True)
    report = {"status": "preparing", "accuracy_evaluated": False, "admission": record(admission_path),
              "admission_scheduler": scheduler, "sources": sources, "mapping": mapping_record,
              "input_fastas": input_records, "job_id": os.environ.get("SLURM_JOB_ID"), "stages": [],
              "semantics": "Cross-species clique pairs from complete partitions; not native reconciled ortholog pairs"}
    try:
        for index, stage in enumerate(stages):
            check(stage["output"])
            directory = output / stage["label"]
            directory.mkdir()
            pairs = directory / "pairs.tsv"
            filtered = directory / "pairs.qfo.tsv"
            command = [sys.executable, str(converter), stage["output"]["path"], str(inputs)]
            with pairs.open("x") as handle, (directory / "conversion.log").open("x") as log:
                subprocess.run(command, stdout=handle, stderr=log, check=True)
            total, retained = filter_pairs(pairs, filtered, valid_ids)
            if total <= 0 or retained <= 0 or retained > total:
                raise ValueError("Invalid converted/retained pair counts")
            report["stages"].append({"index": index, "stage": stage["label"],
                "participant": f"ohmm_checked_v2_{index}", "partition": stage["output"],
                "command": command, "pairs": record(pairs), "filtered_pairs": record(filtered),
                "conversion_log": record(directory / "conversion.log"), "total_pairs": total,
                "retained_pairs": retained, "removed_mapping_pairs": total - retained})
            check(stage["output"])
            (output / "progress.json").write_text(json.dumps({"stages_prepared": index + 1}) + "\n")
        if input_records != [record(p) for p in sorted(inputs.glob("*.fasta"))] or record(mapping) != mapping_record:
            raise ValueError("Input FASTAs or mapping changed during conversion")
        for source in sources:
            if record(source["path"]) != source:
                raise ValueError("Conversion source changed")
        read_frozen(admission_path, ADMISSION_SHA)
        report["status"] = "four_stage_pairs_prepared_unscored"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
