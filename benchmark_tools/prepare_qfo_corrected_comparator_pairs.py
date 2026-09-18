"""Convert admitted corrected Proteinortho/SonicParanoid native predictions."""

import argparse
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.proteinortho_to_pairwise import write_pairs as protein_pairs
from benchmark_tools.sonicparanoid_to_pairwise import write_pairs as sonic_pairs
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs

ENV_SHA = "e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc"
METHODS = {
    "proteinortho": ("corrected_proteinortho_native_graph_admitted_for_conversion", "graph_inventory", "validated_pair_rows"),
    "sonic": ("corrected_sonic_native_tables_admitted_for_conversion", "table_inventory", "distinct_pairs"),
}


def validate_admission(method, admission):
    status, key, count_key = METHODS[method]
    if admission["status"] != status or admission["accuracy_evaluated"] is not False:
        raise ValueError("Require corrected native admission for selected method")
    inventory = admission[key]
    if inventory["input_species"] != 78 or inventory["input_accessions"] != 984137:
        raise ValueError("Not the corrected complete input inventory")
    count = inventory[count_key]
    if type(count) is not int or count <= 0:
        raise ValueError("Invalid admitted native pair count")
    if method == "sonic":
        if inventory["species_pair_tables"] != 3003 or inventory["raw_relations"] != count + inventory["duplicate_relations"]:
            raise ValueError("Invalid native table/relation accounting")
    elif inventory["species_sections"] != 3003:
        raise ValueError("Incomplete native graph")
    return count


def convert(method, admission, output):
    expected = validate_admission(method, admission)
    with output.open("x") as stream:
        if method == "proteinortho":
            actual = protein_pairs(Path(admission["native_outputs"]["pairs"]["path"]), stream)
            duplicates = 0
        else:
            tables, actual, duplicates = sonic_pairs(Path(admission["pair_directory"]), stream)
            if tables != 3003 or duplicates != admission["table_inventory"]["duplicate_relations"]:
                raise ValueError("Native table/duplicate counts changed")
    if actual != expected:
        raise ValueError("Converted pair count differs from admission")
    return actual, duplicates


def prepare(root, method, admission_path, admission_sha):
    admission = read_frozen(admission_path, admission_sha)
    validate_admission(method, admission)
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    mapping = mappings[0]
    sources = [record(Path(__file__).with_name(name)) for name in (
        "proteinortho_to_pairwise.py", "sonicparanoid_to_pairwise.py", "qfo_filter_pairs.py")]
    checked = [record(__file__), record(admission_path), admission["source"], *admission["checked_records"],
               record(environment_path), mapping, *sources]
    for item in checked:
        check(item)
    directory = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1" / method
    directory.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "source": record(__file__), "method": method,
              "participant": "qfo_corrected_" + method, "admission": record(admission_path),
              "mapping": mapping, "checked_records": checked, "accuracy_evaluated": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "started_epoch": time.time(),
              "semantics": "native post-clustering graph" if method == "proteinortho" else "native species-pair relations"}
    try:
        partial = directory / "pairs.partial.tsv"
        total, duplicates = convert(method, admission, partial)
        filtered_partial = directory / "pairs.qfo.partial.tsv"
        observed, retained = filter_pairs(partial, filtered_partial, load_mapping(Path(mapping["path"])))
        if observed != total or retained != total:
            raise ValueError("Unexpected corrected-release mapping loss or count mismatch")
        for item in checked:
            check(item)
        pairs, filtered = directory / "pairs.tsv", directory / "pairs.qfo.tsv"
        partial.rename(pairs)
        filtered_partial.rename(filtered)
        report.update(status="corrected_comparator_pairs_prepared_unscored", pairs=record(pairs),
                      filtered_pairs=record(filtered), total_pairs=total, retained_pairs=retained,
                      removed_mapping_pairs=0, native_duplicate_relations=duplicates)
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        with (directory / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.method, args.admission.resolve(), args.admission_sha256)
