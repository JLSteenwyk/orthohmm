"""Convert admitted corrected FastOMA native pairs to frozen QfO inputs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.fastoma_distinct_pairs import write_pairs
from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMITTER = "857bb5e0d6adad9960ea58d98a3037097ee40d5f"
SEMANTICS = "native phylogenetically inferred pairs; supplied corrected OrthoFinder species tree"


def validate_admission(admission, root):
    if (admission["status"] != "corrected_fastoma_native_evidence_admitted"
            or admission["accuracy_evaluated"] is not False or admission["publication_ready"] is not False):
        raise ValueError("Require corrected FastOMA native admission")
    content, scheduler = admission["content"], admission["scheduler"]
    if (content["input_proteins"] != 984137 or content["species"] != 78
            or type(content["native_pair_rows"]) is not int or content["native_pair_rows"] <= 0):
        raise ValueError("Wrong corrected input scope or pair count")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["AllocCPUS"] != "180" or scheduler["ReqMem"] != "720G" or scheduler["NodeList"] != "bizon"):
        raise ValueError("Wrong native inference accounting")
    pairs = admission["native_pairs"]
    if pairs["path"] != str(root / "benchmarks/results/qfo_corrected_fastoma_v1/output/orthologs.tsv.gz"):
        raise ValueError("Unexpected corrected native pair path")
    inputs = admission["input_fastas"]
    if (len(inputs) != 78 or len({r["path"] for r in inputs}) != 78
            or any(Path(r["path"]).parent != root / "benchmarks/work/qfo_corrected_fastoma_inputs_20260918/proteome"
                   or Path(r["path"]).suffix != ".fa" for r in inputs)):
        raise ValueError("Wrong staged FastOMA input inventory")
    if any(r not in admission["checked_records"] for r in [pairs, *inputs]):
        raise ValueError("Native pair/input records were not admitted")
    return content["native_pair_rows"]


def convert(admission, owners, directory, mapping):
    partial = directory / "pairs.partial.tsv"
    with partial.open("x") as stream:
        counts = write_pairs(Path(admission["native_pairs"]["path"]), owners, stream, directory)
    if counts["native_rows"] != admission["content"]["native_pair_rows"]:
        raise ValueError("Native pair count changed since admission")
    filtered = directory / "pairs.qfo.partial.tsv"
    observed, retained = filter_pairs(partial, filtered, load_mapping(Path(mapping["path"])))
    if observed != counts["distinct_pairs"] or retained != observed:
        raise ValueError("Unexpected corrected-release mapping loss or count mismatch")
    return counts, partial, filtered


def prepare(root, admission_path, admission_sha, admission_job):
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong native admission allocation")
    admission = read_frozen(admission_path, admission_sha)
    validate_admission(admission, root)
    executor = root / "benchmarks/work/publication_qfo_corrected_fastoma_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
        raise ValueError("Native admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_fastoma.py"):
        raise ValueError("Wrong native admission source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen QfO mapping")
    mapping = mappings[0]
    checked = [record(__file__), record(admission_path), admission["source"], *admission["checked_records"],
               *admission["helpers"], record(env_path), mapping,
               *[record(Path(__file__).with_name(name)) for name in (
                   "fastoma_distinct_pairs.py", "fastoma_to_pairwise.py", "qfo_filter_pairs.py")]]
    for item in checked:
        check(item)
    owners = input_owners([Path(r["path"]) for r in admission["input_fastas"]])
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Changed corrected accession universe")
    directory = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1/fastoma"
    directory.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "source": record(__file__), "method": "fastoma",
              "participant": "qfo_corrected_fastoma", "semantics": SEMANTICS,
              "admission": record(admission_path), "admission_scheduler": scheduler,
              "admission_accounting": accounting, "mapping": mapping, "checked_records": checked,
              "accuracy_evaluated": False, "publication_ready": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "started_epoch": time.time()}
    try:
        counts, partial, filtered_partial = convert(admission, owners, directory, mapping)
        for item in checked:
            check(item)
        pairs, filtered = directory / "pairs.tsv", directory / "pairs.qfo.tsv"
        partial.rename(pairs)
        filtered_partial.rename(filtered)
        report.update(status="corrected_fastoma_pairs_prepared_unscored", pairs=record(pairs),
                      filtered_pairs=record(filtered), total_pairs=counts["distinct_pairs"],
                      retained_pairs=counts["distinct_pairs"], removed_mapping_pairs=0,
                      native_pair_rows=counts["native_rows"], native_duplicate_relations=counts["duplicate_relations"])
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
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job)
