"""Compare admitted corrected QfO search hits without benchmark reference labels."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed, validate_conversion, CONVERTER
from benchmark_tools.canonicalize_search_checkpoint import canonicalize
from benchmark_tools.compare_search_hit_coverage import same_species_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_replay import validate_admission, PLAN_SHA
from benchmark_tools.run_qfo_sequence_search_control import verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.summarize_sorted_search_hits import summarize, overlap
from orthohmm.accuracy import load_accuracy_checkpoint

HMM_ADMITTER = "7b5214a5d2169338c4cbe80293fae491ee5b951f"
NUMERIC_ADMITTER = "617588f0f0acef0134b3e5fe18480e00eaa22129"


def frozen(root, name, commit):
    path = root / "benchmarks/work" / name
    if subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Frozen admission revision differs")
    subprocess.run(["git", "-C", str(path), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    return path


def validate_numeric(admission, conversion, conversion_record, source):
    if (admission["status"] != "corrected_qfo_numeric_source_equivalence_admitted"
            or admission["numeric_equivalence"] is not True or admission["accuracy_evaluated"] is not False
            or admission["publication_ready"] is not False or admission["source"] != source
            or admission["conversion"] != conversion_record or admission["genes"] != 984137
            or admission["proteomes"] != 78 or admission["hits"] != conversion["hits"]
            or set(admission["variants"]) != {"all_hits", "top100"}):
        raise ValueError("Wrong corrected numeric equivalence admission")
    for label, cap in (("all_hits", None), ("top100", 100)):
        variant = admission["variants"][label]
        if (variant["status"] != "checkpoint_matches_reconstructed_source_hits"
                or variant["cap"] != cap or type(variant["cap"]) is not type(cap)
                or variant["genes"] != 984137 or type(variant["hits"]) is not int
                or variant["hits"] <= 0 or variant["accuracy_evaluated"] is not False
                or variant["hits"] != conversion["variants"][label]["audit"]["summary"]["hits"]):
            raise ValueError("Wrong reconstructed variant identity/count")


def compare(checkpoints, metadata):
    if set(checkpoints) != {"hmm", "all_hits", "top100"}:
        raise ValueError("Require exactly three search checkpoints")
    names = sorted(metadata)
    labels = sorted({row["species"] for row in metadata.values()})
    ids = {name: i for i, name in enumerate(labels)}
    owners = np.asarray([ids[metadata[name]["species"]] for name in names], dtype=np.int32)
    loaded, reports = {}, {}
    for label in ("hmm", "all_hits", "top100"):
        manifest = checkpoints[label]
        check(manifest)
        actual_names, species, q, t, s = load_accuracy_checkpoint(Path(manifest["path"]).parent, verify=True)
        if actual_names != names:
            raise ValueError("Search gene universe/order differs")
        same_species_partition(species, owners)
        loaded[label] = (q, t, s)
        reports[label] = summarize(q, t, s, owners)
    intersections = {a + "_vs_" + b: overlap(loaded[a], loaded[b], len(names))
                     for a, b in (("hmm", "all_hits"), ("hmm", "top100"), ("all_hits", "top100"))}
    if intersections["all_hits_vs_top100"]["all"]["second_only"]:
        raise ValueError("Top100 has hits absent from all-hit checkpoint")
    for item in checkpoints.values():
        check(item)
        load_accuracy_checkpoint(Path(item["path"]).parent, verify=True)
    return {"species_labels": labels, "searches": reports, "overlaps": intersections}


def run(root, hmm_job, numeric_job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    schedulers = {"hmm_admission": completed(hmm_job, 2, "64G"),
                  "numeric_admission": completed(numeric_job, 2, "192G")}
    hmm_executor = frozen(root, "publication_qfo_corrected_hmm_admission_v1", HMM_ADMITTER)
    numeric_executor = frozen(root, "publication_qfo_sequence_numeric_admission_v1", NUMERIC_ADMITTER)
    converter = frozen(root, "publication_qfo_sequence_numeric_v1", CONVERTER)
    hmm_record = record(root / "benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json")
    numeric_record = record(root / "benchmarks/work/qfo_sequence_numeric_admission_20260918/report.json")
    hmm = read_frozen(Path(hmm_record["path"]), hmm_record["sha256"])
    numeric = read_frozen(Path(numeric_record["path"]), numeric_record["sha256"])
    if hmm["source"] != record(hmm_executor / "benchmark_tools/admit_qfo_corrected_high_sensitivity.py"):
        raise ValueError("HMM admission source differs")
    native_scheduler = completed(hmm["scheduler"]["JobIDRaw"], 32, "192G")
    if any(native_scheduler[k] != v for k, v in hmm["scheduler"].items()):
        raise ValueError("HMM accounting changed")
    plan = verify_plan(root / "benchmarks/work/qfo_sequence_search_control_v1/manifest.json")
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PLAN_SHA)
    checkpoint, checkpoint_sha = validate_admission(hmm, primary, plan["inputs"])
    conversion_record = numeric["conversion"]
    expected_conversion = root / "benchmarks/results/qfo_sequence_numeric_v1/manifest.json"
    if conversion_record["path"] != str(expected_conversion):
        raise ValueError("Wrong numeric conversion location")
    check(conversion_record)
    conversion = read_frozen(expected_conversion, conversion_record["sha256"])
    scheduler = completed(numeric["scheduler"]["JobIDRaw"], 2, "192G")
    if scheduler != numeric["scheduler"]:
        raise ValueError("Conversion accounting changed")
    validate_conversion(conversion, scheduler, record(converter / "benchmark_tools/convert_qfo_sequence_search_control.py"), expected_conversion.parent)
    validate_numeric(numeric, conversion, conversion_record,
                     record(numeric_executor / "benchmark_tools/admit_qfo_sequence_numeric.py"))
    checked = [hmm_record, numeric_record, conversion_record, record(primary_path),
        *hmm["checked_records"], *hmm["content"]["checked_records"], *numeric["checked_records"],
        hmm["source"], numeric["source"], *numeric["helpers"], plan["gene_metadata"], *plan["inputs"]]
    for item in checked:
        check(item)
    metadata = read_frozen(Path(plan["gene_metadata"]["path"]), plan["gene_metadata"]["sha256"])
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "comparing", "source": record(__file__), "checked_records": checked,
        "schedulers": schedulers, "accuracy_evaluated": False, "publication_ready": False,
        "helpers": [record(Path(__file__).with_name(name)) for name in (
            "canonicalize_search_checkpoint.py", "summarize_sorted_search_hits.py",
            "admit_qfo_sequence_numeric.py", "prepare_qfo_corrected_replay.py")],
        "numpy_version": np.__version__}
    try:
        canonical = canonicalize(checkpoint, checkpoint_sha, output / "hmm_sorted")
        checkpoints = {"hmm": canonical["manifest"], **{k: v["manifest"] for k, v in conversion["variants"].items()}}
        result.update(compare(checkpoints, metadata))
        for item in [*checked, result["source"], *result["helpers"]]:
            check(item)
        result.update(status="corrected_qfo_label_free_hit_comparison_complete", checkpoints=checkpoints,
            canonicalization=record(output / "hmm_sorted/manifest.json"),
            limitations=["Hit overlap is not sensitivity against biological truth.",
                "Equal E-value cutoffs do not imply equal calibration, sensitivity or computation.",
                "Top100 is a reporting diagnostic, not an emulation of the HMM prefilter.",
                "Shared-host diagnostic execution is not dedicated inference timing.",
                "Graph replay and scientific accuracy evaluation remain separate."])
    except BaseException as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("hmm-job", "numeric-job"):
        parser.add_argument("--" + flag, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.hmm_job, args.numeric_job, args.output.absolute())
