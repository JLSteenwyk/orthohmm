"""Estimate both admitted corrected-QfO graph inputs without launching inference."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed, validate_conversion, CONVERTER
from benchmark_tools.compare_qfo_search_coverage import frozen, validate_numeric, NUMERIC_ADMITTER
from benchmark_tools.estimate_rbnh_array_payload import estimate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def estimate_variants(conversion, core):
    if set(conversion["variants"]) != {"all_hits", "top100"}:
        raise ValueError("Require both prespecified variants")
    results = {}
    for label in ("all_hits", "top100"):
        variant = conversion["variants"][label]
        manifest = variant["manifest"]
        check(manifest)
        checkpoint = Path(variant["checkpoint"])
        if str(checkpoint / "manifest.json") != manifest["path"]:
            raise ValueError("Checkpoint location differs")
        result = estimate(checkpoint, manifest["sha256"], core)
        expected = variant["audit"]["summary"]
        if any(result["estimate"][key] != expected[key] for key in ("genes", "hits", "self_hits")):
            raise ValueError("Estimator counts differ from admitted checkpoint")
        if result["estimate"]["graph_feasibility_admitted"] is not False:
            raise ValueError("Payload estimate cannot admit graph feasibility")
        results[label] = result
    return results


def run(root, admission_job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    admission_scheduler = completed(admission_job, 2, "192G")
    executor = frozen(root, "publication_qfo_sequence_numeric_admission_v1", NUMERIC_ADMITTER)
    converter = frozen(root, "publication_qfo_sequence_numeric_v1", CONVERTER)
    admission_record = record(root / "benchmarks/work/qfo_sequence_numeric_admission_20260918/report.json")
    admission = read_frozen(Path(admission_record["path"]), admission_record["sha256"])
    conversion_record = admission["conversion"]
    conversion_path = root / "benchmarks/results/qfo_sequence_numeric_v1/manifest.json"
    if conversion_record["path"] != str(conversion_path):
        raise ValueError("Unexpected conversion path")
    conversion = read_frozen(conversion_path, conversion_record["sha256"])
    scheduler = completed(admission["scheduler"]["JobIDRaw"], 2, "192G")
    if scheduler != admission["scheduler"]:
        raise ValueError("Conversion accounting changed")
    validate_conversion(conversion, scheduler,
                        record(converter / "benchmark_tools/convert_qfo_sequence_search_control.py"),
                        conversion_path.parent)
    validate_numeric(admission, conversion, conversion_record,
                     record(executor / "benchmark_tools/admit_qfo_sequence_numeric.py"))
    checked = [admission_record, conversion_record, admission["source"],
               *admission["helpers"], *admission["checked_records"]]
    for item in checked:
        check(item)
    core = root / "benchmarks/work/publication_method_native_v2/orthohmm/accuracy.py"
    result = dict(status="admitted_qfo_graph_payload_estimated", job_id=os.environ.get("SLURM_JOB_ID"),
                  admission_scheduler=admission_scheduler, conversion_scheduler=scheduler,
                  variants=estimate_variants(conversion, core), checked_records=checked,
                  source=record(__file__),
                  helpers=[record(Path(__file__).with_name(name)) for name in (
                      "estimate_rbnh_array_payload.py", "audit_accuracy_checkpoint.py",
                      "admit_qfo_sequence_numeric.py", "compare_qfo_search_coverage.py")],
                  graph_launched=False, graph_feasibility_admitted=False,
                  accuracy_evaluated=False, publication_ready=False,
                  limitations=["Bounds one named-array snapshot, not total peak memory.",
                               "Manual resource review is required before graph submission.",
                               "All-hit failure must not be replaced silently by the top100 diagnostic."])
    for item in [*checked, result["source"], *result["helpers"]]:
        check(item)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--admission-job", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.admission_job, args.output.absolute())
