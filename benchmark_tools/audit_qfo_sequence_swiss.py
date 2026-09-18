"""Bind corrected HMM and sequence-control admissions to exact SwissTrees counts."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_factorial_swiss import BASE_COUNTS_SHA, file_records
from benchmark_tools.audit_qfo_corrected_factorial_swiss import execution_raw
from benchmark_tools.audit_qfo_swiss_counts import read_raw, verify_native, REFERENCE_SHA
from benchmark_tools.bootstrap_qfo_sequence import VARIANTS, validated_values
from benchmark_tools.export_qfo_corrected_factorial import extract
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_sequence_assessment import validate_stage
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_records


def selected_raw(report, execution, conversion, variant):
    if variant == "p0_c0_r0":
        row = extract(report, conversion)
        if row["index"] != 0 or row["cell"] != variant:
            raise ValueError("Require corrected initial HMM without P/C/R")
        return execution_raw(report, execution)
    if (variant not in VARIANTS[1:] or report.get("status") != "corrected_sequence_assessment_admitted"
            or report.get("variant") != variant or report.get("accuracy_admitted") is not True
            or report.get("publication_ready") is not False or report["conversion"] != conversion):
        raise ValueError("Wrong or unadmitted corrected sequence assessment")
    validate_stage(conversion, variant, report["conversion_scheduler"])
    if report["assessment"]["participant"] != conversion["participant"]:
        raise ValueError("Wrong assessment participant")
    scheduler = report["scheduler"]
    if ((scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"])
            != ("COMPLETED", "0:0", "bizon", "8")
            or execution["status"] != "process_succeeded_pending_independent_admission"
            or type(execution["exit_code"]) is not int or execution["exit_code"] != 0
            or execution["job_id"] != scheduler["JobIDRaw"] or execution["variant"] != variant
            or execution["stage"] != conversion or execution["pairs_manifest"] != report["pairs_manifest"]
            or execution["conversion_scheduler"] != report["conversion_scheduler"]):
        raise ValueError("Sequence execution contradicts admission")
    paths = [item for item in execution["outputs"] if Path(item["path"]).parent.name == "SwissTrees"
             and item["path"].endswith("raw.txt.gz")]
    if len(paths) != 1:
        raise ValueError("Require exactly one inventoried SwissTrees raw file")
    return paths[0]


def assemble(entries, baseline):
    if [entry["variant"] for entry in entries] != list(VARIANTS):
        raise ValueError("Require all three predictions in frozen order")
    if baseline["status"] != "raw_swiss_family_counts_verified" or baseline["reference"]["sha256"] != REFERENCE_SHA:
        raise ValueError("Wrong reference audit")
    families = baseline["families"]
    _, expected_truth, expected_members = read_raw(Path(baseline["stages"][0]["raw_file"]["path"]), families)
    orientation, variants = baseline["reference_orientation"], []
    for entry in entries:
        assessment = entry["assessment"]
        if assessment["swiss_reference_families"] != families:
            raise ValueError("Changed native family inventory")
        native = validate_records(assessment["native_assessments"], assessment["participant"], set(families))
        counts, truth, members = read_raw(Path(entry["raw_file"]["path"]), families)
        if truth != expected_truth or members != expected_members:
            raise ValueError("Reference pair identities, labels or members differ")
        if any(sum(counts[f].values()) != orientation[f]["forward_relations"] or
               len(members[f]) != orientation[f]["mapped_proteins"] for f in families):
            raise ValueError("Incomplete reference coverage")
        values, aggregate = verify_native(counts, native)
        score = assessment["endpoints"]["SwissTrees"]["score"]
        if type(score) not in (int, float) or not math.isclose(score, aggregate["F1"], rel_tol=0, abs_tol=1e-6):
            raise ValueError("Reconstructed F1 differs from admitted score")
        variants.append({"variant": entry["variant"], "raw_file": entry["raw_file"], "aggregate": aggregate,
            "families": [{"family": f, "counts_without_prior": dict(counts[f]),
                "statistics_with_prior": values[f], "represented_genes": sorted(members[f])} for f in families]})
    result = dict(status="corrected_qfo_sequence_swiss_counts_verified", publication_ready=False,
        uncertainty_admitted=False, reference=baseline["reference"], families=families, variants=variants,
        reference_relation_count=len(expected_truth), reference_orientation=orientation,
        shared_represented_genes=baseline["shared_represented_genes"],
        count_conversion="native confusion count = raw one-direction relation count / 2 + 1",
        limitations=["Historical evidence anchors reference labels and members only, never prediction counts.",
            "Requires independently admitted corrected assessments; does not rerun inference or scoring.",
            "Development-exposed SwissTrees evidence only; no other endpoint uncertainty is established."])
    validated_values(result)
    return result


def audit(inventory_path, inventory_sha, baseline_path):
    inventory = read_frozen(inventory_path, inventory_sha)
    baseline = read_frozen(baseline_path, BASE_COUNTS_SHA)
    if [row["variant"] for row in inventory["variants"]] != list(VARIANTS):
        raise ValueError("Require complete ordered admission inventory")
    checked = {}

    def verify(items):
        for item in items:
            previous = checked.get(item["path"])
            if previous is not None and previous != item:
                raise ValueError("Conflicting file identities")
            if previous is None:
                check(item)
                checked[item["path"]] = item

    verify([record(inventory_path), record(baseline_path), baseline["reference"],
            baseline["stages"][0]["raw_file"]])
    entries = []
    for entry in inventory["variants"]:
        verify([entry["admission"]])
        report = read_frozen(Path(entry["admission"]["path"]), entry["admission"]["sha256"])
        verify(file_records(report))
        execution_record, pair_record = report["execution_report"], report["pairs_manifest"]
        if execution_record not in report["checked_records"]:
            raise ValueError("Execution was not independently checked")
        execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
        conversion = read_frozen(Path(pair_record["path"]), pair_record["sha256"])
        verify(file_records(execution))
        verify(file_records(conversion))
        raw = selected_raw(report, execution, conversion, entry["variant"])
        verify([raw])
        entries.append({"variant": entry["variant"], "assessment": report["assessment"], "raw_file": raw})
    result = assemble(entries, baseline)
    for item in checked.values():
        check(item)
    result.update(admission_inventory=record(inventory_path), baseline_audit=record(baseline_path),
        checked_inputs=list(checked.values()), source=record(__file__),
        helpers=[record(Path(__file__).with_name(name)) for name in (
            "audit_qfo_factorial_swiss.py", "audit_qfo_corrected_factorial_swiss.py", "audit_qfo_swiss_counts.py",
            "bootstrap_qfo_sequence.py", "bootstrap_qfo_swiss_stages.py", "export_qfo_corrected_factorial.py",
            "run_qfo_sequence_assessment.py", "run_qfo_corrected_factorial_assessment.py",
            "validate_qfo_native_assessment.py")])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inventory", "baseline", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.inventory.resolve(), args.inventory_sha256, args.baseline.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
