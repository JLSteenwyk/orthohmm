"""Verify corrected comparator SwissTrees raw counts, without inference tests."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_factorial_swiss import BASE_COUNTS_SHA
from benchmark_tools.audit_qfo_swiss_counts import read_raw, verify_native, REFERENCE_SHA
from benchmark_tools.export_qfo_corrected_comparison import extract
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_records


def compare(assessment, raw, baseline, anchor):
    families = baseline["families"]
    if (baseline["status"] != "raw_swiss_family_counts_verified"
            or baseline["reference"]["sha256"] != REFERENCE_SHA
            or len(families) != 18 or len(set(families)) != 18
            or assessment["swiss_reference_families"] != families):
        raise ValueError("Changed SwissTrees reference inventory")
    counts, truth, members = raw
    _, expected_truth, expected_members = anchor
    if truth != expected_truth or members != expected_members:
        raise ValueError("Changed reference relation identities, labels or members")
    orientation = baseline["reference_orientation"]
    if set(counts) != set(families) or any(
            sum(counts[f].values()) != orientation[f]["forward_relations"]
            or len(members[f]) != orientation[f]["mapped_proteins"] for f in families):
        raise ValueError("Incomplete raw reference coverage")
    native = validate_records(assessment["native_assessments"], assessment["participant"], set(families))
    values, aggregate = verify_native(counts, native)
    return {"aggregate": aggregate, "reference_relation_count": len(truth),
            "families": [{"family": f, "counts_without_prior": dict(counts[f]),
                          "statistics_with_prior": values[f], "represented_genes": sorted(members[f])}
                         for f in families]}


def audit(admission_path, admission_sha, baseline_path):
    report = read_frozen(admission_path, admission_sha)
    baseline = read_frozen(baseline_path, BASE_COUNTS_SHA)
    checked = [record(admission_path), record(baseline_path), report["pairs_manifest"],
               baseline["reference"], baseline["stages"][0]["raw_file"]]
    for item in checked:
        check(item)
    conversion = json.loads(Path(report["pairs_manifest"]["path"]).read_text())
    row = extract(report, conversion)
    # The independent admission pins the execution report, which inventories raw outputs.
    recovered = report["status"] == "recovered_orthomcl_assessment_admitted"
    suffix = ("/qfo_blast_recovery_assessment_v1/results.json" if recovered else
              f"/qfo_corrected_assessment_v1/{report['method']}/results.json")
    execution_records = [r for r in report["checked_records"] if r["path"].endswith(suffix)]
    if len(execution_records) != 1:
        raise ValueError("Missing or ambiguous admitted execution record")
    execution_record = execution_records[0]
    if recovered and execution_record != report["execution_report"]:
        raise ValueError("Recovered execution record contradicts independent admission")
    check(execution_record)
    execution = json.loads(Path(execution_record["path"]).read_text())
    if (execution["status"] != "process_succeeded_pending_independent_admission"
            or execution["exit_code"] != 0 or execution["method"] != report["method"]
            or execution["job_id"] != report["scheduler"]["JobIDRaw"]
            or report["scheduler"]["State"] != "COMPLETED"
            or report["scheduler"]["ExitCode"] != "0:0"
            or execution["pairs_manifest"] != report["pairs_manifest"]):
        raise ValueError("Execution contradicts independent admission")
    if recovered and execution["stage"] != conversion:
        raise ValueError("Recovered execution stage contradicts admitted conversion")
    raw_records = [r for r in execution["outputs"] if
                   Path(r["path"]).parent.name == "SwissTrees" and r["path"].endswith("raw.txt.gz")]
    if len(raw_records) != 1:
        raise ValueError("Missing or ambiguous inventoried raw evidence")
    raw_record = raw_records[0]
    check(raw_record)
    checked.extend([execution_record, raw_record])
    families = baseline["families"]
    result = compare(report["assessment"], read_raw(Path(raw_record["path"]), families), baseline,
                     read_raw(Path(baseline["stages"][0]["raw_file"]["path"]), families))
    for item in checked:
        check(item)
    result.update(status="corrected_comparator_swiss_counts_verified", method=report["method"],
                  method_key=row["key"], participant=report["assessment"]["participant"],
                  reference=baseline["reference"], reference_orientation=baseline["reference_orientation"],
                  shared_represented_genes=baseline["shared_represented_genes"],
                  admission=record(admission_path), raw_file=raw_record, checked_records=checked,
                  source=record(__file__), helpers=[record(Path(__file__).with_name(n)) for n in (
                      "audit_qfo_swiss_counts.py", "export_qfo_corrected_comparison.py",
                      "validate_qfo_native_assessment.py", "prepare_ob_candidate_neighborhood.py",
                      "run_simulation_methods.py", "audit_qfo_factorial_swiss.py")],
                  publication_ready=False, uncertainty_admitted=False,
                  count_conversion="native confusion count = raw one-direction relation count / 2 + 1",
                  limitations=["Count audit only; no paired intervals, rankings or release-effect tests.",
                      "Historical raw data anchor reference identities and labels only; corrected counts are read afresh.",
                      "Relies on the pinned independent admission for inference and full scoring provenance; does not rerun them.",
                      "Full prespecified eight-method analysis remains required; these are development-exposed families."])
    if recovered:
        result.update(search_recovery=True, query_coverage=report["query_coverage"],
                      group_coverage=report["group_coverage"], group_audit=report["group_audit"],
                      pair_semantics=report["pair_semantics"])
        result["limitations"].append(
            "Recovered BLAST search, not an uninterrupted run; failed-query and group coverage remain in this audit.")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("admission", "baseline", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.admission.resolve(), args.admission_sha256, args.baseline.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
