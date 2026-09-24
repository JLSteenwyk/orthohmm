"""Reconstruct all available corrected comparator families from admitted raw files."""

import math
from pathlib import Path

from benchmark_tools.audit_qfo_corrected_factorial_swiss import execution_raw
from benchmark_tools.audit_qfo_corrected_swiss import audit as audit_comparator, compare
from benchmark_tools.audit_qfo_factorial_swiss import BASE_COUNTS_SHA
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.bootstrap_corrected_swiss_comparators import METHODS, validated_values
from benchmark_tools.export_qfo_complete_comparison import extract, REPLAY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def validate_inventory(comparison):
    if (comparison["status"] != "corrected_qfo_publication_comparison"
            or comparison["publication_ready"] is not False
            or tuple(row["key"] for row in comparison["methods"]) != METHODS
            or type(comparison["admitted_methods"]) is not int
            or comparison["admitted_methods"] != sum(row["status"] == "admitted" for row in comparison["methods"])):
        raise ValueError("Changed corrected comparator inventory")
    for row in comparison["methods"]:
        if row["status"] not in ("admitted", "not_admitted"):
            raise ValueError("Unknown comparator status")
        if row["status"] == "not_admitted" and (row.get("admission") is not None
                or row.get("conversion") is not None or any(v is not None for v in row["scores"].values())):
            raise ValueError("Unavailable comparator has prediction evidence")


def audit(comparison_path, comparison_sha, baseline_path):
    comparison = read_frozen(comparison_path, comparison_sha)
    validate_inventory(comparison)
    baseline = read_frozen(baseline_path, BASE_COUNTS_SHA)
    replay_record = comparison["replay_admission"]
    if replay_record["sha256"] != REPLAY_SHA:
        raise ValueError("Wrong corrected native replay admission")
    replay = read_frozen(Path(replay_record["path"]), REPLAY_SHA)
    if replay["status"] != "corrected_checked_replay_admitted":
        raise ValueError("Corrected replay was not admitted")
    checked = [record(comparison_path), record(baseline_path), replay_record,
               baseline["reference"], baseline["stages"][0]["raw_file"], *comparison["checked_records"]]
    for item in checked:
        check(item)
    anchor = read_raw(Path(baseline["stages"][0]["raw_file"]["path"]), baseline["families"])
    methods = []
    for row in comparison["methods"]:
        if row["status"] == "not_admitted":
            methods.append({"method": row["key"], "status": "not_admitted",
                            "reason": "No corrected score admission supplied in frozen comparison"})
            continue
        admission, conversion_record = row["admission"], row["conversion"]
        check(admission)
        check(conversion_record)
        report = read_frozen(Path(admission["path"]), admission["sha256"])
        conversion = read_frozen(Path(conversion_record["path"]), conversion_record["sha256"])
        if report["pairs_manifest"] != conversion_record:
            raise ValueError("Comparison conversion differs from admission")
        extracted = extract(report, conversion)
        if any(row.get(k) != v for k, v in extracted.items()):
            raise ValueError("Comparison row differs from admitted native scores")
        checked.extend([admission, conversion_record])
        if report["status"] == "corrected_factorial_assessment_admitted":
            execution_record = report["execution_report"]
            if execution_record not in report["checked_records"]:
                raise ValueError("Factorial execution was not independently checked")
            execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
            raw = execution_raw(report, execution)
            check(raw)
            result = compare(report["assessment"], read_raw(Path(raw["path"]), baseline["families"]), baseline, anchor)
            checked.extend([execution_record, raw])
        else:
            result = audit_comparator(Path(admission["path"]), admission["sha256"], baseline_path)
            if result["method_key"] != row["key"]:
                raise ValueError("Raw comparator method differs from publication row")
            raw = result["raw_file"]
            checked.extend(result["checked_records"])
        if not math.isclose(result["aggregate"]["F1"], row["scores"]["SwissTrees"], rel_tol=0, abs_tol=1e-6):
            raise ValueError("Raw SwissTrees aggregate differs from comparison")
        methods.append({"method": row["key"], "status": "counts_verified", "admission": admission,
            "raw_file": raw, "prediction_semantics": row["prediction_semantics"],
            "families": result["families"], "aggregate": result["aggregate"]})
        if report["status"] == "recovered_orthomcl_assessment_admitted":
            methods[-1].update({key: result[key] for key in (
                "search_recovery", "participant", "query_coverage", "group_coverage",
                "group_audit", "pair_semantics")})
    result = {"status": "corrected_comparison_swiss_counts_verified", "reference": baseline["reference"],
        "families": baseline["families"], "reference_relation_count": len(anchor[1]),
        "shared_represented_genes": baseline["shared_represented_genes"],
        "methods": methods, "publication_ready": False, "uncertainty_admitted": False,
        "source": record(__file__), "comparison": record(comparison_path), "checked_inputs": checked}
    validated_values(result)
    for item in checked:
        check(item)
    return result
