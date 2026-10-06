"""Reuse frozen SwissTrees intervals only for identical admitted native counts."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import audit_recovered_native_qfo_swiss_counts as recovered
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.bootstrap_qfo_factorial import CELLS, contrasts, validated_values
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

BOOTSTRAP_SHA = "f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211"
AUDITS = {
    "native_qfo_swiss_family_count_audit_v1": (
        "supplied_native_swiss_family_counts_verified", ordinary, "supplied_native_admission"),
    "recovered_native_qfo_swiss_family_count_audit_v1": (
        "supplied_recovered_native_swiss_family_counts_verified", recovered,
        "supplied_recovered_scientific_admission"),
}


def project(native, retained, bootstrap):
    values = validated_values({**retained, "status": "qfo_factorial_swiss_counts_verified"})
    require(bootstrap["status"] == "paired_corrected_qfo_factorial_swiss_intervals"
            and bootstrap["publication_ready"] is False and bootstrap["replicates"] == 100000
            and bootstrap["seed"] == 20260922 and bootstrap["alpha"] == .05
            and bootstrap["multiplicity_endpoints"] == 42 and bootstrap["quantile_method"] == "linear"
            and bootstrap["families"] == retained["families"], "Changed frozen bootstrap scope")
    points = aggregate(values.mean(axis=1))
    require(bootstrap["point_estimates"] == {
        cell: dict(zip(METRICS, row.tolist())) for cell, row in zip(CELLS, points)},
        "Retained bootstrap point estimates differ from counts")
    expected = contrasts()
    require(len(bootstrap["comparisons"]) == len(expected), "Changed contrast inventory")
    history = {row["cell"]: row for row in retained["cells"]}
    matched, bindings = set(), {}
    for cell, row in native.items():
        require(cell in history, "Unknown native cell")
        same = row["families"] == history[cell]["families"] and row["aggregate"] == history[cell]["aggregate"]
        require(row["retained_family_records_identical"] is (row["families"] == history[cell]["families"])
                and row["retained_aggregate_identical"] is (row["aggregate"] == history[cell]["aggregate"]),
                "Count audit equality flags differ from actual records")
        if same:
            matched.add(cell)
        bindings[cell] = dict(status="native_records_matched" if same else "native_records_differ",
            index=row["index"], native_job_id=row["native_job_id"], admission=row["admission"],
            count_audit=row["count_audit"], native_endpoint_f1=row["native_endpoint_f1"],
            count_aggregate=row["aggregate"])
    rows = []
    for actual, definition in zip(bootstrap["comparisons"], expected):
        require(all(actual[k] == value for k, value in definition.items()), "Changed contrast definition")
        require(set(actual["metrics"]) == set(METRICS), "Incomplete retained metrics")
        weights = np.asarray(definition["weights"])
        point = weights @ points
        require(all(actual["metrics"][m]["difference"] == float(point[j]) for j, m in enumerate(METRICS)),
                "Retained contrast arithmetic differs")
        needed = [cell for cell, weight in zip(CELLS, weights) if weight]
        absent = [cell for cell in needed if cell not in native]
        different = [cell for cell in needed if cell in native and cell not in matched]
        available = not absent and not different
        rows.append(dict(definition, status="native_records_matched" if available else "native_records_unavailable",
            missing_cells=absent, differing_cells=different,
            metrics=actual["metrics"] if available else None,
            family_differences=actual["family_differences"] if available else None))
    return bindings, rows


def bind(snapshot_path, snapshot_sha, audits, retained_path, bootstrap_path):
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["source"] == record(reporter.__file__), "Wrong scientific snapshot/source")
    normal, repairs = [], []
    for row in snapshot["rows"]:
        destination = (normal if row["status"] == "supplied_native_admission" else repairs
                       if row["status"] == "supplied_recovered_scientific_admission" else None)
        if destination is not None:
            destination.append((row["admission"]["path"], row["admission"]["sha256"]))
    replay = reporter.collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], normal, repairs)
    require(all(snapshot[k] == value for k, value in replay.items()), "Scientific snapshot replay differs")
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"]])
    retained, retained_ref = load(retained_path, ordinary.RETAINED_COUNTS_SHA, evidence)
    bootstrap, bootstrap_ref = load(bootstrap_path, BOOTSTRAP_SHA, evidence)
    require(retained["status"] == "corrected_qfo_factorial_swiss_counts_verified"
            and retained["publication_ready"] is False and retained["uncertainty_admitted"] is False,
            "Require retained corrected family counts")
    original_counts, observed = load(bootstrap["counts"]["path"], bootstrap["counts"]["sha256"], evidence)
    require(observed == bootstrap["counts"] and original_counts == retained, "Bootstrap count binding differs")
    evidence.extend([bootstrap["source"], bootstrap["protocol"], bootstrap["corrected_protocol"],
                     *bootstrap["helpers"], retained["reference"]])
    snapshot_rows = {row["cell"]: row for row in snapshot["rows"]}
    native = {}
    require(audits, "No supplied native count audit")
    for path, digest in audits:
        audit, ref = load(path, digest, evidence)
        require(audit.get("schema") in AUDITS, "Unknown count audit schema")
        status, source, row_status = AUDITS[audit["schema"]]
        require(audit["status"] == status and audit["source"] == record(source.__file__)
                and audit["families"] == retained["families"] and audit["retained_counts"] == retained_ref
                and audit["reference"] == retained["reference"]
                and audit["reference_relation_count"] == retained["reference_relation_count"]
                and audit["count_conversion"] == retained["count_conversion"]
                and type(audit["new_bootstrap_draws"]) is int and audit["new_bootstrap_draws"] == 0
                and all(audit[k] is False for k in ("historical_intervals_attached",
                    "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
                "Count audit source, scope or reference differs")
        evidence.extend([*audit["checked_inputs"], *audit["helpers"], audit["source"], audit["snapshot"]])
        for row in audit["cells"]:
            cell = row["cell"]
            require(cell not in native and cell in snapshot_rows, "Duplicate or unknown audited cell")
            target = snapshot_rows[cell]
            require(target["status"] == row_status and row["admission"] == target["admission"]
                    and row["index"] == target["index"] and row["native_job_id"] == target["native_job_id"]
                    and row["native_endpoint_f1"] == target["scores"]["SwissTrees"]
                    and row["raw_file"] in audit["checked_inputs"], "Audit/admission binding differs")
            if row_status == "supplied_recovered_scientific_admission":
                require(row["resources"] is None and row["timing_admitted"] is False
                        and row["timing_eligible"] is False, "Count audit relabels failed timing")
            native[cell] = dict(row, count_audit=ref)
    bindings, rows = project(native, retained, bootstrap)
    for ref in evidence:
        check(ref)
    return dict(schema="native_qfo_retained_swiss_uncertainty_binding_v1", snapshot=snapshot_ref,
        retained_counts=retained_ref, bootstrap=bootstrap_ref, bound_cells=bindings, contrasts=rows,
        families=retained["families"], replicates_reused=bootstrap["replicates"], seed_reused=bootstrap["seed"],
        alpha=.05, multiplicity_endpoints=42, new_bootstrap_draws=0, independent_confirmation=False,
        new_accuracy_or_resource_admission=False, publication_ready=False, evidence=evidence,
        source=record(__file__), limitations=[
            "Intervals reused only for identical audited native family records and aggregates, not new independent evidence.",
            "Unmatched or missing cells remain unavailable; aggregate agreement alone is insufficient.",
            "All 42 planned endpoints remain in the adjustment, even when native contrasts are unavailable.",
            "Development exposure, approximate percentile coverage and family exchangeability limits remain.",
            "SwissTrees only; no intervals for other QfO endpoints or the secondary mean.",
            "P-off retains initial HMM search; R changes clique to resolved-pair prediction semantics.",
            "Failed native timing remains ineligible; no resource admission, inference/scoring or bootstrap rerun.",
            "Direct report/source binding does not repeat raw count audits or transitive scientific admission."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "retained-counts", "bootstrap", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--counts-audit", nargs=2, action="append", required=True, metavar=("PATH", "SHA256"))
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = bind(args.snapshot, args.snapshot_sha256, args.counts_audit, args.retained_counts, args.bootstrap)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(bound_cells=len(result["bound_cells"]), matched_contrasts=sum(
        row["status"] == "native_records_matched" for row in result["contrasts"]))))


if __name__ == "__main__":
    main()
