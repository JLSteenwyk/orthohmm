"""Bind composed-native12 and prior audits to unchanged retained Swiss intervals."""

import argparse
import json
from pathlib import Path

from benchmark_tools import audit_composed_native_qfo_swiss_counts as composed
from benchmark_tools import bind_allocated_native_qfo_swiss_uncertainty as allocated
from benchmark_tools import bind_native_qfo_swiss_uncertainty as original
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

AUDITS = {**allocated.AUDITS, "composed_native_qfo_swiss_family_count_audit_v1": (
    "supplied_composed_native_swiss_family_counts_verified", composed, "supplied_composed_native_admission")}


def bind(snapshot_path, snapshot_sha, audits, retained_path, bootstrap_path):
    evidence = []
    snapshot, snapshot_ref = composed.replay_snapshot(snapshot_path, snapshot_sha, evidence)
    retained, retained_ref = load(retained_path, original.ordinary.RETAINED_COUNTS_SHA, evidence)
    bootstrap, bootstrap_ref = load(bootstrap_path, original.BOOTSTRAP_SHA, evidence)
    require(retained["status"] == "corrected_qfo_factorial_swiss_counts_verified"
        and retained["publication_ready"] is False and retained["uncertainty_admitted"] is False,
        "Require retained corrected family counts")
    original_counts, observed = load(bootstrap["counts"]["path"], bootstrap["counts"]["sha256"], evidence)
    require(observed == bootstrap["counts"] and original_counts == retained, "Bootstrap count binding differs")
    evidence.extend([bootstrap["source"], bootstrap["protocol"], bootstrap["corrected_protocol"],
        *bootstrap["helpers"], retained["reference"]])
    targets = {r["cell"]: r for r in snapshot["rows"]}
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
        if row_status == "supplied_allocated_native_admission":
            require(len(audit["cells"]) == 1 and type(audit["selected_index"]) is int
                and audit["selected_index"] in (10, 11, 12), "Invalid selected allocated audit")
        if row_status == "supplied_composed_native_admission":
            require(len(audit["cells"]) == 1 and type(audit["selected_index"]) is int
                and audit["selected_index"] == 12, "Invalid selected composed audit")
        for row in audit["cells"]:
            cell = row["cell"]
            require(cell not in native and cell in targets, "Duplicate or unknown audited cell")
            target = targets[cell]
            require(target["status"] == row_status and target["accuracy_admitted"] is True
                and row["admission"] == target["admission"] and row["index"] == target["index"]
                and row["native_job_id"] == target["native_job_id"]
                and row["native_endpoint_f1"] == target["scores"]["SwissTrees"]
                and row["raw_file"] in audit["checked_inputs"], "Audit/admission binding differs")
            if row_status == "supplied_allocated_native_admission":
                require(row["index"] == audit["selected_index"], "Allocated selected index differs")
            if row_status == "supplied_composed_native_admission":
                require(type(row["index"]) is int and row["index"] == audit["selected_index"] == 12,
                    "Composed selected index differs")
            if row_status == "supplied_recovered_scientific_admission":
                require(row["resources"] is None and row["timing_admitted"] is False
                    and row["timing_eligible"] is False, "Count audit relabels failed timing")
            native[cell] = dict(row, count_audit=ref)
    bindings, rows = original.project(native, retained, bootstrap)
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "audit_composed_native_qfo_swiss_counts.py", "bind_allocated_native_qfo_swiss_uncertainty.py",
        "bind_native_qfo_swiss_uncertainty.py",
        "bootstrap_qfo_factorial.py", "bootstrap_qfo_swiss_stages.py")]
    source = record(__file__)
    for ref in [*evidence, *helpers, source]:
        check(ref)
    return dict(schema="composed_native_qfo_retained_swiss_uncertainty_binding_v1", snapshot=snapshot_ref,
        retained_counts=retained_ref, bootstrap=bootstrap_ref, bound_cells=bindings, contrasts=rows,
        families=retained["families"], replicates_reused=bootstrap["replicates"], seed_reused=bootstrap["seed"],
        alpha=.05, multiplicity_endpoints=42, new_bootstrap_draws=0, independent_confirmation=False,
        new_accuracy_or_resource_admission=False, publication_ready=False, evidence=evidence,
        source=source, helpers=helpers, limitations=[
            "Unchanged interval projector; complete identical native family records required, not aggregate agreement alone.",
            "All 42 endpoints remain in adjustment; unavailable/differing contrasts retain null metrics.",
            "Development exposure, family exchangeability and approximate percentile coverage limits remain.",
            "SwissTrees only; no uncertainty for other QfO endpoints or the secondary mean.",
            "P-off retains initial HMM search; R changes clique to resolved-pair semantics.",
            "Previous raw audits are reused without recount; failed timing remains ineligible.",
            "No new raw scientific admission, inference, scoring, bootstrap or isolated efficiency claim."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "retained-counts", "bootstrap", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--counts-audit", nargs=2, action="append", required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = bind(args.snapshot, args.snapshot_sha256, args.counts_audit, args.retained_counts, args.bootstrap)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(bound_cells=len(result["bound_cells"]), matched_contrasts=sum(
        r["status"] == "native_records_matched" for r in result["contrasts"]))))


if __name__ == "__main__":
    main()

