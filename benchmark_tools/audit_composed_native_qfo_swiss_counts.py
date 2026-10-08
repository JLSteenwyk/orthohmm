"""Audit the admitted composed-native12 Swiss raw file with unchanged count kernels."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import export_composed_native_qfo_scientific_scores as reporter
from benchmark_tools.bootstrap_qfo_factorial import validated_values
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def replay_snapshot(path, digest, evidence):
    snapshot, ref = load(path, digest, evidence)
    require(snapshot.get("schema") == "composed_native_qfo_scientific_reporting_v1"
        and snapshot.get("source") == record(reporter.__file__), "Wrong composed snapshot/source")
    replay = reporter.collect(snapshot["final_admission"])
    require(all(snapshot.get(k) == value for k, value in replay.items()), "Composed snapshot replay differs")
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"]])
    return snapshot, ref


def audit(snapshot_path, snapshot_sha, retained_path, index):
    require(type(index) is int and index == 12, "Require the composed native12 index")
    evidence = []
    snapshot, snapshot_ref = replay_snapshot(snapshot_path, snapshot_sha, evidence)
    retained, retained_ref = load(retained_path, ordinary.RETAINED_COUNTS_SHA, evidence)
    require(retained["status"] == "corrected_qfo_factorial_swiss_counts_verified"
        and retained["publication_ready"] is False and retained["uncertainty_admitted"] is False,
        "Require retained corrected family-count audit")
    validated_values({**retained, "status": "qfo_factorial_swiss_counts_verified"})
    anchor = retained["cells"][0]["raw_file"]
    evidence.extend([retained["reference"], anchor])
    check(anchor)
    families = retained["families"]
    _, truth, members = ordinary.read_raw(Path(anchor["path"]), families)
    selected = [r for r in snapshot["rows"] if r["index"] == index]
    require(len(selected) == 1 and selected[0]["status"] == "supplied_composed_native_admission"
        and selected[0]["accuracy_admitted"] is True, "Selected index has no composed admission")
    row = selected[0]
    admission, _ = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
    require(admission["execution_report"] in admission["checked_records"], "Execution was not checked")
    execution_ref = admission["execution_report"]
    execution, observed = load(execution_ref["path"], execution_ref["sha256"], evidence)
    require(observed == execution_ref, "Execution record metadata differs")
    raws = [r for r in execution["outputs"] if Path(r["path"]).parent.name == "SwissTrees"
        and r["path"].endswith("raw.txt.gz")]
    require(len(raws) == 1, "Require one inventoried SwissTrees raw file")
    raw = raws[0]
    require(raw in admission["checked_records"], "Raw file was not checked")
    check(raw)
    evidence.append(raw)
    cell = ordinary.family_cell(raw, families, truth, members, retained["reference_orientation"],
        admission["assessment"])
    require(math.isclose(cell["aggregate"]["F1"], row["scores"]["SwissTrees"], rel_tol=0, abs_tol=5e-8),
        "Reconstructed F1 differs from native endpoint")
    old = next(r for r in retained["cells"] if r["cell"] == row["cell"])
    changed = [family for family, fresh, previous in zip(families, cell["families"], old["families"])
        if fresh != previous]
    cell.update(cell=row["cell"], index=index, native_job_id=row["native_job_id"], admission=row["admission"],
        native_endpoint_f1=row["scores"]["SwissTrees"], retained_family_records_identical=not changed,
        differing_families=changed, retained_aggregate_identical=cell["aggregate"] == old["aggregate"])
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "audit_native_qfo_swiss_counts.py", "audit_qfo_swiss_counts.py", "bootstrap_qfo_factorial.py",
        "export_composed_native_qfo_scientific_scores.py", "validate_qfo_native_assessment.py")]
    source = record(__file__)
    for ref in [*evidence, *helpers, source]:
        check(ref)
    return dict(schema="composed_native_qfo_swiss_family_count_audit_v1", snapshot=snapshot_ref,
        retained_counts=retained_ref, status="supplied_composed_native_swiss_family_counts_verified",
        selected_index=index, families=families, reference=retained["reference"], cells=[cell],
        reference_relation_count=len(truth), unavailable_cells=[r["cell"] for r in snapshot["rows"]
            if r["accuracy_admitted"] is not True], count_conversion=retained["count_conversion"],
        checked_inputs=evidence, new_bootstrap_draws=0, historical_intervals_attached=False,
        new_accuracy_or_resource_admission=False, independent_confirmation=False,
        publication_ready=False, source=source, helpers=helpers, limitations=[
            "Only the selected composed-native12 cell is recounted; previous cell audits are reused separately.",
            "Uses unchanged raw-count/native-endpoint kernels and complete retained reference labels/members.",
            "Mismatched native family records are retained, not replaced with cached records or intervals.",
            "No new inference, scoring, FAS sampling, bootstrap draws or resource admission.",
            "Eighteen development-exposed SwissTrees families are not independent confirmation."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "retained-counts", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = audit(args.snapshot, args.snapshot_sha256, args.retained_counts, args.index)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(cells=1, family_records=len(result["families"]), selected_index=args.index)))


if __name__ == "__main__":
    main()
