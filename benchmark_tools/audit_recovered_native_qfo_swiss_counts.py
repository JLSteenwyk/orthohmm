"""Verify recovered native SwissTrees counts without admitting failed timing."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.bootstrap_qfo_factorial import validated_values
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def audit(snapshot_path, snapshot_sha, retained_path):
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["source"] == record(reporter.__file__), "Wrong scientific snapshot/source")
    normal, recovered = [], []
    for row in snapshot["rows"]:
        destination = (normal if row["status"] == "supplied_native_admission" else recovered
                       if row["status"] == "supplied_recovered_scientific_admission" else None)
        if destination is not None:
            destination.append((row["admission"]["path"], row["admission"]["sha256"]))
    require(recovered, "No supplied recovered accuracy admission")
    replay = reporter.collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], normal, recovered)
    require(all(snapshot[k] == value for k, value in replay.items()), "Scientific snapshot replay differs")
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"]])
    retained, retained_ref = load(retained_path, ordinary.RETAINED_COUNTS_SHA, evidence)
    require(retained["status"] == "corrected_qfo_factorial_swiss_counts_verified"
            and retained["publication_ready"] is False and retained["uncertainty_admitted"] is False,
            "Require retained corrected family-count audit")
    validated_values({**retained, "status": "qfo_factorial_swiss_counts_verified"})
    check(retained["reference"])
    families = retained["families"]
    anchor = retained["cells"][0]
    check(anchor["raw_file"])
    evidence.extend([retained["reference"], anchor["raw_file"]])
    _, truth, members = read_raw(Path(anchor["raw_file"]["path"]), families)
    retained_cells = {row["cell"]: row for row in retained["cells"]}
    cells, excluded, unavailable = [], [], []
    for row in snapshot["rows"]:
        if row["status"] == "supplied_native_admission":
            excluded.append(row["cell"])
            continue
        if row["status"] != "supplied_recovered_scientific_admission":
            unavailable.append(row["cell"])
            continue
        admission, _ = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
        execution_ref = admission["execution_report"]
        require(execution_ref in admission["checked_records"], "Execution was not checked by admission")
        execution, observed = load(execution_ref["path"], execution_ref["sha256"], evidence)
        require(observed == execution_ref, "Execution record metadata differs")
        raw_files = [ref for ref in execution["outputs"] if Path(ref["path"]).parent.name == "SwissTrees"
                     and ref["path"].endswith("raw.txt.gz")]
        require(len(raw_files) == 1, "Require one inventoried recovered SwissTrees raw file")
        raw = raw_files[0]
        require(raw in admission["checked_records"], "Raw file was not checked by admission")
        check(raw)
        evidence.append(raw)
        cell = ordinary.family_cell(raw, families, truth, members, retained["reference_orientation"],
                                    admission["assessment"])
        require(math.isclose(cell["aggregate"]["F1"], row["scores"]["SwissTrees"],
                             rel_tol=0, abs_tol=5e-8), "Recovered counts differ from native endpoint")
        historical = retained_cells[row["cell"]]
        changed = [family for family, fresh, old in zip(families, cell["families"], historical["families"])
                   if fresh != old]
        cell.update(cell=row["cell"], index=row["index"], native_job_id=row["native_job_id"],
                    admission=row["admission"], native_endpoint_f1=row["scores"]["SwissTrees"],
                    measurement_status=row["measurement_status"], resources=None, timing_admitted=False,
                    timing_eligible=False, scientific_recovery=row["scientific_recovery"],
                    retained_family_records_identical=not changed, differing_families=changed,
                    retained_aggregate_identical=cell["aggregate"] == historical["aggregate"])
        cells.append(cell)
    for ref in evidence:
        check(ref)
    return dict(schema="recovered_native_qfo_swiss_family_count_audit_v1", snapshot=snapshot_ref,
        retained_counts=retained_ref, status="supplied_recovered_native_swiss_family_counts_verified",
        families=families, reference=retained["reference"], cells=cells, reference_relation_count=len(truth),
        unavailable_cells=unavailable, successful_native_cells_not_recounted=excluded,
        count_conversion=retained["count_conversion"], checked_inputs=evidence, new_bootstrap_draws=0,
        historical_intervals_attached=False, new_accuracy_or_resource_admission=False,
        independent_confirmation=False, publication_ready=False, source=record(__file__), helpers=[
            record(Path(__file__).with_name(name)) for name in (
                "audit_native_qfo_swiss_counts.py", "audit_qfo_swiss_counts.py", "bootstrap_qfo_factorial.py",
                "export_native_qfo_scientific_scores.py", "validate_qfo_native_assessment.py")],
        limitations=[
            "Reads actual recovered raw counts only after independent accuracy admission; missing cells are not imputed.",
            "Successful native cells use their separate retained count audit and are not recounted here.",
            "Failed native timing remains inadmissible; no new measurement or raw resource replay.",
            "Reference truth and members are anchored to retained corrected-reference evidence.",
            "Identical family sufficient statistics do not establish identical pair decisions or partitions.",
            "No contrast, new bootstrap interval or cached interval is attached by this count audit.",
            "Count-derived F1 uses full precision; native decimal endpoint values remain unchanged.",
            "SwissTrees only: no uncertainty for other QfO endpoints or the secondary mean.",
            "Eighteen development-exposed families may not be exchangeable; not independent confirmation.",
            "Initial HMM search remains on; P-off is not a total-HMM control.",
            "R changes clique to resolved-pair semantics, not only group splitting.",
            "No inference/scoring/FAS rerun or isolated-efficiency claim; shared-host effects remain unknown."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--retained-counts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = audit(args.snapshot, args.snapshot_sha256, args.retained_counts)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(dict(cells=len(result["cells"]), family_records=len(result["cells"]) * len(result["families"]))))


if __name__ == "__main__":
    main()
