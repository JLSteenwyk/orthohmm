"""Read native SwissTrees family counts without inheriting cached intervals."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import read_raw, verify_native
from benchmark_tools.bootstrap_qfo_factorial import validated_values
from benchmark_tools.export_native_qfo_factorial_scores import collect
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.validate_qfo_native_assessment import validate_records


RETAINED_COUNTS_SHA = "c9d8bf02ef6f287c56d073fa61f165c982caf834a94fc6e44e6589f03e46eba6"


def family_cell(raw, families, truth, members, orientation, assessment):
    require(assessment["swiss_reference_families"] == families, "Changed family inventory")
    native = validate_records(assessment["native_assessments"], assessment["participant"], set(families))
    counts, observed_truth, observed_members = read_raw(Path(raw["path"]), families)
    require(observed_truth == truth and observed_members == members,
            "Reference pair labels or represented members differ")
    require(all(sum(counts[f].values()) == orientation[f]["forward_relations"]
                and len(observed_members[f]) == orientation[f]["mapped_proteins"] for f in families),
            "Incomplete raw reference universe")
    values, aggregate = verify_native(counts, native)
    return {"raw_file": raw, "aggregate": aggregate, "families": [
        {"family": f, "counts_without_prior": dict(counts[f]),
         "statistics_with_prior": values[f], "represented_genes": sorted(observed_members[f])}
        for f in families]}


def audit(snapshot_path, snapshot_sha, retained_path):
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    retained, retained_ref = load(retained_path, RETAINED_COUNTS_SHA, evidence)
    require(snapshot["schema"] == "native_qfo_reporting_snapshot_v1"
            and snapshot["source"] == record(collect.__code__.co_filename), "Wrong native snapshot/source")
    admissions = [(r["admission"]["path"], r["admission"]["sha256"])
                  for r in snapshot["rows"] if r["status"] == "supplied_native_admission"]
    require(admissions, "No supplied native admission")
    replay = collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], admissions)
    require(all(snapshot[k] == value for k, value in replay.items()), "Native snapshot replay differs")
    require(retained["status"] == "corrected_qfo_factorial_swiss_counts_verified"
            and retained["publication_ready"] is False and retained["uncertainty_admitted"] is False,
            "Require retained corrected family-count audit")
    validated_values({**retained, "status": "qfo_factorial_swiss_counts_verified"})
    check(retained["reference"])
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"], retained["reference"]])
    families = retained["families"]
    anchor = retained["cells"][0]
    check(anchor["raw_file"])
    evidence.append(anchor["raw_file"])
    _, truth, members = read_raw(Path(anchor["raw_file"]["path"]), families)
    retained_cells = {r["cell"]: r for r in retained["cells"]}
    cells, unavailable = [], []
    for row in snapshot["rows"]:
        if row["status"] != "supplied_native_admission":
            unavailable.append(row["cell"])
            continue
        admission, _ = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
        execution_ref = admission["execution_report"]
        require(execution_ref in admission["checked_records"], "Execution was not checked by admission")
        execution, observed = load(execution_ref["path"], execution_ref["sha256"], evidence)
        require(observed == execution_ref, "Execution record metadata differs")
        raw_files = [r for r in execution["outputs"] if Path(r["path"]).parent.name == "SwissTrees"
                     and r["path"].endswith("raw.txt.gz")]
        require(len(raw_files) == 1, "Require one inventoried SwissTrees raw file")
        raw = raw_files[0]
        require(raw in admission["checked_records"], "Raw file was not checked by admission")
        check(raw)
        evidence.append(raw)
        cell = family_cell(raw, families, truth, members, retained["reference_orientation"], admission["assessment"])
        require(math.isclose(cell["aggregate"]["F1"], row["scores"]["SwissTrees"],
                             rel_tol=0, abs_tol=5e-8), "Reconstructed F1 differs from native endpoint")
        historical = retained_cells[row["cell"]]
        changed = [f for f, fresh, old in zip(families, cell["families"], historical["families"])
                   if fresh != old]
        cell.update(cell=row["cell"], index=row["index"], native_job_id=row["native_job_id"],
                    admission=row["admission"], native_endpoint_f1=row["scores"]["SwissTrees"],
                    retained_family_records_identical=not changed, differing_families=changed,
                    retained_aggregate_identical=cell["aggregate"] == historical["aggregate"])
        cells.append(cell)
    for ref in evidence:
        check(ref)
    return {"schema": "native_qfo_swiss_family_count_audit_v1", "snapshot": snapshot_ref,
            "retained_counts": retained_ref, "status": "supplied_native_swiss_family_counts_verified",
            "families": families, "reference": retained["reference"], "cells": cells,
            "reference_relation_count": len(truth), "unavailable_cells": unavailable,
            "count_conversion": retained["count_conversion"], "checked_inputs": evidence,
            "new_bootstrap_draws": 0, "historical_intervals_attached": False,
            "new_accuracy_or_resource_admission": False, "independent_confirmation": False,
            "publication_ready": False, "source": record(__file__), "helpers": [
                record(Path(__file__).with_name(name)) for name in (
                    "audit_qfo_swiss_counts.py", "bootstrap_qfo_factorial.py",
                    "export_native_qfo_factorial_scores.py", "validate_qfo_native_assessment.py")],
            "limitations": [
                "Reads actual raw counts only for supplied native admissions; missing native cells are not imputed.",
                "Exact reference truth and members are anchored to retained corrected-reference raw evidence.",
                "Identical family sufficient statistics do not prove identical pair decisions or whole partitions.",
                "No contrast, bootstrap interval or cached interval is attached by this count audit.",
                "Count-derived F1 uses full precision; native decimal endpoint values remain unchanged.",
                "SwissTrees only: no uncertainty for other QfO endpoints or the secondary mean.",
                "Eighteen development-exposed families may not be exchangeable; this is not independent confirmation.",
                "Initial HMM search remains on; downstream P-off is not a total-HMM control.",
                "R changes group-clique to resolved-pair semantics, not only group splitting.",
                "No inference/scoring/FAS rerun, raw-resource replay, new job or isolated efficiency claim."]}


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
    print(json.dumps({"cells": len(result["cells"]), "family_records": len(result["cells"]) * len(result["families"])}))


if __name__ == "__main__":
    main()
