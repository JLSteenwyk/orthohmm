"""Assemble eight admitted SwissTrees cells with exact shared-reference checks."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import read_raw, verify_native, REFERENCE_SHA
from benchmark_tools.bootstrap_qfo_factorial import CELLS, validated_values
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_factorial_assessment import ADMITTED_SHA, verify_reuse_binding
from benchmark_tools.validate_qfo_native_assessment import validate_records

BASE_COUNTS_SHA = "546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183"


def file_records(value):
    if isinstance(value, dict):
        if set(value) == {"path", "bytes", "sha256"}:
            yield value
        else:
            for child in value.values():
                yield from file_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from file_records(child)


def selected_admission(report, index):
    if report["cell"] != CELLS[index] or report["index"] != index or report["accuracy_admitted"] is not True:
        raise ValueError("Wrong cell or unadmitted assessment")
    if index in (0, 4):
        if (report["status"] != "admitted_reused_assessment" or report["scoring_rerun"] is not False or
                report["reused_admission"]["sha256"] != ADMITTED_SHA):
            raise ValueError("Baseline reuse not established")
        native = report["admitted_stage"]
        verify_reuse_binding(report["stage"], native, index)
        if report["original_participant"] != native["participant"]:
            raise ValueError("Wrong original participant")
    else:
        if report["status"] != "fresh_factorial_assessment_admitted":
            raise ValueError("Fresh assessment not independently admitted")
        native = report
        scheduler = report["scheduler"]
        if scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0":
            raise ValueError("Assessment did not finish successfully")
        if report["assessment"]["participant"] != f"ohmm_qfo_factorial_{CELLS[index]}":
            raise ValueError("Wrong fresh participant")
    return native


def assemble(entries, baseline):
    if [entry["cell"] for entry in entries] != list(CELLS):
        raise ValueError("Require all eight cells in frozen order")
    if baseline["status"] != "raw_swiss_family_counts_verified" or baseline["reference"]["sha256"] != REFERENCE_SHA:
        raise ValueError("Unverified reference count audit")
    families = baseline["families"]
    anchor = baseline["stages"][0]
    _, expected_truth, expected_members = read_raw(Path(anchor["raw_file"]["path"]), families)
    orientation = baseline["reference_orientation"]
    cells = []
    for entry in entries:
        assessment = entry["assessment"]
        if assessment["swiss_reference_families"] != families:
            raise ValueError("Changed family inventory")
        native = validate_records(assessment["native_assessments"], assessment["participant"], set(families))
        counts, truth, members = read_raw(Path(entry["raw_file"]["path"]), families)
        if truth != expected_truth or members != expected_members:
            raise ValueError("Reference pair identities, labels or members differ")
        if any(sum(counts[f].values()) != orientation[f]["forward_relations"] or
               len(members[f]) != orientation[f]["mapped_proteins"] for f in families):
            raise ValueError("Incomplete native reference coverage")
        values, aggregate = verify_native(counts, native)
        cells.append({"cell": entry["cell"], "raw_file": entry["raw_file"], "aggregate": aggregate,
                      "families": [{"family": f, "counts_without_prior": dict(counts[f]),
                          "statistics_with_prior": values[f], "represented_genes": sorted(members[f])} for f in families]})
    report = {"status": "qfo_factorial_swiss_counts_verified", "reference": baseline["reference"],
              "families": families, "cells": cells, "reference_relation_count": len(expected_truth),
              "shared_represented_genes": baseline["shared_represented_genes"], "publication_ready": False,
              "reference_orientation": orientation,
              "count_conversion": "native confusion count = raw one-direction relation count / 2 + 1",
              "limitations": ["Requires independently admitted assessments; does not rerun inference or scoring.",
                  "Exact reference comparison is against the pinned audited native raw reference universe.",
                  "SwissTrees only; no paired uncertainty for other endpoints or the secondary mean.",
                  "Development-exposed families may be dependent; this audit alone establishes no accuracy advantage."]}
    validated_values(report)
    return report


def audit(manifest_path, manifest_sha, baseline_path):
    manifest_record, baseline_record = record(manifest_path), record(baseline_path)
    if manifest_record["sha256"] != manifest_sha or baseline_record["sha256"] != BASE_COUNTS_SHA:
        raise ValueError("Changed admission inventory or baseline audit")
    manifest, baseline = json.loads(manifest_path.read_text()), json.loads(baseline_path.read_text())
    if [row["cell"] for row in manifest["cells"]] != list(CELLS):
        raise ValueError("Incomplete or reordered admissions")
    checked = {}

    def verify(items):
        for item in items:
            old = checked.get(item["path"])
            if old is not None and old != item:
                raise ValueError("Conflicting identities for one file")
            if old is None:
                check(item)
                checked[item["path"]] = item

    verify([manifest_record, baseline_record, *file_records(baseline)])
    entries = []
    for index, row in enumerate(manifest["cells"]):
        verify([row["admission"]])
        report = json.loads(Path(row["admission"]["path"]).read_text())
        native = selected_admission(report, index)
        verify(file_records(report))
        if index in (0, 4):
            original = json.loads(Path(report["reused_admission"]["path"]).read_text())
            if native != original["records"][1 if index == 0 else 3]:
                raise ValueError("Reused native admission differs from frozen original")
        execution = json.loads(Path(native["execution_report"]["path"]).read_text())
        verify(file_records(execution))
        if execution["status"] != "process_succeeded_pending_independent_admission" or execution["exit_code"] != 0:
            raise ValueError("Execution report contradicts admission")
        paths = [item for item in execution["outputs"] if
                 Path(item["path"]).parent.name == "SwissTrees" and item["path"].endswith("raw.txt.gz")]
        if len(paths) != 1:
            raise ValueError("Ambiguous or missing inventoried raw evidence")
        entries.append({"cell": row["cell"], "assessment": native["assessment"], "raw_file": paths[0]})
    result = assemble(entries, baseline)
    for item in checked.values():
        check(item)
    result.update(admission_inventory=manifest_record, baseline_audit=baseline_record,
                  checked_inputs=list(checked.values()), source=record(__file__), helpers=[record(Path(__file__).with_name(n))
                      for n in ("audit_qfo_swiss_counts.py", "bootstrap_qfo_factorial.py", "run_qfo_factorial_assessment.py",
                                "validate_qfo_native_assessment.py")])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("manifest", "baseline-counts", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.manifest, args.manifest_sha256, args.baseline_counts)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
