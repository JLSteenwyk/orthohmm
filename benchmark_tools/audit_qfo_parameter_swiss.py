"""Bind parameter score admissions to freshly reconstructed SwissTrees counts."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_factorial_swiss import BASE_COUNTS_SHA, file_records
from benchmark_tools.audit_qfo_swiss_counts import REFERENCE_SHA, read_raw, verify_native
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import ARMS, validated_values
from benchmark_tools.export_qfo_corrected_factorial import extract
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_cpm_assessment import ARMS as CPM_ARMS, validate_stage as validate_cpm_stage
from benchmark_tools.run_qfo_parameter_assessment import VARIANTS, validate_stage
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_records

CONTROL_SHA = "49b7d837b2ba2928b0676c9974e2ec16db5211f1086c366c7bc11805a2ee751f"
PARAMETER_ADMITTER_SHA = "6574b40cc97624c937eb72e2ed27c0822aff0b6ff69156815b0255055d62d81c"
CPM_ADMITTER_SHA = "907e5abbf7ec5b60781dfd48d95452bb67578662b1aefa44af8be131b7542c66"


def validate_admission(arm, report, conversion):
    if arm not in ARMS:
        raise ValueError("Unknown parameter arm")
    if arm == "control":
        row = extract(report, conversion)
        if row["cell"] != "p1_c1_r1" or row["index"] != 7:
            raise ValueError("Wrong full-pipeline control")
    elif arm in CPM_ARMS:
        index = CPM_ARMS.index(arm)
        if (report["status"] != "cpm_assessment_admitted"
                or report["arm"] != arm or type(report["index"]) is not int or report["index"] != index
                or report["source"]["sha256"] != CPM_ADMITTER_SHA
                or report["accuracy_admitted"] is not True or report["publication_ready"] is not False
                or report["conversion"] != conversion or report["context"] != conversion["context"]):
            raise ValueError("Wrong CPM score admission")
        validate_cpm_stage(conversion, index, report["conversion_scheduler"])
    else:
        index = VARIANTS.index(arm)
        if (report["status"] != "corrected_parameter_assessment_admitted"
                or report["variant"] != arm or type(report["index"]) is not int or report["index"] != index
                or report["source"]["sha256"] != PARAMETER_ADMITTER_SHA
                or report["accuracy_admitted"] is not True or report["publication_ready"] is not False
                or report["conversion"] != conversion):
            raise ValueError("Wrong parameter score admission")
        validate_stage(conversion, index, report["conversion_scheduler"])
    if report["assessment"]["participant"] != conversion["participant"]:
        raise ValueError("Assessment participant differs from conversion")


def raw_from_execution(arm, report, execution):
    if arm not in ARMS:
        raise ValueError("Unknown parameter arm")
    scheduler = report["scheduler"]
    if (tuple(scheduler[k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != ("COMPLETED", "0:0", "bizon", "8")
            or execution["status"] != "process_succeeded_pending_independent_admission"
            or type(execution["exit_code"]) is not int or execution["exit_code"] != 0
            or execution["job_id"] != scheduler["JobIDRaw"]
            or type(execution["index"]) is not int or execution["index"] != report["index"]
            or execution["pairs_manifest"] != report["pairs_manifest"]
            or execution["stage"] != report["conversion"]):
        raise ValueError("Execution contradicts score admission")
    key = "cell" if arm == "control" else "arm" if arm in CPM_ARMS else "variant"
    if execution.get(key) != ("p1_c1_r1" if arm == "control" else arm):
        raise ValueError("Wrong execution arm")
    if arm in CPM_ARMS and execution["context"] != report["context"]:
        raise ValueError("CPM execution context differs from admission")
    raw = [r for r in execution["outputs"] if Path(r["path"]).parent.name == "SwissTrees"
           and r["path"].endswith("raw.txt.gz")]
    if len(raw) != 1:
        raise ValueError("Require exactly one inventoried SwissTrees raw file")
    return raw[0]


def assemble(entries, baseline):
    if [r["arm"] for r in entries] != list(ARMS):
        raise ValueError("Require complete ordered parameter inventory")
    if baseline["status"] != "raw_swiss_family_counts_verified" or baseline["reference"]["sha256"] != REFERENCE_SHA:
        raise ValueError("Require pinned audited reference universe")
    families = baseline["families"]
    _, expected_truth, expected_members = read_raw(Path(baseline["stages"][0]["raw_file"]["path"]), families)
    orientation = baseline["reference_orientation"]
    arms = []
    for entry in entries:
        if entry["status"] == "not_admitted":
            if set(entry) != {"arm", "status", "reason"}:
                raise ValueError("Unavailable arm cannot contain prediction evidence")
            arms.append(dict(entry))
            continue
        if entry["status"] != "admitted":
            raise ValueError("Unknown inventory status")
        assessment = entry["assessment"]
        if assessment["swiss_reference_families"] != families:
            raise ValueError("Changed SwissTrees families")
        native = validate_records(assessment["native_assessments"], assessment["participant"], set(families))
        counts, truth, members = read_raw(Path(entry["raw_file"]["path"]), families)
        if truth != expected_truth or members != expected_members:
            raise ValueError("Reference pair identities, labels or members differ")
        if any(sum(counts[f].values()) != orientation[f]["forward_relations"] or
               len(members[f]) != orientation[f]["mapped_proteins"] for f in families):
            raise ValueError("Incomplete reference coverage")
        values, aggregate = verify_native(counts, native)
        if not math.isclose(assessment["endpoints"]["SwissTrees"]["score"], aggregate["F1"], rel_tol=0, abs_tol=1e-6):
            raise ValueError("Raw counts differ from admitted SwissTrees F1")
        arms.append({"arm": entry["arm"], "status": "counts_verified", "raw_file": entry["raw_file"],
            "aggregate": aggregate, "families": [{"family": f, "counts_without_prior": dict(counts[f]),
                "statistics_with_prior": values[f], "represented_genes": sorted(members[f])} for f in families]})
    result = {"status": "corrected_qfo_parameter_swiss_counts_verified", "reference": baseline["reference"],
        "families": families, "arms": arms, "reference_relation_count": len(expected_truth),
        "reference_orientation": orientation, "shared_represented_genes": baseline["shared_represented_genes"],
        "publication_ready": False, "uncertainty_admitted": False,
        "count_conversion": "native confusion count = raw one-direction relation count / 2 + 1",
        "limitations": ["Reference anchor supplies labels/members only, not historical prediction counts.",
            "Only admitted corrected raw predictions contribute counts; missing arms are not imputed.",
            "SwissTrees only; this count audit establishes neither uncertainty nor superiority."]}
    validated_values(result)
    return result


def audit(inventory_path, inventory_sha, baseline_path):
    inventory = read_frozen(inventory_path, inventory_sha)
    baseline = read_frozen(baseline_path, BASE_COUNTS_SHA)
    if (inventory["status"] != "qfo_parameter_score_admission_inventory"
            or [r["arm"] for r in inventory["arms"]] != list(ARMS)):
        raise ValueError("Wrong parameter admission inventory")
    checked = {}
    def verify(items):
        for item in items:
            prior = checked.get(item["path"])
            if prior is not None and prior != item:
                raise ValueError("Conflicting file identities")
            if prior is None:
                check(item)
                checked[item["path"]] = item
    source = record(__file__)
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith("benchmark_tools.") and getattr(m, "__file__", None)]
    verify([source, *helpers, record(inventory_path), record(baseline_path), baseline["reference"], baseline["stages"][0]["raw_file"]])
    entries = []
    for row in inventory["arms"]:
        if row["status"] == "not_admitted":
            entries.append(dict(row))
            continue
        if row["status"] != "admitted":
            raise ValueError("Unknown admission inventory status")
        admission = row["admission"]
        if row["arm"] == "control" and admission["sha256"] != CONTROL_SHA:
            raise ValueError("Control is not the frozen corrected full-pipeline admission")
        verify([admission])
        report = read_frozen(Path(admission["path"]), admission["sha256"])
        pair_record, execution_record = report["pairs_manifest"], report["execution_report"]
        if execution_record not in report["checked_records"]:
            raise ValueError("Execution was not independently checked")
        verify(file_records(report))
        conversion = read_frozen(Path(pair_record["path"]), pair_record["sha256"])
        validate_admission(row["arm"], report, conversion)
        execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
        raw = raw_from_execution(row["arm"], report, execution)
        verify([raw])
        entries.append({"arm": row["arm"], "status": "admitted", "assessment": report["assessment"], "raw_file": raw})
    result = assemble(entries, baseline)
    for item in checked.values():
        check(item)
    result.update(source=source, helpers=helpers, admission_inventory=record(inventory_path),
                  baseline_audit=record(baseline_path), checked_inputs=list(checked.values()))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inventory", "baseline", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--inventory-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.inventory.resolve(), args.inventory_sha256, args.baseline.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
