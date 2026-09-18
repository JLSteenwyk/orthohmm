"""Bind eight corrected-release admissions to fresh SwissTrees family counts."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_factorial_swiss import BASE_COUNTS_SHA, assemble
from benchmark_tools.export_qfo_corrected_factorial import CELLS, extract
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def execution_raw(report, execution):
    scheduler = report["scheduler"]
    if ((scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"])
            != ("COMPLETED", "0:0", "bizon", "8")
            or execution["status"] != "process_succeeded_pending_independent_admission"
            or type(execution["exit_code"]) is not int or execution["exit_code"] != 0
            or execution["job_id"] != scheduler["JobIDRaw"]
            or execution["index"] != report["index"] or execution["cell"] != report["cell"]
            or execution["pairs_manifest"] != report["pairs_manifest"]
            or execution["stage"] != report["conversion"]):
        raise ValueError("Corrected execution contradicts admission")
    paths = [r for r in execution["outputs"] if Path(r["path"]).parent.name == "SwissTrees"
             and r["path"].endswith("raw.txt.gz")]
    if len(paths) != 1:
        raise ValueError("Require one inventoried SwissTrees raw file")
    return paths[0]


def audit(inventory_path, inventory_sha, baseline_path):
    inventory = read_frozen(inventory_path, inventory_sha)
    baseline = read_frozen(baseline_path, BASE_COUNTS_SHA)
    if [r["cell"] for r in inventory["cells"]] != list(CELLS):
        raise ValueError("Require all eight corrected admissions in frozen order")
    checked = [record(inventory_path), record(baseline_path), baseline["reference"],
               baseline["stages"][0]["raw_file"]]
    for item in checked:
        check(item)
    entries, rows = [], []
    for index, entry in enumerate(inventory["cells"]):
        admission = entry["admission"]
        check(admission)
        report = read_frozen(Path(admission["path"]), admission["sha256"])
        pair_record, execution_record = report["pairs_manifest"], report["execution_report"]
        if execution_record not in report["checked_records"]:
            raise ValueError("Execution record was not independently checked")
        check(pair_record)
        check(execution_record)
        conversion = read_frozen(Path(pair_record["path"]), pair_record["sha256"])
        row = extract(report, conversion)
        if row["index"] != index:
            raise ValueError("Admission index differs from inventory")
        execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
        raw = execution_raw(report, execution)
        check(raw)
        entries.append({"cell": row["cell"], "assessment": report["assessment"], "raw_file": raw})
        rows.append(row)
        checked.extend([admission, pair_record, execution_record, raw])
    result = assemble(entries, baseline)
    for row, cell in zip(rows, result["cells"]):
        if not math.isclose(row["scores"]["SwissTrees"], cell["aggregate"]["F1"], rel_tol=0, abs_tol=1e-6):
            raise ValueError("Reconstructed SwissTrees F1 differs from admitted endpoint")
    for item in checked:
        check(item)
    result.update(status="corrected_qfo_factorial_swiss_counts_verified", uncertainty_admitted=False,
        admission_inventory=record(inventory_path), baseline_audit=record(baseline_path),
        checked_inputs=checked, source=record(__file__), helpers=[record(Path(__file__).with_name(n)) for n in (
            "audit_qfo_factorial_swiss.py", "audit_qfo_swiss_counts.py", "bootstrap_qfo_factorial.py",
            "export_qfo_corrected_factorial.py", "run_qfo_corrected_factorial_assessment.py")])
    result["limitations"].append("Historical raw evidence anchors reference labels/members only; every prediction count is newly read from corrected execution output.")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("inventory", "baseline", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.inventory.resolve(), args.inventory_sha256, args.baseline.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
