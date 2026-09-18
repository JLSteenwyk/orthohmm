"""Validate the archived 27-run DGX native output panel, not its timings."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_dgx_scientific_metadata import SPEC_SHA, PLAN_SHA, PROJECT
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.validate_scaling_outputs import validate


def validate_review(review, spec):
    if (review["status"] != "retained_panel_review_complete_not_timing_admission"
            or review["review_failures"] != 0 or len(review["runs"]) != 27
            or review["scientific_timings_admitted"] != 0
            or review["controlled_workload_verified"] is not False):
        raise ValueError("Require complete metadata review without timing admission")
    for index, (row, run) in enumerate(zip(review["runs"], spec["runs"])):
        if (row["index"] != index or run["index"] != index
                or row["method"] != run["native_method"] or row["proteomes"] != run["proteomes"]
                or row["repeat"] != run["repeat"]
                or row["status"] != "metadata_and_payload_hashes_verified_not_admitted"):
            raise ValueError("Metadata review run identities differ")


def run(archive, results, inventory_path):
    if inventory_path.exists():
        raise FileExistsError(inventory_path)
    spec_path, plan_path, review_path = [results / name for name in (
        "dgx_scientific_execution_20260917.json", "dgx_scaling_commands_20260917.json",
        "dgx_completed_panel_review_20260918.json")]
    spec, plan = read_pinned(spec_path, SPEC_SHA), read_pinned(plan_path, PLAN_SHA)
    if spec["runs"] != plan["runs"] or len(spec["runs"]) != 27:
        raise ValueError("Frozen run panel differs")
    review_record = record(review_path)
    review = read_pinned(review_path, review_record["sha256"])
    validate_review(review, spec)
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "validate_scaling_outputs.py", "simulation_method_outputs.py", "validate_simulation_outputs.py",
        "report_ygob_validation.py", "score_ygob_groups.py", "gnu_time_companion.py")]
    inventory = [record(p) for p in sorted(archive.rglob("*")) if p.is_file()]
    lookup = {r["path"]: r for r in inventory}
    with inventory_path.open("x") as stream:
        json.dump({"archive_root": str(archive), "files": inventory}, stream, indent=2, sort_keys=True)
        stream.write("\n")
    # Metadata is relocated, not regenerated: require the previously reviewed bytes.
    old_root = Path(review["transferred_inventory"][0]["path"])
    old_root = next(parent for parent in old_root.parents if parent.name == "dgx_completed_evidence_20260918")
    for item in review["transferred_inventory"]:
        local = archive / "scaling_native_v1" / Path(item["path"]).relative_to(old_root)
        if lookup.get(str(local)) != {**item, "path": str(local)}:
            raise ValueError("Archived metadata differs from reviewed evidence")
    reports = []
    for frozen in spec["runs"]:
        index = frozen["index"]
        directory = archive / "scaling_native_v1" / f"run_{index:02d}"
        prepared = json.loads((directory / "preparation.json").read_text())
        measured = json.loads((directory / "measurement/results.json").read_text())
        row = {"index": index, "method": frozen["native_method"], "proteomes": frozen["proteomes"], "repeat": frozen["repeat"]}
        try:
            result = validate(prepared["run"], measured, {PROJECT: archive})
            for item in result["checked_files"]:
                if lookup.get(item["path"]) != item:
                    raise ValueError("Validated output differs from pre-validation archive inventory")
            row.update(status="native_outputs_validated", validation=result)
        except Exception as error:
            row.update(status="native_validation_failed", error_type=type(error).__name__, error=str(error))
        reports.append(row)
    for item in [*inventory, *helpers, review_record]:
        check(item)
    if {str(p) for p in archive.rglob("*") if p.is_file()} != set(lookup):
        raise ValueError("Archived native file membership changed")
    return {"status": "native_panel_validation_complete_not_timing_admission", "source": record(__file__),
        "metadata_review": review_record, "specification": record(spec_path), "plan": record(plan_path), "helpers": helpers,
        "archive_inventory": record(inventory_path), "archive_files": len(inventory),
        "archive_bytes": sum(r["bytes"] for r in inventory), "runs": reports,
        "failures": sum(r["status"] == "native_validation_failed" for r in reports),
        "scientific_timings_admitted": 0, "accuracy_evaluated": False,
        "limitations": ["Native output validity is not prediction accuracy or controlled timing.",
            "Historical host-isolation and sampled RSS limitations remain unchanged.",
            "Archive inventory hashes are acquired after completion; reviewed metadata/input hashes supply earlier bindings."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("archive", "results", "inventory", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = run(args.archive.resolve(), args.results.resolve(), args.inventory.resolve())
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"runs": len(report["runs"]), "failures": report["failures"]}))
