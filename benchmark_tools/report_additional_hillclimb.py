"""Incrementally admit the separate QfO continuation and shared-hit trials."""

import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import verified_native
from benchmark_tools.report_current_hillclimb import canonical_cell, upsert_cell
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, screen, snapshot
from benchmark_tools.run_full_search_hillclimb import runtime, validate_native, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import DATASETS, write_json
from benchmark_tools.run_shared_hit_pipeline import selected_control


def check_plan(repo, root, report):
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    if Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root:
        raise ValueError("Producer plan belongs to another repository or output")
    pins = [report["execution_plan"], *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    if runtime() != plan["runtime"]:
        raise ValueError("Scientific runtime differs from the pinned producer")
    return plan, pins


def fresh_cell(repo, dataset, row, native, folder, participant, revision):
    scores, details, pins = readback_scores(repo, dataset, row, native, folder, participant)
    checked = [*native["checked"], *pins]
    for pin in checked:
        check(pin)
    return {"dataset": dataset, "status": "independent_content_readback_verified",
            "scores": endpoint_vector(dataset, scores), "details": details,
            "products": native["products"], "checked": checked,
            "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
            "source_revision": revision, "resources": row["execution"],
            "execution_kind": "fresh_full_pipeline_with_separate_scoring",
            "scope": "Fresh variant inference. Shared-host cost is not comparable to cached control replay cost."}


def continuation_cell(repo, plan, report):
    case = plan["case"]["fresh_case"]
    if (case["dataset"] != "QfO" or case["label"] != "cap400"
            or case["search_max_candidates_per_query"] != 400 or report["execution"]["returncode"] != 0):
        raise ValueError("Require the successful original cap400 continuation")
    verify_input_inventory(case)
    native = verified_native(case)
    cell = fresh_cell(repo, "QfO", report, native, Path(plan["case"]["scoring_output"]),
                      "ohmm_hillclimb_r1_" + plan["case"]["scoring_id"], plan["scientific_revision"])
    for pin in report["control"]["checked"]:
        check(pin)
    cell["execution_kind"] = "fresh_original_arrival_cap400_continuation"
    return cell


def shared_cell(repo, root, plan, report, row):
    dataset = row["dataset"]
    cases = [case for case in plan["cases"] if case["dataset"] == dataset and case["label"] == "shared_hit1"]
    if len(cases) != 1 or report["source_revision"] != plan["source_revision"]:
        raise ValueError("Ambiguous shared-hit case or changed source revision")
    case = cases[0]
    if (case["search_min_total_hits"] != 1 or case["search_max_candidates_per_query"] != 100
            or row["execution"]["returncode"] != 0 or row["control"]["initial_min_total_hits"] != 4):
        raise ValueError("Wrong shared-hit intervention or control threshold")
    verify_input_inventory(case)
    native = validate_native(case, 16)
    metadata = native["metrics"]["metadata"]
    if metadata.get("search_result_order") != "canonical" or metadata.get("search_min_total_hits") != 1:
        raise ValueError("Actual native shared-hit mechanism differs from the plan")
    harness = native["metrics"]["harness"]
    if harness.get("git_commit") != plan["source_revision"] or harness.get("exit_code") != 0:
        raise ValueError("Shared-hit native harness source or exit differs")
    source = Path(plan["source_root"])
    expected_sources = [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    expected_sources.append(record(source / "benchmark_tools/benchmark_production.py"))
    actual_sources = {r["path"]: (r["bytes"], r["sha256"]) for r in harness["source_manifest"]}
    expected = {str(Path(r["path"]).relative_to(source)): (r["bytes"], r["sha256"]) for r in expected_sources}
    expected_inputs = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in case["inputs"]}
    if (len(actual_sources) != len(harness["source_manifest"]) or actual_sources != expected
            or len(harness["input_manifest"]) != len(case["inputs"])
            or {r["path"]: (r["bytes"], r["sha256"]) for r in harness["input_manifest"]} != expected_inputs):
        raise ValueError("Shared-hit native source/input manifests differ")
    control_path = report["control_report"]
    check(control_path)
    control_report = json.loads(Path(control_path["path"]).read_text())
    control_root = Path(plan["control_root"])
    control_pin = control_report["execution_plan"]
    if (Path(control_path["path"]).resolve() != control_root / "results.json"
            or Path(control_pin["path"]).resolve() != control_root / "execution_plan.json"
            or control_pin not in plan["inputs"]):
        raise ValueError("Canonical control is not bound to the frozen trial plan")
    check(control_pin)
    control_plan = json.loads(Path(control_pin["path"]).read_text())
    if not (control_plan["source_revision"] == control_report["source_revision"] == row["control"]["source_revision"]):
        raise ValueError("Canonical control source revision differs")
    control = selected_control(control_report, dataset)
    control_cell = canonical_cell(repo, control_root, control, control_report)
    if (endpoint_vector(dataset, row["control"]["scores"]) != control_cell["scores"]
            or (native["genes"], native["species"]) !=
            (control_cell["coverage"]["genes"], control_cell["coverage"]["species"])):
        raise ValueError("Independent canonical control scores or universe differ")
    folder = Path(case["output"])
    native.update(prediction=str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                  phylogeny=str(folder / "output/orthohmm_phylogeny"),
                  checked=[*native["products"].values(), record(folder / "metrics.json"),
                           record(folder / "output/orthohmm_phylogeny/provenance_manifest.json"),
                           record(folder / "output/orthohmm_phylogeny/reconciliation_summary.json"),
                           *[record(p) for p in sorted((folder / "output/orthohmm_working_res/high_sensitivity_checkpoint").iterdir())]])
    nonce = hashlib.sha256(str(root).encode()).hexdigest()[:16]
    cell = fresh_cell(repo, dataset, row, native, folder, "ohmm_hillclimb_r1_sharedhit_" + nonce,
                      plan["source_revision"])
    cell["matched_control"] = {"report": control_path, "cell": control_cell}
    cell["execution_kind"] = "fresh_canonical_initial_shared_hit1_cap100"
    return cell


def build(repo, parent_path, output, continuation_root=None, shared_root=None):
    if output.exists() or continuation_root is None and shared_root is None:
        raise ValueError("Require a fresh output and at least one producer namespace")
    parent = json.loads(parent_path.read_text())
    check(parent["score_table"])
    methods = deepcopy(parent["methods"])
    output.mkdir()
    captures = output / "snapshots"
    captures.mkdir()
    observations = []
    for kind, root in (("continuation", continuation_root), ("shared_hit", shared_root)):
        if root is None:
            continue
        path = root / "results.json"
        if not path.exists():
            observations.append({"kind": kind, "path": str(path), "status": "report_missing_execution_state_not_inferred"})
            continue
        report, capture = snapshot(path, captures, kind + ".json")
        plan, pins = check_plan(repo, root, report)
        observations.append({"kind": kind, "producer_status": report["status"], "snapshot": capture["snapshot"]})
        rows = [("QfO", report)] if kind == "continuation" else [(row["dataset"], row) for row in report["arms"]]
        if any(dataset not in DATASETS for dataset, _ in rows) or len({dataset for dataset, _ in rows}) != len(rows):
            raise ValueError("Wrong or duplicate producer dataset inventory")
        for dataset, row in rows:
            key = "full_search:cap400" if kind == "continuation" else "shared_hit:shared_hit1"
            cell = {"dataset": dataset, "status": row["status"], "producer_status": row["status"],
                    "scores": {}, "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            cell.update({k: row[k] for k in ("error", "error_type") if k in row})
            success = ("original_cap400_native_and_all_six_scores_verified" if kind == "continuation"
                       else "independent_content_and_score_readback_verified")
            if row["status"] == success:
                cell.update(continuation_cell(repo, plan, row) if kind == "continuation"
                            else shared_cell(repo, root, plan, report, row))
                if kind == "shared_hit":
                    control = {**cell["matched_control"]["cell"], "dataset": dataset}
                    upsert_cell(methods, "canonical_collection:canonical", "OrthoHMM canonical collection (experimental)", control)
            upsert_cell(methods, key, key, cell)
        for pin in pins:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Scientific runtime changed during admission")
    allowed = {"full_search:cap400", "shared_hit:shared_hit1", "canonical_collection:canonical"}
    for previous in parent["methods"]:
        if previous["key"] not in allowed and previous != next(m for m in methods if m["key"] == previous["key"]):
            raise ValueError("Unrelated method changed")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", *ENDPOINTS, "complete_eight_endpoints"])
        for method in methods:
            writer.writerow([method["label"], *(method["scores"].get(k, "NA") for k in ENDPOINTS),
                             set(method["scores"]) == set(ENDPOINTS)])
    result = deepcopy(parent)
    # Earlier reporting-code identities remain in the pinned parent, not this run.
    result.pop("reporter_revision", None)
    result["reporter_files"] = [record(Path(__file__).parent / name) for name in (
        "report_additional_hillclimb.py", "report_current_hillclimb.py", "report_hillclimb.py",
        "recover_full_qfo_search.py", "run_full_search_hillclimb.py", "run_membership_hillclimb.py",
        "run_shared_hit_pipeline.py")]
    result.update(status="incremental_separate_trial_readback_no_default_promotion", parent_readback=record(parent_path),
                  reporter_source=record(__file__), methods=methods, score_table=record(output / "scores.tsv"),
                  additional_trial_observations=observations, defaults_promoted=False, goal_completed=False)
    vectors = {m["key"]: m["scores"] for m in methods
               if m["key"] in {"shared_hit:shared_hit1", "canonical_collection:canonical"}}
    result["round_screening"]["shared_hit"] = screen(vectors, "canonical_collection:canonical")
    full = {m["key"].split(":", 1)[1]: m["scores"] for m in methods if m["key"].startswith("full_search:")}
    result["round_screening"]["full_search"] = screen(full, "control100")
    result["limitations"].extend(["Newly completed trial cells are independently read back; unrelated earlier admissions are inherited, not rescored.",
        "Pending, failed, absent and unverified cells do not supply scores. Cell admission is not a claim that the producer job is terminal.",
        "Fresh shared-hit inference and cached canonical controls have different timing scopes. FAS remains unseeded."])
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "parent", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("continuation-root", "shared-root"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    result = build(args.repo.resolve(), args.parent.resolve(), args.output.resolve(),
                   args.continuation_root.resolve() if args.continuation_root else None,
                   args.shared_root.resolve() if args.shared_root else None)
    print(result["status"], len(result["methods"]))
