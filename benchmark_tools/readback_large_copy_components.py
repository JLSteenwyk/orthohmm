"""Reconstruct component candidates and independently read completed native scores."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import run_large_copy_components as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import EXPECTED_GENES, endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_membership_hillclimb import write_json


def readback(plan_path, output):
    if output.exists():
        raise ValueError("Require an untouched independent readback namespace")
    plan = json.loads(plan_path.read_text())
    source, root, repo = (Path(plan[k]) for k in ("source_root", "output", "repo"))
    if Path(__file__).resolve() != source / "benchmark_tools/readback_large_copy_components.py":
        raise ValueError("Use explicitly bound frozen independent reader source")
    producer_path = root / "results.json"
    producer_record = record(producer_path)
    producer = json.loads(producer_path.read_text())
    if (producer["status"] != "independent_content_readback_verified" or producer["execution_plan"] != record(plan_path)
            or producer["source_revision"] != plan["source_revision"] or producer["dataset"] != plan["dataset"]
            or producer.get("inputs_reverified_unchanged") is not True):
        raise ValueError("Require the exact successful component producer and postflight")
    accounting = subprocess.check_output(["sacct", "-j", producer["job_id"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = runner.terminal_job(accounting, producer["job_id"])
    if scheduler["State"] != "COMPLETED":
        raise ValueError("Successful component producer has no successful terminal scheduler record")
    checked = [record(plan_path), producer_record]
    runner.verify(plan, checked)
    output.mkdir(parents=True)
    rows, diagnostics = runner.construct(plan, output)
    if diagnostics["component_counts"] != producer["component_counts"] or diagnostics["ownership"] != producer["ownership"]:
        raise ValueError("Independent component activation or ownership differs from producer")
    for label in ("control", "components500"):
        previous, rebuilt = producer["structural_controls"][label], rows[label]
        for key in ("seed", "copy_split_trace"):
            check(previous[key])
            if previous[key]["sha256"] != rebuilt[key]["sha256"]:
                raise ValueError("Reconstructed component seed or executed trace differs")
            checked.extend((previous[key], rebuilt[key]))
        for key in ("seed_partition", "candidate_partition", "membership_constraints"):
            check(previous["candidate"][key])
            if previous["candidate"][key]["sha256"] != rebuilt["candidate"][key]["sha256"]:
                raise ValueError("Reconstructed native candidate or constraints differ")
            checked.extend((previous["candidate"][key], rebuilt["candidate"][key]))
    dataset = plan["dataset"]
    reused = runner.native_inputs_unchanged(rows)
    if dataset == "QfO" and not reused:
        if (producer["execution_kind"] != "new_cached_graph_component_inference_own_tree"
                or producer["native_execution"]["returncode"] != 0):
            raise ValueError("Changed QfO candidates need successful native own-tree inference")
        native = replay_native(producer, plan["case"])
        scores, details, files = readback_scores(repo, dataset, producer, native, Path(plan["case"]["output"]),
            "ohmm_hillclimb_r1_" + plan["case"]["scoring_id"])
    else:
        expected = ("exact_component_candidate_reuse_no_native_rerun" if dataset == "QfO"
                    else "exact_inactive_candidate_reuse_no_native_rerun")
        if not reused or producer["execution_kind"] != expected or "native_execution" in producer:
            raise ValueError("Exact native-input reuse is mislabeled or has changed candidates")
        native = replay_native(plan["parent_native_row"], plan["parent_native_case"])
        participant = ("ohmm_hillclimb_r1_" + plan["parent_native_case"]["scoring_id"] if dataset == "QfO" else "unused")
        scores, details, files = readback_scores(repo, dataset, plan["parent_native_row"], native,
            Path(plan["parent_native_case"]["output"]), participant)
        if endpoint_vector(dataset, scores) != plan["matched_control"]["scores"]:
            raise ValueError("Independent reuse scores disagree with matched control")
    runner.require_markers(native["summary"])
    vector = endpoint_vector(dataset, scores)
    coverage = {k: native[k] for k in ("genes", "species", "groups", "pairs")}
    if (native["genes"] != EXPECTED_GENES[dataset] or native["species"] != (78 if dataset == "QfO" else 12)
            or producer["coverage"] != coverage or producer["native_products"] != native["products"]
            or producer["native_summary"] != native["summary"] or producer["scores"] != vector):
        raise ValueError("Independent native coverage, products or endpoint scores disagree")
    checked.extend([*native["checked"], *files])
    runner.verify(plan, checked)
    result = {"status": "terminal_component_candidates_and_raw_scores_independently_verified", "dataset": dataset,
        "method_key": "large_copy_components:components500", "scores": vector, "details": details,
        "coverage": coverage, "products": native["products"], "native_summary": native["summary"],
        "component_counts": diagnostics["component_counts"], "execution_kind": producer["execution_kind"],
        "scheduler": scheduler, "execution_plan": record(plan_path), "producer": producer_record,
        "checked": checked, "resources": {k: producer[k] for k in ("native_execution",) if k in producer},
        "scope": "Independent cached candidate reconstruction and raw-score readback; not repeated native inference or fresh pipeline timing.",
        "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = readback(args.execution_plan.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "scores": result["scores"]}))
