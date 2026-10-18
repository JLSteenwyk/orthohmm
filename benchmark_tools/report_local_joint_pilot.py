"""Separate complete-proposal/native/accuracy admission for the local pilot."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.local_rescue_readback import scientific_compatibility
from benchmark_tools.native_pair_turnover import native_pair_turnover
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.report_profile_evidence_pilot import execution_receipt, retained_context, terminal, uncertainty, write_table
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_local_joint_pilot import read_proposals, final_effects
from benchmark_tools.run_membership_hillclimb import PRODUCTS, compare_golden, validate_partition, write_json
from benchmark_tools.run_profile_weighting_pilot import advancement, normalized
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, candidate_control_matches
from orthohmm.accuracy import load_accuracy_checkpoint


LAUNCH = "a5b4f9a3ab716b360659d5e3568a5a7be8dc8c74"
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
PROTOCOL = "hmm_local_joint_postrun_20261018.json"
UNIVERSE = (251378, 12)
LABELS = ["off", "local_joint_v1"]


def launch(repo):
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + LEDGER], cwd=repo)
    entry = json.loads(blob)["local_joint_final_pilot"]
    pin = {key: entry["execution_plan"][key] for key in ("path", "bytes", "sha256")}
    pin["path"] = str((repo / pin["path"]).resolve())
    return entry, pin, {"revision": LAUNCH, "path": LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}


def require_completed(report, plan, job_id):
    if (report.get("status") != "finished_local_readback_and_native_pilot_verified"
            or report.get("job_id") != str(job_id) or report.get("source_revision") != plan["source_revision"]
            or plan.get("cpu") != 16 or plan.get("status") != "prepared_local_joint_readback_and_native_unrun"
            or [row.get("label") for row in report.get("arms", [])] != LABELS
            or any(row.get("status") != "independent_content_readback_verified" or row.get("dataset") != "OrthoBench"
                   or "error" in row for row in report["arms"])
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require the complete successful ordered local pilot, not partial results")


def content(pin):
    return {key: pin[key] for key in ("bytes", "sha256")}


def cache_signature(cache):
    if cache["status"] != "authenticated_composite_cache_built":
        raise ValueError("Composite cache is unfinished")
    root = Path(cache["path"]).resolve()
    result = {}
    for copied in cache["files"]:
        check(copied["source"])
        check(copied["copy"])
        relative = str(Path(copied["copy"]["path"]).resolve().relative_to(root))
        if relative in result or content(copied["source"]) != content(copied["copy"]):
            raise ValueError("Duplicate or altered composite cache copy")
        result[relative] = {"source": copied["source"], **content(copied["copy"])}
    return result


def recheck_proposals(plan, report, output):
    original_folder = Path(plan["output"]) / "readback"
    if report["proposal_readback"] != record(original_folder / "summary.json"):
        raise ValueError("Original full proposal readback identity differs")
    original = json.loads((original_folder / "summary.json").read_text())
    expected_status = "all_original_proposals_and_effective_partition_read_back"
    if original["status"] != expected_status or original["accuracy_evaluated"] is not False:
        raise ValueError("Original proposal verification is incomplete or misclassified")
    read_proposals({**plan, "cpu": 1}, output)
    # Compare persisted structures: species-class keys are integers in memory.
    repeated = json.loads((output / "summary.json").read_text())
    if (repeated["status"] != expected_status or repeated["contexts"] != original["contexts"]
            or repeated["ownership"] != original["ownership"] or repeated["nomination"] != original["nomination"]
            or repeated["effective_partition"] != original["effective_partition"]
            or any(content(repeated[key]) != content(original[key]) for key in ("candidate", "constraints"))
            or original["composite_cache"]["accepted_contexts"] != repeated["composite_cache"]["accepted_contexts"]
            or Path(original["composite_cache"]["path"]).resolve() != (original_folder / "composite_cache").resolve()
            or cache_signature(original["composite_cache"]) != cache_signature(repeated["composite_cache"])):
        raise ValueError("Whole-panel proposal/candidate/cache verification does not reproduce")
    checked = [report["proposal_readback"], record(output / "summary.json")]
    for context in repeated["contexts"]:
        checked.extend(context["checked"])
        if context["tree_inferred"]:
            family = context["context_id"]
            relatives = [Path("gene_trees") / (family + suffix) for suffix in (".raw.nwk", ".rooted.nwk", ".reconciled.nwk")]
            relatives.append(Path("checkpoints") / (family + ".json"))
            for relative in relatives:
                before = record(original_folder / "contexts" / family / relative)
                after = record(output / "contexts" / family / relative)
                if content(before) != content(after):
                    raise ValueError("Producer cached reconciliation artifacts differ from repeated readback")
                checked.extend((before, after))
    for cache in (original["composite_cache"], repeated["composite_cache"]):
        for copied in cache["files"]: checked.extend(copied.values())
    return original, repeated, list({p["path"]: p for p in checked}.values())


def admit_arm(repo, plan, row, original, repeated, names):
    source, root, case = Path(plan["source_root"]), Path(plan["output"]), plan["case"]
    mode = row["label"]
    if mode not in LABELS: raise ValueError("Unknown local pilot arm")
    folder = root / mode
    working = folder / "candidate/orthohmm_working_res"
    cache = case["baseline_output"] if mode == "off" else original["composite_cache"]["path"]
    _, command = commands(source, plan, {**case, "baseline_output": cache}, folder,
                          Path(case["checkpoint"]), case["checkpoint_manifest"]["sha256"])
    preflight = folder / "preflight.json"
    if json.loads(preflight.read_text()) != {"native": command, "candidate": row["candidate"], "cache_source": cache}:
        raise ValueError("Local final-native preflight differs from checked source/cache/policies")
    checked = [record(preflight), *execution_receipt(row["native_execution"], command, folder, "native")]
    origin = case if mode == "off" else repeated
    for key, filename, origin_key in (("candidate_partition", "orthohmm_edges_clustered.txt", "candidate"),
                                      ("membership_constraints", "phylogeny_candidate_merges.json", "constraints")):
        pin = row["candidate"][key]
        if pin != record(working / filename) or content(pin) != content(origin[origin_key]):
            raise ValueError("Local native input differs from independently checked candidates/constraints")
        checked.append(pin)
    candidate, constraints = working / "orthohmm_edges_clustered.txt", working / "phylogeny_candidate_merges.json"
    read_partition(candidate, set(names))
    load_membership_constraints(constraints, candidate)
    prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    if validate_partition(candidate, prediction) != row["coverage"]:
        raise ValueError("Local native complete coverage differs")
    native = replay_native(row, {"command": command})
    if (native["genes"], native["species"]) != UNIVERSE or native["summary"] != row["native_summary"]:
        raise ValueError("Local native universe/summary differs")
    scores, details, pins = readback_scores(repo, "OrthoBench", row, native, folder, "local_postrun")
    if scores != row["scores"] or details != row["details"] or endpoint_vector("OrthoBench", scores) != row["metric_vector"]:
        raise ValueError("Local aggregate/per-family scoring does not reproduce")
    cell = {"label": mode, "dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": endpoint_vector("OrthoBench", scores), "details": details,
            "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
            "products": native["products"], "source_revision": plan["source_revision"], "candidate": row["candidate"],
            "resources": {"native_execution": row["native_execution"]}}
    if mode == "off":
        matches = compare_golden(Path(native["phylogeny"]), case["golden"])
        candidate_control_matches(candidate, constraints, case["candidate"], case["constraints"])
        if matches != row["golden_matches"] or normalized(row) != case["expected"]:
            raise ValueError("Local off control differs from canonical products/scores")
    else:
        phylo = Path(native["phylogeny"])
        turnover = native_pair_turnover(Path(case["golden"][PRODUCTS[2]]["path"]), phylo / PRODUCTS[2])
        effects = final_effects(case, prediction, phylo / PRODUCTS[2], repeated["contexts"])
        changed = record(phylo / "species_tree.rooted.nwk")["sha256"] != case["golden"]["species_tree.rooted.nwk"]["sha256"]
        accepted = repeated["composite_cache"]["accepted_contexts"]
        if (not accepted or row["accepted_contexts"] != accepted or row["final_effects"] != effects
                or row["final_species_tree_changed"] != changed
                or {k: v for k, v in turnover.items() if k != "wall_s"} != {k: v for k, v in row["native_pair_turnover"].items() if k != "wall_s"}):
            raise ValueError("Local activation/membership/pair/species-tree effects differ from actual final outputs")
        cell.update(accepted_contexts=accepted, native_pair_turnover=turnover, final_effects=effects,
                    final_species_tree_changed=changed)
        checked.extend(turnover["files"])
    checked += [*native["checked"], *pins]
    for pin in checked: check(pin)
    cell["checked"] = list({p["path"]: p for p in checked}.values())
    return cell


def local_uncertainty(cells):
    protocol = json.loads(Path(__file__).with_name(PROTOCOL).read_text())["uncertainty"]
    expected = {"replicates": 100000, "seed": 20261018, "alpha": .05, "families": 70,
                "metrics": ["f_score", "precision", "recall"], "multiplicity_endpoints": 3}
    if ([cell["label"] for cell in cells] != LABELS
            or any(protocol.get(key) != value for key, value in expected.items())):
        raise ValueError("Changed local paired-family uncertainty protocol or contrast")
    return uncertainty(cells, multiplicity_endpoints=3)


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every local post-run readback namespace")
    entry, plan_pin, ledger_pin = launch(repo)
    check(plan_pin)
    plan = json.loads(Path(plan_pin["path"]).read_text())
    source, root = Path(plan["source_root"]), Path(plan["output"])
    producer = record(root / "results.json")
    report = json.loads(Path(producer["path"]).read_text())
    if (report["execution_plan"] != plan_pin or plan["source_revision"] != entry["runner_commit"]
            or Path(plan_pin["path"]).parent != root or Path(plan["repo"]) != repo):
        raise ValueError("Local pilot differs from immutable final launch")
    require_completed(report, plan, entry["inference"]["job_id"])
    scheduler = terminal(entry["inference"]["job_id"])
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan["source_revision"]:
        raise ValueError("Local scientific checkout moved")
    if sources(source) + [record(source / "benchmark_tools/run_local_joint_pilot.slurm")] != plan["sources"] or runtime() != plan["runtime"]:
        raise ValueError("Local source inventory/runtime differs from frozen pilot")
    reader = Path(__file__).resolve().parents[1]
    for relative in ("benchmark_tools/run_local_joint_pilot.py", "benchmark_tools/local_rescue_readback.py"):
        if content(record(reader / relative)) != content(record(source / relative)):
            raise ValueError("Proposal verification helpers changed since the source pilot")
    scientific_compatibility(plan["original_source"], reader)
    checked = [plan_pin, producer, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked: check(pin)
    names, species, *_ = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    if (len(names), len(set(species))) != UNIVERSE or len(set(names)) != len(names):
        raise ValueError("Incorrect complete local OrthoBench checkpoint universe")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=reader, text=True).strip()
    reader_files = sources(reader) + [record(reader / "benchmark_tools" / name) for name in (PROTOCOL, "hmm_profile_evidence_readback_20261018.json")]
    for pin in reader_files:
        relative = Path(pin["path"]).relative_to(reader)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=reader)
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Local post-run reader differs from committed source")
    context, context_pins = retained_context(repo)
    checked += [*reader_files, *context_pins]
    output.mkdir(parents=True)
    result = {"status": "running", "producer": producer, "execution_plan": plan_pin, "launch_ledger": ledger_pin,
        "scheduler": scheduler, "source_revision": plan["source_revision"], "reporter_revision": revision,
        "protocol": record(reader / "benchmark_tools" / PROTOCOL), "reporter_files": reader_files,
        "cells": [], "context": context, "accuracy_result_admitted": False, "defaults_promoted": False, "goal_completed": False,
        "scope": "Separate repeated shared-algorithm proposal/native/score readback, not independent biological confirmation. Shared-host cached costs; functional means are not F1 and assessed populations can differ."}
    write_json(output / "summary.json", result)
    try:
        original, repeated, pins = recheck_proposals(plan, report, output / "proposal_recheck")
        checked += pins
        result["proposal_recheck"] = record(output / "proposal_recheck/summary.json")
        for row in report["arms"]:
            cell = admit_arm(repo, plan, row, original, repeated, names)
            if row["label"] == "local_joint_v1":
                off = result["cells"][0]
                control = {"f1": off["scores"]["OrthoBench: OrthoBench F1"], "precision": off["details"]["precision"], "recall": off["details"]["recall"]}
                actual = {"f1": cell["scores"]["OrthoBench: OrthoBench F1"], "precision": cell["details"]["precision"], "recall": cell["details"]["recall"]}
                cell["advancement"] = advancement(control, actual)
                if cell["advancement"] != row["advancement"]: raise ValueError("Local original point gate differs")
            result["cells"].append(cell)
            checked += cell["checked"]
            write_json(output / "summary.json", result)
        result["uncertainty"] = local_uncertainty(result["cells"])
        write_table(output / "scores.tsv", result["cells"], context)
        result["score_table"] = record(output / "scores.tsv")
        for pin in checked: check(pin)
        if runtime() != plan["runtime"]: raise ValueError("Runtime changed during local post-run admission")
        result["checked"] = list({p["path"]: p for p in checked}.values())
        result["advancing_labels"] = [c["label"] for c in result["cells"] if c.get("advancement", {}).get("eligible", False)]
        result.update(status="local_joint_pilot_independently_read_back", accuracy_result_admitted=True)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
