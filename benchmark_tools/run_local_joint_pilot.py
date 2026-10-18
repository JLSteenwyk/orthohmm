"""Authenticate local proposals and run matched ordinary-native OG evaluation."""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.local_rescue_readback import scientific_compatibility, read_context, build_composite_cache
from benchmark_tools.native_pair_turnover import native_pair_turnover
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import replay_native, readback_scores, endpoint_vector
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_local_rescue_probe import context_baseline_pairs, validate_plan as validate_probe_plan
from benchmark_tools.run_membership_hillclimb import PRODUCTS, compare_golden, score_orthobench, scoring_records, validate_partition, write_json
from benchmark_tools.run_profile_weighting_pilot import selected_baseline, normalized, advancement
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage, validate_ownership
from orthohmm.accuracy import load_accuracy_checkpoint
from orthohmm.local_rescue import nominate_local_contexts, write_effective_candidates
from orthohmm.phylogeny import PhylogenyConfig, canonical_species_name
from orthohmm.phylogeny_pipeline import _load_sequence_data
from orthohmm.search.profile_expansion import load_global_sequence_database


LAUNCH = "a80e10d45b3914d418e89ad810411ef4baaedd65"
ORIGINAL = "8c4ddfb7bb7dbff7f56fabed91788df725f33367"
PROTOCOL = "benchmark_tools/hmm_local_joint_readback_20261018.json"
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
OBSERVATION = "benchmark_tools/results/local_joint_probe_observation_20261018.json"
UNIVERSE = (251378, 12)


def committed(repo, path):
    payload = subprocess.check_output(["git", "show", LAUNCH + ":" + path], cwd=repo)
    if path == LEDGER:
        return json.loads(payload), {"revision": LAUNCH, "path": path,
                                     "git_blob_sha256": hashlib.sha256(payload).hexdigest()}
    if hashlib.sha256(payload).hexdigest() != record(repo / path)["sha256"]:
        raise ValueError("Changed committed local prerequisite: " + path)
    return json.loads(payload), record(repo / path)


def terminal(job):
    command = ["sacct", "-j", str(job), "-X", "-P", "-o", "JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList"]
    rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
    expected = {"JobIDRaw": str(job), "State": "COMPLETED", "ExitCode": "0:0", "AllocCPUS": "16", "ReqMem": "128G", "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Original local probe is not uniquely terminal successful")
    return expected


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require a fresh output and frozen scientific reader source")
    ledger, ledger_pin = committed(repo, LEDGER)
    observation, observation_pin = committed(repo, OBSERVATION)
    entry = ledger["local_rescue_recovery"]
    pin = {key: entry["execution_plan"][key] for key in ("path", "bytes", "sha256")}
    pin["path"] = str((repo / pin["path"]).resolve())
    check(pin)
    original = json.loads(Path(pin["path"]).read_text())
    validate_probe_plan(original)
    if original["source_revision"] != ORIGINAL or entry["fix_revision"] != ORIGINAL:
        raise ValueError("Wrong original probe scientific revision")
    old_source = Path(original["source_root"])
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=old_source, text=True).strip() != ORIGINAL:
        raise ValueError("Original scientific checkout moved")
    original_sources = sources(old_source) + [record(old_source / "benchmark_tools/run_local_rescue_probe.slurm")]
    if original_sources != original["sources"] or runtime() != original["runtime"]:
        raise ValueError("Original probe source inventory/runtime changed")
    scheduler = terminal(entry["inference"]["job_id"])
    for key in ("producer", "case_report"):
        evidence_pin = dict(observation[key])
        evidence_pin["path"] = str((repo / evidence_pin["path"]).resolve())
        check(evidence_pin)
    producer_path = Path(original["output"]) / "results.json"
    producer = json.loads(producer_path.read_text())
    case_path = Path(original["output"]) / "OrthoBench/probe.json"
    case_report = json.loads(case_path.read_text())
    if (producer["status"] != "applicability_found_accuracy_evaluation_pending" or producer["job_id"] != scheduler["JobIDRaw"]
            or producer["source_revision"] != ORIGINAL or producer["execution_plan"] != pin
            or len(producer["cases"]) != 1 or producer["cases"][0]["report"] != record(case_path)
            or record(producer_path)["sha256"] != observation["producer"]["sha256"]
            or record(case_path)["sha256"] != observation["case_report"]["sha256"]
            or case_report["status"] != "analyzable_contexts_found" or case_report["accepted_contexts"] <= 0):
        raise ValueError("Require the complete positive original OrthoBench applicability result")
    case = selected_baseline(repo)
    probe_case = original["cases"][0]
    for key in ("fasta_directory", "baseline_output", "checkpoint", "checkpoint_manifest", "fastas", "candidate", "constraints", "golden", "aligner", "tree_builder"):
        if case[key] != probe_case[key]:
            raise ValueError("Final-native canonical control differs from original probe: " + key)
    compatibility = scientific_compatibility(old_source, source)
    inputs = [observation_pin, pin, record(producer_path), record(case_path),
              *original["inputs"], *original["sources"], *original["runtime"]["files"],
              *scoring_records(repo, "OrthoBench"), *compatibility["checked"], record(source / PROTOCOL),
              record(Path(case["baseline_output"]) / "orthohmm_working_res/orthohmm_edges_clustered.txt")]
    for directory in (case_path.parent, Path(case["baseline_output"]) / "orthohmm_phylogeny"):
        inputs.extend(record(p) for p in sorted(directory.rglob("*")) if p.is_file())
    inputs = list({p["path"]: p for p in inputs}.values())
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    source_pins = sources(source) + [record(source / "benchmark_tools/run_local_joint_pilot.slurm")]
    for source_pin in source_pins:
        relative = Path(source_pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=source)
            if hashlib.sha256(blob).hexdigest() != source_pin["sha256"]:
                raise ValueError("Reader source differs from committed bytes")
    for name in ("hmm_viterbi.so", "pair_align.so", "kmer_prefilter.so"):
        if record(source / "orthohmm/search/csrc" / name)["sha256"] != record(old_source / "orthohmm/search/csrc" / name)["sha256"]:
            raise ValueError("Require original unchanged CPU libraries")
    for input_pin in [*inputs, *source_pins]:
        check(input_pin)
    plan = {"status": "prepared_local_joint_readback_and_native_unrun", "source_root": str(source), "source_revision": revision,
        "repo": str(repo), "output": str(output), "cpu": 16, "case": case,
        "probe_plan": pin, "producer": record(producer_path), "case_report": record(case_path),
        "launch_ledger": ledger_pin,
        "original_source": str(old_source), "original_scheduler": scheduler, "compatibility": compatibility,
        "inputs": inputs, "sources": source_pins, "runtime": original["runtime"], "tool_path": original["tool_path"],
        "scope": "Repeat original local nominations/profiles/cached reconciliation, then matched ordinary-native OG pilot; shared-host cached incremental costs, not isolated fresh inference timing.",
        "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def read_proposals(plan, folder):
    case = plan["case"]
    probe = json.loads(Path(plan["case_report"]["path"]).read_text())
    names, species, q, t, scores = load_accuracy_checkpoint(Path(case["checkpoint"]), verify=True)
    files = [Path(pin["path"]).name for pin in case["fastas"]]
    sequences, gene_species = _load_sequence_data(case["fasta_directory"], files)
    by_species = {canonical_species_name(f): f for f in files}
    owners = {gene: by_species[taxon] for gene, taxon in gene_species.items()}
    ownership = validate_ownership(names, species, owners)
    if (len(names), len(set(owners.values()))) != UNIVERSE:
        raise ValueError("Readback requires the complete original OrthoBench universe")
    groups = read_partition(Path(case["candidate"]["path"]), set(names))
    indices = {gene: i for i, gene in enumerate(names)}
    nomination = nominate_local_contexts([sorted(indices[g] for g in group) for group in groups], names, species, q, t, scores)
    check(probe["nomination"])
    if nomination != json.loads(Path(probe["nomination"]["path"]).read_text()):
        raise ValueError("Complete original nomination/exclusions/order do not reproduce")
    contexts = nomination["contexts"]
    if [row["context_id"] for row in probe["contexts"]] != [c["context_id"] for c in contexts]:
        raise ValueError("Original context panel is incomplete or reordered")
    constraints = load_membership_constraints(Path(case["constraints"]["path"]), Path(case["candidate"]["path"]))
    baseline = context_baseline_pairs(Path(case["golden"][PRODUCTS[2]]["path"]), contexts, gene_species)
    database = load_global_sequence_database(case["fasta_directory"], files, names)
    config = PhylogenyConfig(mode="reconcile", species_tree_mode="infer", species_tree_rooting="min_variance",
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy", aligner=case["aligner"], tree_builder=case["tree_builder"])
    result = {"status": "running", "contexts": [], "ownership": ownership, "nomination": probe["nomination"], "accuracy_evaluated": False}
    folder.mkdir(parents=True, exist_ok=False)
    proposals = []
    for context, original in zip(contexts, probe["contexts"]):
        check(original["proposal"])
        expected_path = Path(plan["case_report"]["path"]).parent / "contexts" / context["context_id"] / "proposal.json"
        if original["proposal"] != record(expected_path):
            raise ValueError("Original proposal is outside its fixed context namespace")
        actual = read_context(expected_path, context, names, database, gene_species, sequences,
            case["golden"]["species_tree.rooted.nwk"], files, config, constraints, baseline[context["context_id"]],
            folder / "contexts" / context["context_id"], cpu=plan["cpu"])
        if any(actual[key] != original[key] for key in ("accepted", "tree_inferred", "profile_accepted")):
            raise ValueError("Original context activity/admission differs from readback")
        result["contexts"].append(actual)
        proposals.append(json.loads(expected_path.read_text()))
        write_json(folder / "summary.json", result)
    working = folder / "candidate/orthohmm_working_res"
    candidate = working / "orthohmm_edges_clustered.txt"
    result["effective_partition"] = write_effective_candidates(case["candidate"]["path"], names, proposals, candidate, constraints)
    if (record(candidate)["sha256"] != probe["candidate"]["sha256"]
            or result["effective_partition"] != probe["effective_partition"]
            or sum(c["accepted"] for c in result["contexts"]) != probe["accepted_contexts"]
            or sum(c["tree_inferred"] for c in result["contexts"]) != probe["tree_contexts"]):
        raise ValueError("Rebuilt complete partition or context counts differ from original")
    shutil.copy2(case["constraints"]["path"], working / "phylogeny_candidate_merges.json")
    load_membership_constraints(working / "phylogeny_candidate_merges.json", candidate)
    result["candidate"] = record(candidate)
    result["constraints"] = record(working / "phylogeny_candidate_merges.json")
    result["composite_cache"] = build_composite_cache(Path(case["baseline_output"]) / "orthohmm_phylogeny",
        result["contexts"], folder / "composite_cache")
    result["status"] = "all_original_proposals_and_effective_partition_read_back"
    write_json(folder / "summary.json", result)
    return result


def final_effects(case, prediction, confidence, contexts):
    old_path = Path(case["baseline_output"]) / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    old_genes = old_path.read_text().split()
    old_groups = [frozenset(g) for g in read_partition(old_path, set(old_genes))]
    new_groups = [frozenset(g) for g in read_partition(prediction, set(old_genes))]
    old_owner = {g: group for group in old_groups for g in group}
    new_owner = {g: group for group in new_groups for g in group}
    accepted = [row for row in contexts if row["accepted"]]
    context_by_gene, pairs = {}, {row["context_id"]: {} for row in accepted}
    for row in accepted:
        for gene in row["genes"]:
            if gene in context_by_gene: raise ValueError("Final local contexts overlap")
            context_by_gene[gene] = row["context_id"]
    with Path(confidence).open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            a, b = row["gene_a"], row["gene_b"]
            context = context_by_gene.get(a)
            if context is not None and context_by_gene.get(b) == context:
                pair = tuple(sorted((a, b)))
                if pair in pairs[context]: raise ValueError("Duplicate final context pair")
                pairs[context][pair] = row["confidence"]
    details = []
    for row in accepted:
        final = pairs[row["context_id"]]
        baseline = {tuple(sorted(p)) for p in row["baseline_high_pairs"]}
        novel = {tuple(sorted(p)) for p in row["joint_admission"]["novel_cross_high_pairs"]}
        details.append({"context_id": row["context_id"], "baseline_high_pairs": len(baseline),
            "baseline_high_pairs_absent_final": len(baseline - final.keys()),
            "baseline_high_pairs_downgraded_final": sum(final[p] != "high" for p in baseline & final.keys()),
            "provisional_novel_cross_high_pairs": len(novel),
            "provisional_novel_pairs_present_final": len(novel & final.keys()),
            "provisional_novel_pairs_high_final": sum(final[p] == "high" for p in novel & final.keys())})
    return {"baseline_root_ogs": len(old_groups), "final_root_ogs": len(new_groups),
            "unchanged_root_ogs": len(set(old_groups) & set(new_groups)),
            "genes_with_changed_root_og_membership": sum(old_owner[g] != new_owner[g] for g in old_owner),
            "contexts": details, "scope": "Membership and provisional-versus-final native confidence turnover, not reference accuracy. Provisional pair preservation does not guarantee final preservation."}


def run(path):
    plan = json.loads(path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if (Path(__file__).resolve() != source / "benchmark_tools/run_local_joint_pilot.py" or plan["cpu"] != 16
            or (output / "results.json").exists() or plan["status"] != "prepared_local_joint_readback_and_native_unrun"):
        raise ValueError("Require the fixed frozen local runner and a fresh attempt")
    pins = [record(path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins: check(pin)
    _, ledger = committed(repo, LEDGER)
    if ledger != plan["launch_ledger"]:
        raise ValueError("Immutable original launch binding differs")
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed scientific runtime or CUDA visibility")
    scientific_compatibility(plan["original_source"], source)
    report = {"status": "running", "execution_plan": record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
        "source_revision": plan["source_revision"], "arms": [], "scope": plan["scope"], "host_before": host_activity(),
        "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    try:
        started = time.monotonic()
        readback = read_proposals(plan, output / "readback")
        report["proposal_readback"] = record(output / "readback/summary.json")
        report["readback_wall_s"] = time.monotonic() - started
        case = plan["case"]
        env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
               "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
        control = None
        for mode in ("off", "local_joint_v1"):
            row = {"dataset": "OrthoBench", "label": mode, "status": "running", "scores": {}}
            report["arms"].append(row)
            folder = output / mode
            working = folder / "candidate/orthohmm_working_res"
            working.mkdir(parents=True, exist_ok=False)
            origin = case if mode == "off" else readback
            shutil.copy2(origin["candidate"]["path"], working / "orthohmm_edges_clustered.txt")
            shutil.copy2(origin["constraints"]["path"], working / "phylogeny_candidate_merges.json")
            row["candidate"] = {"candidate_partition": record(working / "orthohmm_edges_clustered.txt"),
                                "membership_constraints": record(working / "phylogeny_candidate_merges.json")}
            if (row["candidate"]["candidate_partition"]["sha256"] != origin["candidate"]["sha256"]
                    or row["candidate"]["membership_constraints"]["sha256"] != origin["constraints"]["sha256"]):
                raise ValueError("Final native candidate/constraint copy differs from checked source")
            cache = case["baseline_output"] if mode == "off" else readback["composite_cache"]["path"]
            native_case = {**case, "baseline_output": cache}
            _, command = commands(source, plan, native_case, folder, Path(case["checkpoint"]), case["checkpoint_manifest"]["sha256"])
            write_json(folder / "preflight.json", {"native": command, "candidate": row["candidate"], "cache_source": cache})
            if mode != "off":
                for copied in readback["composite_cache"]["files"]:
                    check(copied["copy"])
                    check(copied["source"])
            row["native_execution"] = timed_stage(command, folder, "native", source, env)
            if row["native_execution"]["returncode"]: raise RuntimeError("Final ordinary native inference failed")
            phylo = folder / "native/orthohmm_phylogeny"
            prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["outputs"] = [record(phylo / name) for name in PRODUCTS]
            row["coverage"] = validate_partition(working / "orthohmm_edges_clustered.txt", prediction)
            native = replay_native(row, {"command": command})
            if (native["genes"], native["species"]) != UNIVERSE: raise ValueError("Final native universe differs")
            row["native_summary"] = native["summary"]
            row.update(score_orthobench(repo, prediction))
            scores, details, score_pins = readback_scores(repo, "OrthoBench", row, native, folder, "local_joint")
            if scores != row["scores"]: raise ValueError("Independent local pilot scores differ")
            row.update(details=details, scoring_readback=score_pins, metric_vector=endpoint_vector("OrthoBench", scores))
            for pin in [*pins, *native["checked"], *score_pins]: check(pin)
            if mode == "off":
                row["golden_matches"] = compare_golden(phylo, case["golden"])
                if normalized(row) != case["expected"]: raise ValueError("Final off control scores differ from canonical")
                control = normalized(row)
            else:
                row["advancement"] = advancement(control, normalized(row))
                row["accepted_contexts"] = readback["composite_cache"]["accepted_contexts"]
                if not row["accepted_contexts"]: raise ValueError("Positive original local admissions disappeared")
                row["native_pair_turnover"] = native_pair_turnover(Path(case["golden"][PRODUCTS[2]]["path"]), phylo / PRODUCTS[2])
                row["final_effects"] = final_effects(case, prediction, phylo / PRODUCTS[2], readback["contexts"])
                row["final_species_tree_changed"] = record(phylo / "species_tree.rooted.nwk")["sha256"] != case["golden"]["species_tree.rooted.nwk"]["sha256"]
            row["status"] = "independent_content_readback_verified"
            if mode != "off":
                for copied in readback["composite_cache"]["files"]:
                    check(copied["copy"])
                    check(copied["source"])
            write_json(output / "results.json", report)
        if runtime() != plan["runtime"]: raise ValueError("Scientific runtime changed during final evaluation")
        report["status"] = "finished_local_readback_and_native_pilot_verified"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        if report["arms"] and report["arms"][-1]["status"] == "running":
            report["arms"][-1].update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["host_after"] = host_activity()
        write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for key in ("repo", "source", "output"): prep.add_argument("--" + key, required=True, type=Path)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])


if __name__ == "__main__":
    main()
