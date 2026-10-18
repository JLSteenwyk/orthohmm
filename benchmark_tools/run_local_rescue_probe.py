"""Freeze and probe label-free local-rescue applicability before accuracy scoring."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.plan_singleton_tie_followup import verified_cell
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage, validate_ownership
from benchmark_tools.validate_profile_runtime import require_profile_runtime
from orthohmm.accuracy import load_accuracy_checkpoint
from orthohmm.local_rescue import evaluate_local_context, nominate_local_contexts, write_effective_candidates
from orthohmm.phylogeny import PhylogenyConfig, canonical_species_name
from orthohmm.phylogeny_pipeline import _load_sequence_data
from orthohmm.search.profile_expansion import load_global_sequence_database


UNIVERSES = {"OrthoBench": (251378, 12), "QfO": (984137, 78), "ThreeKingdoms": (443217, 12)}
PROTOCOL = "benchmark_tools/hmm_local_joint_rescue_20261018.json"
CANONICAL = "benchmarks/results/canonical_downstream_r1_v1_20261014"
ADMISSION = "benchmark_tools/results/hillclimb_comparison_readback_v27_20261016/summary.json"
PARAMETERS = {"minimum_family_genes": 3, "maximum_joint_genes": 200, "minimum_family_species": 2,
              "minimum_directional_query_genes": 2, "minimum_directional_target_genes": 2,
              "maximum_contexts": 128, "initial_hit_chunk_size": 500000, "profile_evalue": 0.0001,
              "profile_band": 64, "aligned_coverage": 0.8, "minimum_profile_supported_genes_per_side": 2,
              "minimum_profile_supported_species_per_side": 2, "minimum_novel_high_confidence_genes_per_side": 2,
              "minimum_novel_high_confidence_species_per_side": 2, "post_result_parameter_retuning": False}


def selected_cases(repo):
    root = repo / CANONICAL
    plan_path, report_path, admission_path = root / "execution_plan.json", root / "results.json", repo / ADMISSION
    plan, report, admission = [json.loads(path.read_text()) for path in (plan_path, report_path, admission_path)]
    check(report["execution_plan"])
    if report["execution_plan"] != record(plan_path):
        raise ValueError("Canonical report is not bound to its execution plan")
    cases = []
    for dataset, universe in UNIVERSES.items():
        admitted = verified_cell(admission, "canonical_collection:canonical", dataset)
        if (admitted["coverage"]["genes"], admitted["coverage"]["species"]) != universe:
            raise ValueError("Admitted canonical universe differs")
        rows = [row for row in report["arms"] if (row["dataset"], row["label"]) == (dataset, "canonical")]
        configs = [case for case in plan["cases"] if case["dataset"] == dataset]
        if len(rows) != 1 or len(configs) != 1:
            raise ValueError("Require a unique original canonical case")
        row, config = rows[0], configs[0]
        native = root / dataset / "canonical/native"
        golden = {name: record(native / "orthohmm_phylogeny" / name)
                  for name in (*PRODUCTS, "species_tree.rooted.nwk")}
        for name in PRODUCTS:
            original = admitted["products"][name]
            if any(original[key] != golden[name][key] for key in ("bytes", "sha256")):
                raise ValueError("Retained canonical products differ from independent admission")
        if golden["species_tree.rooted.nwk"] not in admitted["checked"]:
            raise ValueError("Baseline species tree differs from its independent admission")
        manifest = row["ordering"]["manifest"]
        check(manifest)
        candidate, constraints = (row["candidate"][key] for key in ("candidate_partition", "membership_constraints"))
        for pin in (candidate, constraints, *config["fastas"]):
            check(pin)
        cases.append({"dataset": dataset, "universe": list(universe), "fasta_directory": config["fasta_directory"],
                      "fastas": config["fastas"], "checkpoint": str(Path(manifest["path"]).parent),
                      "checkpoint_manifest": manifest, "candidate": candidate, "constraints": constraints,
                      "baseline_output": str(native), "golden": golden,
                      "original_producer_status": row["status"],
                      "aligner": config["aligner"], "tree_builder": config["tree_builder"]})
    return cases, [record(path) for path in (plan_path, report_path, admission_path)]


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require frozen scientific source and a fresh output namespace")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    protocol = source / PROTOCOL
    spec = json.loads(protocol.read_text())
    if spec["arms"] != ["off", "local_joint_v1"] or spec["parameters"] != PARAMETERS:
        raise ValueError("Unexpected local-rescue protocol")
    cases, evidence = selected_cases(repo)
    inputs = [record(protocol), *evidence]
    for case in cases:
        inputs.extend([*case["fastas"], case["candidate"], case["constraints"], *case["golden"].values(),
                       *[record(path) for path in sorted(Path(case["checkpoint"]).iterdir()) if path.is_file()],
                       *[record(shutil.which(case[key]) or case[key]) for key in ("aligner", "tree_builder")]])
    inputs = list({pin["path"]: pin for pin in inputs}.values())
    source_pins = sources(source) + [record(source / "benchmark_tools/run_local_rescue_probe.slurm")]
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=source)
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Frozen source differs from committed bytes: " + str(relative))
    original = repo / "benchmarks/work/full_search_hillclimb_r1_v3_20261012/orthohmm/search/csrc"
    for name in ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so"):
        if record(source / "orthohmm/search/csrc" / name)["sha256"] != record(original / name)["sha256"]:
            raise ValueError("Require unchanged original CPU libraries")
    scientific_runtime = runtime()
    for pin in (*inputs, *source_pins, *scientific_runtime["files"]):
        check(pin)
    result = {"status": "prepared_label_free_probe_unrun", "source_root": str(source), "source_revision": revision,
              "repo": str(repo), "output": str(output), "cpu": 16, "cases": cases, "protocol": record(protocol),
              "inputs": inputs, "sources": source_pins, "runtime": scientific_runtime,
              "profile_runtime": require_profile_runtime(source), "tool_path": os.environ["PATH"],
              "scope": "Label-free HMM/local-tree applicability only; no benchmark scoring or final native inference. Shared-host incremental timing, not isolated full-pipeline runtime.",
              "defaults_promoted": False, "accuracy_evaluated": False, "goal_completed": False}
    validate_plan(result)
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def context_baseline_pairs(path, contexts, gene_species):
    context_by_gene = {gene: row["context_id"] for row in contexts for key in ("genes_a", "genes_b") for gene in row[key]}
    if len(context_by_gene) != sum(len(row[key]) for row in contexts for key in ("genes_a", "genes_b")):
        raise ValueError("Selected local contexts overlap")
    pairs = {row["context_id"]: [] for row in contexts}
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            a, b = row["gene_a"], row["gene_b"]
            left, right = context_by_gene.get(a), context_by_gene.get(b)
            if left is None and right is None:
                continue
            if (left != right or a == b or gene_species.get(a) != row["species_a"]
                    or gene_species.get(b) != row["species_b"] or row["species_a"] == row["species_b"]
                    or row["confidence"] not in ("high", "medium", "low")):
                raise ValueError("Baseline pairs violate complete local-context ownership")
            if row["confidence"] == "high":
                pairs[left].append((a, b))
    return pairs


def probe_case(plan, case, folder):
    folder.mkdir(parents=True, exist_ok=False)
    result_path = folder / "probe.json"
    report = {"status": "running", "dataset": case["dataset"], "contexts": [], "accuracy_evaluated": False,
              "defaults_promoted": False, "final_native_inference": False}
    write_json(result_path, report)
    started, before = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
    try:
        checkpoint = Path(case["checkpoint"])
        check(case["checkpoint_manifest"])
        if record(checkpoint / "manifest.json") != case["checkpoint_manifest"]:
            raise ValueError("Wrong canonical initial-hit checkpoint")
        names, species, q, t, scores = load_accuracy_checkpoint(checkpoint, verify=True)
        fasta = Path(case["fasta_directory"])
        files = [Path(pin["path"]).name for pin in case["fastas"]]
        sequences, gene_species = _load_sequence_data(str(fasta), files)
        filenames = {canonical_species_name(filename): filename for filename in files}
        owners = {gene: filenames[taxon] for gene, taxon in gene_species.items()}
        report["ownership"] = validate_ownership(names, species, owners)
        if [len(names), len(set(owners.values()))] != case["universe"]:
            raise ValueError("Probe does not cover the complete dataset")
        groups = read_partition(Path(case["candidate"]["path"]), set(names))
        indices = {gene: index for index, gene in enumerate(names)}
        clusters = [sorted(indices[gene] for gene in group) for group in groups]
        nomination = nominate_local_contexts(clusters, names, species, q, t, scores)
        write_json(folder / "nomination.json", nomination)
        report["nomination"] = record(folder / "nomination.json")
        report["nomination_counts"] = nomination["counts"]
        write_json(result_path, report)
        contexts = nomination["contexts"]
        constraints = load_membership_constraints(Path(case["constraints"]["path"]), Path(case["candidate"]["path"]))
        database = load_global_sequence_database(str(fasta), files, names) if contexts else None
        pairs = context_baseline_pairs(Path(case["golden"][PRODUCTS[2]]["path"]), contexts, gene_species) if contexts else {}
        config = PhylogenyConfig(mode="reconcile", species_tree_mode="infer", aligner=case["aligner"],
            tree_builder=case["tree_builder"], root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
            species_tree_rooting="min_variance", satellite_membership_rule="high_confidence_pair", gene_tree_support_mode="legacy")
        proposals = []
        for context in contexts:
            directory = folder / "contexts" / context["context_id"]
            row = {"context_id": context["context_id"], "status": "running", "directory": str(directory)}
            report["contexts"].append(row)
            write_json(result_path, report)
            try:
                proposal = evaluate_local_context(context, names, database, gene_species, sequences,
                    case["golden"]["species_tree.rooted.nwk"]["path"], files, config, directory, constraints,
                    pairs[context["context_id"]], cpu=plan["cpu"])
            except Exception:
                row["status"] = "failed"
                if (directory / "proposal.json").is_file():
                    row["proposal"] = record(directory / "proposal.json")
                raise
            proposals.append(proposal)
            row.update({"status": proposal["status"], "proposal": record(directory / "proposal.json"),
                "profile_accepted": proposal["profile_validation"]["accepted"],
                "tree_inferred": proposal["tree_inferred"], "accepted": proposal["accepted"]})
            write_json(result_path, report)
        working = folder / "candidate/orthohmm_working_res"
        candidate = working / "orthohmm_edges_clustered.txt"
        report["effective_partition"] = write_effective_candidates(case["candidate"]["path"], names, proposals, candidate, constraints)
        shutil.copyfile(case["constraints"]["path"], working / "phylogeny_candidate_merges.json")
        load_membership_constraints(working / "phylogeny_candidate_merges.json", candidate)
        report["candidate"] = record(candidate)
        report["constraints"] = record(working / "phylogeny_candidate_merges.json")
        report["tree_contexts"] = sum(row["tree_inferred"] for row in report["contexts"])
        report["accepted_contexts"] = sum(row["accepted"] for row in report["contexts"])
        report["status"] = ("analyzable_contexts_found" if report["tree_contexts"] else "no_analyzable_contexts")
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        after = resource.getrusage(resource.RUSAGE_SELF)
        report["parent_resources"] = {"wall_s": time.monotonic() - started,
            "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
            "lifetime_peak_rss_kib": after.ru_maxrss, "scope": "Parent process only; not a simultaneous process-tree peak"}
        write_json(result_path, report)
    return report


def validate_plan(plan):
    if (plan["cpu"] != 16 or [case["dataset"] for case in plan["cases"]] != list(UNIVERSES)
            or any(case["universe"] != list(UNIVERSES[case["dataset"]]) for case in plan["cases"])
            or plan["status"] != "prepared_label_free_probe_unrun"):
        raise ValueError("Unexpected local-rescue dataset order, universe or resource budget")
    if "inputs" not in plan:
        raise ValueError("Missing authenticated probe inputs")
    for case in plan["cases"]:
        pins = [*case["fastas"], case["candidate"], case["constraints"], case["checkpoint_manifest"], *case["golden"].values()]
        if (set(case["golden"]) != {*PRODUCTS, "species_tree.rooted.nwk"}
                or any(pin not in plan["inputs"] for pin in pins)
                or any(Path(pin["path"]).parent != Path(case["fasta_directory"]) for pin in case["fastas"])
                or Path(case["checkpoint_manifest"]["path"]) != Path(case["checkpoint"]) / "manifest.json"
                or any(Path(pin["path"]).parent != Path(case["baseline_output"]) / "orthohmm_phylogeny"
                       for pin in case["golden"].values())):
            raise ValueError("Probe case paths are not bound to authenticated inputs")


def run(path):
    plan = json.loads(path.read_text())
    validate_plan(plan)
    source, output = (Path(plan[key]) for key in ("source_root", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_local_rescue_probe.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve every existing attempt")
    pins = [record(path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
    report = {"status": "running", "execution_plan": record(path), "cases": [],
              "job_id": os.environ.get("SLURM_JOB_ID"), "source_revision": plan["source_revision"],
              "scope": plan["scope"], "host_before": host_activity(), "accuracy_evaluated": False, "defaults_promoted": False}
    write_json(output / "results.json", report)
    for case in plan["cases"]:
        report["active_dataset"] = case["dataset"]
        write_json(output / "results.json", report)
        folder = output / case["dataset"]
        command = [sys.executable, str(Path(__file__).resolve()), "case", "--execution-plan", str(path), "--dataset", case["dataset"]]
        execution = timed_stage(command, output, "probe_" + case["dataset"], source, env)
        row = {"dataset": case["dataset"], "execution": execution}
        report["cases"].append(row)
        if (folder / "probe.json").is_file():
            result = json.loads((folder / "probe.json").read_text())
            row.update({"report": record(folder / "probe.json"), "status": result["status"],
                        "tree_contexts": result.get("tree_contexts"), "accepted_contexts": result.get("accepted_contexts")})
        else:
            row["status"] = "failed_before_report"
        write_json(output / "results.json", report)
        if execution["returncode"] or row["status"] not in ("analyzable_contexts_found", "no_analyzable_contexts"):
            report["status"] = "failed_probe_no_fallback_on_failure"
            break
        if row["status"] == "analyzable_contexts_found":
            report["status"] = ("applicability_found_accuracy_evaluation_pending" if row["accepted_contexts"]
                                else "analyzable_but_zero_admitted_no_accuracy_claim")
            break
    else:
        report["status"] = "all_three_datasets_inactive_no_accuracy_claim"
    try:
        for pin in pins:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Scientific runtime changed during the probe")
    except Exception as error:
        report.update(status="failed_final_verification", error_type=type(error).__name__, error=str(error))
    report["active_dataset"] = None
    report["unrun_datasets"] = [case["dataset"] for case in plan["cases"]
                               if case["dataset"] not in {row["dataset"] for row in report["cases"]}]
    report["host_after"] = host_activity()
    write_json(output / "results.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    for name in ("run", "case"):
        command = commands.add_parser(name)
        command.add_argument("--execution-plan", required=True, type=Path)
        if name == "case":
            command.add_argument("--dataset", choices=UNIVERSES, required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
    elif args.command == "run":
        result = run(args.execution_plan.resolve())
    else:
        plan = json.loads(args.execution_plan.read_text())
        validate_plan(plan)
        if Path(__file__).resolve() != Path(plan["source_root"]) / "benchmark_tools/run_local_rescue_probe.py":
            raise ValueError("Execute the frozen local-rescue worker")
        for pin in [*plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]:
            check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Changed scientific worker runtime or CUDA visibility")
        coordinator = json.loads((Path(plan["output"]) / "results.json").read_text())
        if (coordinator["status"] != "running" or coordinator["active_dataset"] != args.dataset
                or coordinator["execution_plan"] != record(args.execution_plan)):
            raise ValueError("Dataset worker is not the coordinator's active priority case")
        case = next(case for case in plan["cases"] if case["dataset"] == args.dataset)
        result = probe_case(plan, case, Path(plan["output"]) / args.dataset)
    print(result["status"])
    return 1 if result["status"].startswith("failed") else 0


if __name__ == "__main__":
    raise SystemExit(main())
