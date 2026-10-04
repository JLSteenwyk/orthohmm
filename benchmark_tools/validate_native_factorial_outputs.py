"""Check completed full-native factorial outputs, not accuracy or timing admission.

The production entry point requires the original request and a fresh terminal
scheduler observation. The semantic kernel also accepts explicit fixture
contexts for tests; that does not authorize production execution or admission.
"""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo, SeqIO
import numpy as np

from benchmark_tools.native_factorial_adapter import CORE_SHA, factors
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_native_factorial_cost import ROOT, SCOPE, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.simulation_method_outputs import orthohmm_pairs
from benchmark_tools.validate_factorial_partition import validate_partition


SATELLITE = dict(max_component_genes=500, max_satellite_genes=12,
    max_satellite_to_anchor_ratio=.75, min_margin=1.5, iteration_margin_increment=.5,
    max_satellites_per_anchor=4, max_species_overlap_fraction=1., min_avg_score=0.,
    min_max_score=0., min_coverage=.5, min_norm=.03, max_iterations=2)
CHECKPOINT_FILES = {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy",
                    "hit_targets.npy", "hit_scores.npy"}
TOOL_VERSIONS = {"aligner": "v7.525 (2024/Mar/13)",
                 "tree_builder": "FastTree Version 2.2.0 Double precision"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, name):
    require(type(value) is int and value >= 0, "Invalid integer count: " + name)
    return value


class Evidence:
    def __init__(self):
        self.files = {}

    def bind(self, path):
        path = Path(path)
        require(path.is_absolute() and path.resolve() == path and not path.is_symlink(),
                "Require direct absolute native evidence path")
        ref = record(path)
        require(path.is_file(), "Require a native evidence file")
        require(str(path) not in self.files or self.files[str(path)] == ref,
                "Native evidence changed during review")
        self.files[str(path)] = ref
        return path

    def json(self, path):
        return json.loads(self.bind(path).read_text())

    def finish(self):
        for ref in self.files.values():
            check(ref)
        return list(self.files.values())


def expected_factors(cell):
    flags = factors(cell)
    original = dict(name="high_sensitivity", kmer_k=4, max_candidates_per_query=100,
                    multipass_graph=True, profile_expansion=True, leiden_seed=4)
    return dict(flags, frozen_pipeline_sha256=CORE_SHA, original_accuracy_profile=original,
        active_accuracy_profile=dict(original, profile_expansion=flags["profile_expansion"]),
        candidate_gate_decoupled=flags["candidate_expansion"] and not flags["reconciliation"],
        production_source_modified=False, production_globals_modified=False)


def universe(context, evidence):
    originals = context["inputs"]
    by_name = {Path(ref["path"]).name: ref for ref in originals}
    order = context["native_order"]
    require(len(by_name) == len(originals) == context["proteomes"]
            and len(order) == len(set(order)) == len(by_name) and set(order) == set(by_name),
            "Input inventory or native enumeration differs")
    source = Path(context["input_directory"])
    entries = list(source.iterdir())
    require({p.name for p in entries} == set(by_name), "Prepared input inventory differs")
    owners, counts, digest = {}, {}, hashlib.sha256()
    for name in order:
        check(by_name[name])
        original = evidence.bind(by_name[name]["path"])
        copied = evidence.bind(source / name)
        require(all(evidence.files[str(original)][k] == evidence.files[str(copied)][k]
                    for k in ("bytes", "sha256")), "Prepared input bytes differ")
        taxon = Path(name).stem
        require(taxon not in counts, "Duplicate input taxon")
        counts[taxon] = 0
        for sequence in SeqIO.parse(copied, "fasta"):
            require(sequence.id and sequence.id not in owners and sequence.seq,
                    "Empty sequence or duplicate input gene")
            owners[sequence.id] = taxon
            counts[taxon] += 1
            digest.update(json.dumps([name, sequence.id], separators=(",", ":")).encode() + b"\n")
        require(counts[taxon] > 0, "Empty proteome")
    require(len(owners) == context["genes"], "Input gene count differs")
    return owners, digest.hexdigest(), {name: counts[Path(name).stem] for name in order}


def checkpoint(path, owners, metrics, evidence, chunk_size=1000000):
    require(type(chunk_size) is int and chunk_size > 0, "Invalid checkpoint chunk size")
    manifest = evidence.json(path / "manifest.json")
    require(manifest.get("schema_version") == 1 and manifest.get("complete") is True
            and manifest.get("genes") == len(owners) and set(manifest["files"]) == CHECKPOINT_FILES,
            "Checkpoint completion, universe or inventory differs")
    require({p.name for p in path.iterdir()} == CHECKPOINT_FILES | {"manifest.json"},
            "Checkpoint on-disk inventory differs")
    hits = integer(manifest["hits"], "checkpoint hits")
    require(hits <= metrics["counts"]["significant_hits"], "Checkpoint exceeds significant search hits")
    for name, pin in manifest["files"].items():
        file = evidence.bind(path / name)
        require(pin == {k: evidence.files[str(file)][k] for k in ("bytes", "sha256")},
                "Checkpoint checksum differs")
    names = (path / "gene_names.txt").read_text().splitlines()
    require(names == sorted(owners), "Checkpoint gene order differs from frozen lexical order")
    species_ids = {}
    expected = []
    for gene in names:
        species_ids.setdefault(owners[gene], len(species_ids))
        expected.append(species_ids[owners[gene]])
    specifications = {"gene_to_species.npy": (len(names), np.dtype("int32")),
        "hit_queries.npy": (hits, np.dtype("int32")), "hit_targets.npy": (hits, np.dtype("int32")),
        "hit_scores.npy": (hits, np.dtype("float64"))}
    # Memory-map large hit files and check bounded slices; never allocate an NxN graph.
    for name, (length, dtype) in specifications.items():
        values = np.load(path / name, mmap_mode="r", allow_pickle=False)
        require(values.shape == (length,) and values.dtype == dtype, "Checkpoint shape/dtype differs: " + name)
        for start in range(0, length, chunk_size):
            block = values[start:start + chunk_size]
            if name == "gene_to_species.npy":
                require(np.array_equal(block, expected[start:start + chunk_size]), "Checkpoint gene ownership differs")
            elif name == "hit_scores.npy":
                require(np.isfinite(block).all(), "Nonfinite checkpoint hit score")
            else:
                require(((block >= 0) & (block < len(names))).all(), "Checkpoint hit endpoint outside universe")
        del values
    return dict(genes=len(names), hits=hits, order="lexical_gene_ids",
                validation="mmap_bounded_slices", exact_historical_hit_equivalence_established=False)


def groups(path, format, owners, evidence):
    path = evidence.bind(path)
    if format == "space_separated_groups":
        value = {f"Family{i:07d}": line.split() for i, line in
                 enumerate(line for line in path.read_text().splitlines() if line.strip())}
    else:
        value = read_predictions(path, format)
    require(set(membership(value)) == set(owners), "Native partition loses or adds input genes")
    return value


def candidate(meta, counts, working, owners, evidence):
    profile = meta["phylogeny_candidate_profile"]
    paths = dict(candidate_checkpoint=working / "phylogeny_candidate_superfamilies.txt",
                 seed_sidecar=working / "phylogeny_candidate_seeds.tsv",
                 merge_trace_sidecar=working / "phylogeny_candidate_merges.json")
    require(profile.get("profile") == "satellite_v2" and profile.get("parameters") == SATELLITE
            and profile.get("membership_policy") == "high_confidence_pair"
            and all(profile.get(k) == str(v) for k, v in paths.items()), "Candidate settings or paths differ")
    partitions = groups(paths["candidate_checkpoint"], "space_separated_groups", owners, evidence)
    require(integer(profile["candidate_families"], "candidate families") == len(partitions)
            and integer(profile["seed_families"], "seed families") == counts["phylogeny_seed_families"]
            and integer(profile["merges"], "candidate merges") == counts["phylogeny_candidate_merges"]
            and profile["seed_families"] - profile["merges"] == len(partitions)
            and integer(profile["iterations"], "candidate iterations") <= 2,
            "Candidate completion counts differ")
    seeds = []
    with evidence.bind(paths["seed_sidecar"]).open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        require(rows.fieldnames == ["candidate_family", "seed_families"], "Candidate seed header differs")
        names = []
        for row in rows:
            require(None not in row and all(v is not None for v in row.values()), "Malformed seed sidecar")
            names.append(row["candidate_family"])
            ids = row["seed_families"].split(",")
            require(ids and ids == sorted(set(ids)), "Invalid candidate seed list")
            seeds.extend(ids)
    require(names == list(partitions) and sorted(seeds) ==
            [f"Seed{i:07d}" for i in range(profile["seed_families"])], "Candidate seed coverage differs")
    trace = evidence.json(paths["merge_trace_sidecar"])
    require(isinstance(trace, list) and len(trace) == profile["merges"], "Candidate trace count differs")
    candidate_owner = membership(partitions)
    for item in trace:
        left, right = item["source_genes"], item["target_genes"]
        require(left and right and len(set(left)) == len(left) and len(set(right)) == len(right)
                and not set(left) & set(right) and set(left) | set(right) <= set(owners),
                "Invalid candidate constraint endpoints")
        require(len({candidate_owner[g] for g in (*left, *right)}) == 1,
                "Candidate constraint crosses candidate families")
    return paths["candidate_checkpoint"], partitions, len(trace)


def phylogeny(output, owners, counts, context, evidence, candidate_path, constraints):
    directory = output / "orthohmm_phylogeny"
    native = evidence.json(directory / "provenance_manifest.json")
    summary = evidence.json(directory / "reconciliation_summary.json")
    expected = dict(mode="reconcile", cpu_budget=context["cpu"], species_tree_mode="infer",
        species_tree_rooting="min_variance", root_duplication_rule="species_overlap",
        pair_orthology_rule="positive_paralogy",
        species_tree_source="internally_inferred_from_orthohmm_single_copy_families")
    require(all(native.get(k) == v for k, v in expected.items()), "Native phylogeny settings differ")
    require(native["results"] == summary and summary["schema_version"] == 1
            and summary["output_directory"] == str(directory), "Native phylogeny completion records differ")
    for key, value in summary.items():
        if key not in {"schema_version", "output_directory"}:
            require(counts.get("phylogeny_" + key) == value, "Phylogeny metric count differs: " + key)
    require(summary["checkpoint_hits"] == summary["remapped_checkpoint_hits"] == 0
            and summary["species_tree_checkpoint_hit"] is False
            and native["species_tree_inference"]["checkpoint_hit"] is False,
            "Native phylogeny reused checkpoints")
    inputs = [{"filename": name, "taxon": Path(name).stem,
               "sha256": evidence.files[str(Path(context["input_directory"]) / name)]["sha256"]}
              for name in context["native_order"]]
    require(native["input_proteomes"] == inputs, "Phylogeny input provenance differs")
    for role in TOOL_VERSIONS:
        require(native["tools"][role] == dict(path=context[role], version=TOOL_VERSIONS[role]),
                "Native phylogeny tool provenance differs")
    tree_path = evidence.bind(directory / "species_tree.rooted.nwk")
    require(native["species_tree_sha256"] == evidence.files[str(tree_path)]["sha256"], "Species tree checksum differs")
    tree = Phylo.read(tree_path, "newick")
    leaves = [leaf.name for leaf in tree.get_terminals()]
    taxa = sorted(set(owners.values()))
    require(sorted(leaves) == native["species_tree_taxa"] == taxa, "Species tree taxon coverage differs")
    require(all(node.branch_length is None or math.isfinite(node.branch_length)
                for node in tree.find_clades()), "Nonfinite species tree branch")
    root_path = directory / "orthohmm_root_hogs.tsv"
    roots = groups(root_path, "root_hogs", owners, evidence)
    if candidate_path is not None:
        require(native["input_cluster_sha256"] == evidence.files[str(candidate_path)]["sha256"],
                "Phylogeny candidate hash differs")
        validate_partition(candidate_path, root_path, set(owners), summary)
    else:
        # C-off reconciliation overwrites its input. Reconstruct only the
        # canonical original families from source labels, not inferred truth.
        original = defaultdict(list)
        with root_path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                original[row["source_family"]].extend(row["genes"].split(","))
        require(sorted(original) == [f"Family{i:07d}" for i in range(summary["candidate_families"])],
                "Unknown phylogeny source family")
        payload = "".join(" ".join(sorted(original[name])) + "\n" for name in sorted(original))
        require(hashlib.sha256(payload.encode()).hexdigest() == native["input_cluster_sha256"],
                "Reconstructed canonical candidate hash differs")
    require(list(roots) == [f"RootHOG{i:07d}" for i in range(len(roots))]
            and len(roots) == counts["phylogeny_root_hogs"]
            and summary["bypassed_families"] + summary["reconciled_families"] == summary["candidate_families"],
            "Phylogeny family or root completion differs")
    membership_review = native["membership_reconciliation"]
    if constraints:
        require(isinstance(membership_review, dict) and membership_review.get("policy") == "high_confidence_pair"
                and membership_review.get("constraints") == constraints
                and membership_review.get("supported_constraints", -1) +
                    membership_review.get("detached_constraints", -1) == constraints,
                "Phylogeny membership constraint accounting differs")
    else:
        require(membership_review is None, "Unexpected phylogeny membership filtering")
    pair_path = evidence.bind(directory / "orthohmm_pairwise_orthologs.tsv")
    previous, pair_count = None, 0
    for pair in orthohmm_pairs(pair_path, owners):
        require(pair[0] < pair[1] and (previous is None or previous < pair),
                "Native pairs are not canonical, unique and sorted")
        previous, pair_count = pair, pair_count + 1
    require(pair_count == counts["phylogeny_ortholog_pairs"], "Native pair count differs")
    return roots, dict(root_hogs=len(roots), native_pair_rows=pair_count, membership=membership_review)


def validate_semantics(context):
    """Validate exact caller-bound context; no execution or resource authorization."""
    evidence = Evidence()
    owners, owner_digest, per_species = universe(context, evidence)
    root = Path(context["output_root"])
    output, working = root / "native", root / "native/orthohmm_working_res"
    metrics = evidence.json(root / "metrics.json")
    flags = expected_factors(context["cell"])
    meta, counts = metrics["metadata"], metrics["counts"]
    common = dict(accuracy_profile="high_sensitivity", clustering="leiden", cpm_resolution=.1,
        cpu_budget=context["cpu"], evalue_threshold=.0001, leiden_seed=4, search_kmer_k=4,
        search_max_candidates_per_query=100, search_mode="builtin", search_total_threads=context["cpu"],
        search_threads_per_worker=context["threads_per_worker"],
        search_workers=context["cpu"] // context["threads_per_worker"], substitution_matrix="BLOSUM62",
        fasta_directory=context["input_directory"], output_directory=str(output), native_factorial=flags,
        high_sensitivity_checkpoint=str(working / "high_sensitivity_checkpoint"))
    phylogeny_settings = dict(phylogeny="reconcile",
        phylogeny_candidates="satellite_v2" if flags["candidate_expansion"] else "seed",
        species_tree_mode="infer", species_tree_rooting="min_variance", species_tree=None,
        phylogeny_root_rule="species_overlap", phylogeny_pair_rule="positive_paralogy",
        aligner=context["aligner"], tree_builder=context["tree_builder"])
    if flags["reconciliation"]:
        common.update(phylogeny_settings)
    else:
        require(not set(phylogeny_settings) & set(meta), "Unexpected reconciliation settings in R-off cell")
    require(metrics.get("status") == "complete" and metrics.get("command") == context["command"]
            and metrics.get("cwd") == context["cwd"] and all(meta.get(k) == v for k, v in common.items()),
            "Native completion, command or frozen settings differ")
    stages = {"search", "edge_thresholds", "network_edges", "clustering", "refinement", "orthogroup_materialization"}
    for enabled, stage in ((flags["profile_expansion"], "profile_expansion"),
                           (flags["candidate_expansion"], "phylogeny_candidates"),
                           (flags["reconciliation"], "phylogeny")):
        if enabled:
            stages.add(stage)
    require(set(metrics["stages"]) == stages, "Native stage set differs")
    for stage in metrics["stages"].values():
        require(type(stage.get("wall_s")) in (int, float) and math.isfinite(stage["wall_s"])
                and stage["wall_s"] > 0, "Invalid native stage duration")
    for name, value in counts.items():
        if name == "phylogeny_species_tree_checkpoint_hit":
            require(value is False, "Native species tree checkpoint reused")
        else:
            integer(value, name)
    require(counts["genes"] == len(owners) and counts["species"] == context["proteomes"], "Native universe counts differ")
    require(flags["candidate_expansion"] == ("phylogeny_candidate_profile" in meta), "Unexpected candidate metadata")
    require(flags["profile_expansion"] == ("high_sensitivity_profiles" in counts), "Unexpected profile counts")
    checkpoint_review = checkpoint(working / "high_sensitivity_checkpoint", owners, metrics, evidence)
    materialized = groups(output / "orthohmm_orthogroups.txt", "named_groups", owners, evidence)
    require(len(materialized) == counts["orthogroups"], "Materialized group count differs")
    clustered = groups(working / "orthohmm_edges_clustered.txt", "space_separated_groups", owners, evidence)
    require({frozenset(g) for g in materialized.values()} == {frozenset(g) for g in clustered.values()},
            "Materialized partition differs from final clustering")
    candidate_path, constraint_count = None, 0
    if flags["candidate_expansion"]:
        candidate_path, candidates, constraint_count = candidate(meta, counts, working, owners, evidence)
        if not flags["reconciliation"]:
            require({frozenset(g) for g in candidates.values()} == {frozenset(g) for g in clustered.values()},
                    "Candidate-only partition differs from final clustering")
    phylogeny_review = None
    if flags["reconciliation"]:
        roots, phylogeny_review = phylogeny(output, owners, counts, context, evidence, candidate_path, constraint_count)
        require({frozenset(g) for g in roots.values()} == {frozenset(g) for g in materialized.values()},
                "Materialized partition differs from phylogenetic root HOGs")
    else:
        require(not (output / "orthohmm_phylogeny").exists()
                and not any(k.startswith("phylogeny_") and k not in
                    {"phylogeny_seed_families", "phylogeny_candidate_merges"} for k in counts),
                "Unexpected reconciliation output in R-off cell")
    return dict(schema="native_factorial_output_review_v1", status="native_outputs_validated",
        cell=context["cell"], factors=flags, input_genes=len(owners), orthogroups=len(materialized),
        gene_ownership_sha256=owner_digest, per_species_counts=per_species,
        checkpoint=checkpoint_review, phylogeny=phylogeny_review, checked_files=evidence.finish(),
        native_outputs_validated=True, accuracy_evaluated=False, resource_measurements_admitted=False,
        next_identity_authorized=False, source=record(__file__),
        limitations=["Semantic output review is not accuracy evaluation, runtime integrity or resource admission.",
            "Species tree provenance/coverage are checked; no independent tree reconstruction is claimed.",
            "Stage durations are diagnostics, not full native resource costs or causal component overheads."])


def validate(request_ref):
    request = read(request_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    runs = validate_plan(plan)
    validate_request(request, plan_ref, request["job_id"])
    terminal = verify_terminal(request["job_id"])
    state = terminal["verified"].get("fields", terminal["verified"])
    require(state.get("JobState", state.get("State")) == "COMPLETED"
            and state["ExitCode"] == "0:0", "Require successful terminal allocation for completed output review")
    run = runs[request["index"]]
    baseline = read(plan["baseline"])
    evidence = [request_ref, plan_ref, plan["baseline"], *plan["helper_sources"], *run["inputs"]]
    for ref in evidence:
        check(ref)
    python = baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    command = [python, str(ROOT / "benchmark_tools/run_native_factorial_cost.py"), "--native", "--plan",
               plan_ref["path"], "--plan-sha256", plan_ref["sha256"], "--index", str(run["index"])]
    context = dict(run, input_directory=str(Path(run["output_root"]) / "input"), cpu=32,
        threads_per_worker=4, command=command, cwd=baseline["core_root"],
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    result = validate_semantics(context)
    preparation_ref = record(Path(run["output_root"]) / "preparation.json")
    preparation = read(preparation_ref)
    execution_ref = record(Path(run["output_root"]) / "native_execution.json")
    execution = read(execution_ref)
    require(preparation["status"] == "fresh_factorial_inputs_prepared"
            and preparation["gene_ownership_sha256"] == result["gene_ownership_sha256"]
            and preparation["per_species_counts"] == result["per_species_counts"]
            and preparation["genes"] == run["genes"], "Preparation gene provenance differs")
    require(execution["status"] == "native_factorial_completed_pending_output_review"
            and execution["plan"] == plan_ref and execution["index"] == run["index"]
            and execution["cell"] == run["cell"] and execution["factors"] == result["factors"]
            and execution["native_order"] == run["native_order"]
            and execution["automatic_retry"] is False, "Native execution receipt differs")
    evidence += [preparation_ref, execution_ref]
    for ref in evidence:
        check(ref)
    return dict(result, index=run["index"], job_id=request["job_id"], plan=plan_ref, request=request_ref,
        scheduler=terminal, execution_scope=SCOPE, evidence=evidence, terminal_reviewed=False,
        terminal_scheduler_confirmed=True, uncontended_timing=False,
        contention_distortion="unknown_potentially_tool_dependent")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output review already exists")
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Request checksum differs")
    report = validate(ref)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
