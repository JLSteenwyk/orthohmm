"""Trace every frozen OrthoBench reference family through retained checkpoints."""

import argparse
from collections import Counter
import csv
import itertools
import json
import math
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_ob_stratified_errors import STRATA_SHA, validate_strata
from benchmark_tools.assemble_orthobench_factorial import load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_orthobench_factorial import indexed_species
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.replay_high_sensitivity import load_hits
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.score_ygob_groups import membership, read_predictions

FACTORIAL_SHA = "6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141"
STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined", "candidates", "root_hogs")
TRANSITIONS = (("multipass", "multipass_refined", "profile_off_refinement"),
               ("multipass", "strict_profiles", "profile_expansion_and_reclustering"),
               ("strict_profiles", "strict_profiles_refined", "profile_on_refinement"),
               ("multipass_refined", "strict_profiles_refined", "matched_refined_branch_comparison_not_direct_step"),
               ("strict_profiles_refined", "candidates", "candidate_expansion"),
               ("candidates", "root_hogs", "tree_inference_reconciliation_and_constraints"))


def reference_inventory(references):
    owners = {}
    for name, genes in sorted(references.items()):
        for gene in sorted(genes):
            owners.setdefault(gene, []).append(name)
    return {"families": len(references), "memberships": sum(map(len, references.values())),
            "unique_genes": len(owners), "shared_genes": {g: names for g, names in owners.items() if len(names) > 1}}


def partition(path, format_name, universe):
    if format_name == "plain":
        with path.open() as handle:
            groups = {str(i): line.split() for i, line in enumerate(handle) if line.strip()}
    else:
        groups = read_predictions(path, format_name)
    index = membership(groups)
    if set(index) != universe:
        raise ValueError("Checkpoint does not partition the complete input universe")
    return {key: set(value) for key, value in groups.items()}, index


def family_group_summary(genes, groups, index, reference_universe):
    touched = sorted({index[g] for g in genes})
    parts = [groups[key] & genes for key in touched]
    return {"groups_touching_family": len(parts), "within_family_pairs": sum(len(p) * (len(p) - 1) // 2 for p in parts),
            "largest_family_component": max(map(len, parts)),
            "cross_reference_pairs_incident": sum(len(groups[key] & genes) * len((groups[key] & reference_universe) - genes) for key in touched),
            "unlabelled_pairs_incident": sum(len(groups[key] & genes) * len(groups[key] - reference_universe) for key in touched),
            "groups": [{"id": key, "total_genes": len(groups[key]), "family_genes": sorted(groups[key] & genes),
                        "other_reference_genes": sorted((groups[key] & reference_universe) - genes),
                        "unlabelled_gene_count": len(groups[key] - reference_universe)} for key in touched]}


def pair_trace(genes, hits, indices, owners):
    rows = []
    for left, right in itertools.combinations(sorted(genes), 2):
        forward, reverse = hits.get((left, right)), hits.get((right, left))
        if any(value is not None and (not math.isfinite(value) or value <= 0) for value in (forward, reverse)):
            raise ValueError("Invalid retained normalized hit score")
        row = {"left": left, "right": right, "same_species": owners[left] == owners[right],
               "forward_normalized_hit": forward, "reverse_normalized_hit": reverse}
        row.update({stage: indices[stage][left] == indices[stage][right] for stage in STAGES})
        rows.append(row)
    return rows


def transitions(rows):
    result = {}
    for before, after, _ in TRANSITIONS:
        counts = Counter("retained" if row[before] and row[after] else "lost" if row[before]
                         else "gained" if row[after] else "absent_both" for row in rows)
        result[before + "_to_" + after] = {key: counts[key] for key in ("retained", "lost", "gained", "absent_both")}
    return result


def validate_merge_reconstruction(events, seeds, candidates, reference_genes):
    index = membership({k: list(v) for k, v in seeds.items()})
    universe = set(index)
    current = {key: set(value) for key, value in seeds.items()}
    parent = {key: key for key in seeds}
    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key
    retained = []
    previous_iteration = -1
    for i, event in enumerate(events):
        source, target = (set(event[key]) for key in ("source_genes", "target_genes"))
        if (not source or not target or source & target or event["iteration"] < previous_iteration
                or len(source) != len(event["source_genes"]) or len(target) != len(event["target_genes"])
                or len(source) != event["source_size"] or len(target) != event["target_size"]
                or not (source | target) <= universe):
            raise ValueError("Malformed or unordered merge evidence")
        previous_iteration = event["iteration"]
        roots = [{find(index[g]) for g in side} for side in (source, target)]
        if any(len(value) != 1 for value in roots) or roots[0] == roots[1]:
            raise ValueError("Merge sides are not distinct current components")
        a, b = next(iter(roots[0])), next(iter(roots[1]))
        if (source | target) & reference_genes:
            retained.append({"event_index": i, "iteration": event["iteration"],
                             "source_reference_genes": sorted(source & reference_genes),
                             "target_reference_genes": sorted(target & reference_genes),
                             "source_size": len(source), "target_size": len(target),
                             "support": event["support"], "margin": event["margin"],
                             "note": "Logged sides describe iteration-start clusters, not necessarily current union components."})
        if len(current[a]) > len(current[b]):
            a, b = b, a
        current[b].update(current.pop(a))
        parent[a] = b
    if {frozenset(value) for value in current.values()} != {frozenset(value) for value in candidates.values()}:
        raise ValueError("Logged merges do not reconstruct the retained candidate partition")
    return retained


def render(report):
    lines = ["# All-Family OrthoBench Stage Trace", "",
             "Descriptive reference co-membership counts, not official benchmark recall or independent biological evidence.",
             "Direct cached hits are not orthology predictions. Missing hits do not identify a prefilter or scoring failure.", "",
             "| RefOG | Genes | Cached Hit Pairs | Multipass | Refined | Profiles | Profiles Refined | Candidates | Root HOGs |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name, row in report["families"].items():
        numbers = [row["genes"], row["search"]["either_direction_pairs"],
                   *[row["stages"][stage]["within_family_pairs"] for stage in STAGES]]
        lines.append("| " + name + " | " + " | ".join(str(n) for n in numbers) + " |")
    lines += ["", "## Evidence Limits", "", *["- " + note for note in report["limitations"]]]
    return "\n".join(lines) + "\n"


def assemble(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    factorial = read_frozen(results / "orthobench_factorial_results_20260916.json", FACTORIAL_SHA)
    strata = read_frozen(results / "ob_error_strata_prepared_20260916.json", STRATA_SHA)
    validate_strata(strata)
    references, uncertain, _, reference_records = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    if set(references) != set(strata["families"]):
        raise ValueError("Reference family panel differs")
    reference_genes = set().union(*references.values())
    inventory = reference_inventory(references)
    records = [prepared["cache"], prepared["replay_verification"], *prepared["fasta_inputs"], *reference_records]
    for item in records:
        verify_file(Path(item["path"]), item)
    replay = json.loads(Path(prepared["replay_verification"]["path"]).read_text())
    if replay["status"] != "equivalent" or any(not row["byte_equal"] or not row["partition_equal"] for row in replay["stages"].values()):
        raise ValueError("Retained replay was not equivalent")
    sources = {stage: replay["stages"][stage]["output"] for stage in STAGES[:4]}
    sources.update(candidates=factorial["predictions"]["p1_c1_r0"], root_hogs=factorial["predictions"]["p1_c1_r1"])
    arm = prepared["candidate_arms"]["p1_c1"]
    records += [*sources.values(), arm["membership_constraints"]]
    owners = {}
    for item in prepared["fasta_inputs"]:
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in owners:
                raise ValueError("Duplicate input protein")
            owners[protein.id] = Path(item["path"]).name
    universe = set(owners)
    if not reference_genes <= universe:
        raise ValueError("Reference genes absent from inference inputs")
    groups, indices, scores = {}, {}, {}
    for stage in STAGES:
        item = sources[stage]
        verify_file(Path(item["path"]), item)
        groups[stage], indices[stage] = partition(Path(item["path"]), "root_hogs" if stage == "root_hogs" else "plain", universe)
        scores[stage] = score_partition(list(groups[stage].values()), references, uncertain)
    for stage, cell in (("multipass_refined", "p0_c0_r0"), ("strict_profiles_refined", "p1_c0_r0"),
                        ("candidates", "p1_c1_r0"), ("root_hogs", "p1_c1_r1")):
        if scores[stage]["refog_records"] != factorial["scores"][cell]["refog_records"]:
            raise ValueError("Fresh stage scores differ from frozen sufficient statistics")
    if any(len({indices["candidates"][gene] for gene in group}) != 1 for group in groups["root_hogs"].values()):
        raise ValueError("Root HOG crosses candidate-family boundaries")
    verify_file(Path(arm["membership_constraints"]["path"]), arm["membership_constraints"])
    events = json.loads(Path(arm["membership_constraints"]["path"]).read_text())
    merge_evidence = validate_merge_reconstruction(events, groups["strict_profiles_refined"], groups["candidates"], reference_genes)
    if len(events) != arm["expansion"]["merges"]:
        raise ValueError("Merge count differs from frozen candidate preparation")
    payload = load_hits(Path(prepared["cache"]["path"]))
    indexed_species(payload, owners)
    output.mkdir(parents=True)
    families = {}
    pair_path = output / "reference_pair_trace.tsv"
    with pair_path.open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=["refog", "left", "right", "same_species", "forward_normalized_hit", "reverse_normalized_hit", *STAGES], delimiter="\t")
        writer.writeheader()
        for name, genes in sorted(references.items()):
            pairs = pair_trace(genes, payload["all_hits"], indices, owners)
            for pair in pairs:
                writer.writerow({"refog": name, **{key: "NA" if value is None else value for key, value in pair.items()}})
            families[name] = {"genes": len(genes), "possible_pairs": len(pairs),
                "search": {"either_direction_pairs": sum(row["forward_normalized_hit"] is not None or row["reverse_normalized_hit"] is not None for row in pairs),
                           "both_direction_pairs": sum(row["forward_normalized_hit"] is not None and row["reverse_normalized_hit"] is not None for row in pairs)},
                "stages": {stage: family_group_summary(genes, groups[stage], indices[stage], reference_genes) for stage in STAGES},
                "transitions": transitions(pairs),
                "merge_event_indices": [event["event_index"] for event in merge_evidence
                    if genes & (set(event["source_reference_genes"]) | set(event["target_reference_genes"]))]}
    for item in records:
        verify_file(Path(item["path"]), item)
    report = {"status": "retained_stage_trace_complete", "publication_ready": False, "job_id": os.environ.get("SLURM_JOB_ID"),
        "source": file_provenance(Path(__file__)), "inputs": records, "families": families, "scores": scores,
        "transition_specification": TRANSITIONS,
        "reference_inventory": inventory,
        "frozen_snapshots": [file_provenance(results / name) for name in
            ("orthobench_factorial_prepared_20260916.json", "orthobench_factorial_results_20260916.json", "ob_error_strata_prepared_20260916.json")],
        "merge_events_total": len(events), "reference_incident_merge_events": merge_evidence,
        "candidate_partition_reconstructed": True, "pair_trace": file_provenance(pair_path),
        "illustrative_families_by_stratum": strata["illustrative_families_by_stratum"],
        "limitations": ["All70 development-exposed families retained; illustrations follow the frozen feature-based hash ranking.",
            "Reference memberships can overlap; shared genes retain every reference assignment. Family-level counts are not disjoint totals.",
            "Pair counts include within-species pairs and low-certainty members; these descriptive counts are not official benchmark precision or recall.",
            "Cross-reference and unlabelled incident pairs are raw membership descriptors, not official false-positive counts.",
            "The normalized-hit cache lacks rejected candidates and prefilter logs; absent hits cannot distinguish prefilter rejection from scoring rejection.",
            "Initial RBNH edges and added profile-hit edge identities are not traced here; grouping changes alone do not identify their causal edges.",
            "Root-HOG membership changes include tree inference, reconciliation and membership constraints; distinguishing their causes needs per-family tree and constraint inspection.",
            "Merge events are iteration-start cluster snapshots; union reconstruction validates final membership, not calibrated biological support.",
            "These traces do not supply independent duplication, fragment or domain annotations or replace the prespecified biological application."]}
    (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "results.md").write_text(render(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
