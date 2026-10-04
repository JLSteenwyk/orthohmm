"""Link retained native scaling observations to original factorial configurations."""

import argparse
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.collect_factorial_resources import CELLS, FIXED_INPUTS, finite, require
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.score_ygob_groups import membership, read_predictions


PANEL = ("benchmark_tools/results/threadripper_shared_panel_snapshot_20261004_v27/panel.json",
         "1a586d5d85262267ddfbdd4fccc43a2d43208dae5b028368964d710d2b587187")
PLAN = ("benchmark_tools/results/threadripper_private_commands_20260928.json",
        "a9358ac4f3c2f3eb9c1d7ce6a32528f8dd4bffef25d43cf95ac9765c5f2f3d05")
SELECTED = {7: "p1_c0_r0", 8: "p1_c1_r1", 15: "p1_c0_r0",
            16: "p1_c1_r1", 24: "p1_c1_r1", 26: "p1_c0_r0"}
COMMON = {"accuracy_profile": "high_sensitivity", "clustering": "leiden",
          "cpm_resolution": .1, "cpu_budget": 32, "evalue_threshold": .0001,
          "leiden_seed": 4, "search_kmer_k": 4, "search_max_candidates_per_query": 100,
          "search_mode": "builtin", "search_threads_per_worker": 4,
          "search_total_threads": 32, "search_workers": 8, "substitution_matrix": "BLOSUM62"}
STAGES = {"search", "edge_thresholds", "network_edges", "clustering",
          "profile_expansion", "refinement", "orthogroup_materialization"}
SCOPES = {"cpu_seconds": "native_task_subtree_cpu_stat_bracket_including_wrapper",
          "peak_memory_bytes": "native_step_lifetime_memory_peak_including_launcher",
          "wall_seconds": "native_command_monotonic_interval"}


def partition(path, format):
    if format == "space_separated_groups":
        groups = {str(i): line.split() for i, line in enumerate(path.read_text().splitlines())
                  if line.strip()}
    else:
        groups = read_predictions(path, format)
    owners = membership(groups)
    require(bool(owners), "Empty partition")
    return set(owners), {frozenset(genes) for genes in groups.values()}


def compare(old, new, reference_genes):
    require(old[0] == new[0], "Prediction universe differs")
    removed, added = old[1] - new[1], new[1] - old[1]
    changed = set().union(*removed, *added)
    return {"genes": len(new[0]), "groups": len(new[1]),
            "partition_equal": not removed and not added,
            "original_only_groups": [sorted(g) for g in sorted(removed, key=lambda g: sorted(g))],
            "scaling_only_groups": [sorted(g) for g in sorted(added, key=lambda g: sorted(g))],
            "genes_in_changed_groups": len(changed),
            "reference_genes_in_changed_groups": sorted(changed & reference_genes),
            "reference_touching_groups_equal":
                {g for g in old[1] if g & reference_genes} ==
                {g for g in new[1] if g & reference_genes}}


def validate_metrics(metrics, run, expansion):
    require(metrics["status"] == "complete", "Incomplete native metrics")
    require(all(metrics["metadata"].get(k) == v for k, v in COMMON.items()),
            "Native settings differ from frozen configuration")
    argv = run["native_argv"]
    require(metrics["command"] == [argv[0], str(Path(run["cwd"]) / "orthohmm/__main__.py"), *argv[3:]]
            and metrics["cwd"] == run["cwd"],
            "Executed command differs from native plan")
    require(argv[:3][1:] == ["-m", "orthohmm"]
            and argv[argv.index("--refinement_profile") + 1] == "default"
            and argv[argv.index("--stop") + 1] == "infer", "Wrong native execution mode")
    expected = STAGES | ({"phylogeny_candidates", "phylogeny"} if expansion else set())
    require(set(metrics["stages"]) == expected, "Missing or unexpected inference stage")
    for stage in metrics["stages"].values():
        finite(stage["wall_s"], "native stage wall", positive=True)
    require(metrics["counts"]["genes"] == 251378 and metrics["counts"]["species"] == 12,
            "Wrong native input size")
    if expansion:
        meta = metrics["metadata"]
        for key, value in {"phylogeny": "reconcile", "species_tree_mode": "infer",
                           "species_tree_rooting": "min_variance", "phylogeny_candidates": "satellite_v2",
                           "phylogeny_root_rule": "species_overlap",
                           "phylogeny_pair_rule": "positive_paralogy"}.items():
            require(meta.get(key) == value, "Wrong reconciliation setting: " + key)
        for key in ("parameters", "profile", "membership_policy", "candidate_families",
                    "seed_families", "merges", "iterations"):
            require(meta["phylogeny_candidate_profile"][key] == expansion[key],
                    "Candidate expansion differs: " + key)
        require(metrics["counts"]["phylogeny_checkpoint_hits"] == 0
                and metrics["counts"]["phylogeny_remapped_checkpoint_hits"] == 0
                and metrics["counts"]["phylogeny_species_tree_checkpoint_hit"] is False,
                "Reconciliation checkpoint reused")


def validate_attempt(row, summary, scopes):
    require(row["comparative_timing_eligible"] is True
            and row["status"] == "reviewed_shared_observation", "Unreviewed panel observation")
    for key in ("index", "job_id", "method", "proteomes", "repeat", "resources",
                "preflight_foreign_average_cores", "whole_run_maximum_foreign_average_cores"):
        require(summary[key] == row[key], "Summary differs from panel: " + key)
    require(summary["scheduler_state"] == "COMPLETED" and summary["scheduler_exit_code"] == "0:0"
            and summary["primary_resources_replayed"] is True
            and summary["shared_host_resources_reviewed"] is True
            and summary["execution_scope"] == "shared_host_matched_resources"
            and summary["resource_scopes"] == scopes == SCOPES,
            "Wrong resource admission or scope")
    for key, value in row["resources"].items():
        finite(value, key, positive=True)
    require(type(row["resources"]["peak_memory_bytes"]) is int, "Memory bytes must be integer")


def summaries(points):
    rows = []
    for cell in sorted(set(SELECTED.values())):
        selected = [p for p in points if p["cell"] == cell]
        require(len(selected) == 3 and {p["repeat"] for p in selected} == {0, 1, 2},
                "Incomplete or duplicate native repeats")
        rows.append({"dataset": "OrthoBench", "cell": cell, "repeats": 3,
                     "exact_partition_repeats": sum(p["partition"]["partition_equal"] for p in selected),
                     "reference_touching_groups_equal_all_repeats":
                         all(p["partition"]["reference_touching_groups_equal"] for p in selected),
                     "resources": {k: {"median": statistics.median(p["resources"][k] for p in selected),
                                        "minimum": min(p["resources"][k] for p in selected),
                                        "maximum": max(p["resources"][k] for p in selected)} for k in SCOPES}})
    return rows


def collect(root):
    root = Path(root).resolve()
    inputs = {}
    def checked(reference):
        observed = record(reference["path"])
        require(observed == reference, "Input checksum mismatch: " + reference["path"])
        inputs[observed["path"]] = observed
        return Path(observed["path"])
    def fixed(spec):
        ref = record(root / spec[0])
        require(ref["sha256"] == spec[1], "Fixed source checksum mismatch")
        return json.loads(checked(ref).read_text())
    ob = fixed(FIXED_INPUTS["ob_results"])
    prep = fixed(FIXED_INPUTS["ob_preparation"])
    panel, plan = fixed(PANEL), fixed(PLAN)
    require(prep["core_source_equivalence_verified"] is True, "Unverified original core")
    for ref in prep["core_sources"]:
        checked(ref)
    universe = set()
    for ref in prep["fasta_inputs"]:
        path = checked(ref)
        with path.open() as stream:
            for line in stream:
                if line.startswith(">"):
                    fields = line[1:].split()
                    require(bool(fields), "Empty FASTA identifier")
                    gene = fields[0]
                    require(gene not in universe, "Duplicate FASTA identifier")
                    universe.add(gene)
    require(len(universe) == 251378, "Wrong source gene universe")
    references = set()
    refogs = [ref for ref in ob["references"] if "low_certainty_assignments" not in ref["path"]]
    require(len(refogs) == 70, "Wrong reference inventory")
    for ref in refogs:
        references.update(checked(ref).read_text().split())
    require(references <= universe, "Reference genes absent from input universe")
    old = {}
    for cell in set(SELECTED.values()):
        path = checked(ob["predictions"][cell])
        old[cell] = partition(path, "root_hogs" if cell.endswith("r1") else "space_separated_groups")
        require(old[cell][0] == universe, "Original predictions differ from input universe")
    candidate = partition(checked(prep["candidate_arms"]["p1_c1"]["candidate_partition"]),
                          "space_separated_groups")
    original_metrics = json.loads(checked(ob["native_validation"]["p1_c1_r1"]["native_metrics"]).read_text())
    require(all(original_metrics["metadata"].get(k) == v for k, v in {
        "cpu_budget": 32, "species_tree_mode": "infer", "species_tree_rooting": "min_variance",
        "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy"}.items()),
        "Original reconciliation settings differ")
    evidence = {ref["path"]: ref for ref in panel["evidence"]}
    points = []
    for index, cell in SELECTED.items():
        row = next(r for r in panel["runs"] if r["index"] == index)
        run = plan["runs"][index]
        require(run["dataset"]["inputs"] == prep["fasta_inputs"]
                and run["proteomes"] == 12 and run["repeat"] == row["repeat"]
                and run["method"] == row["method"], "Input or run identities differ")
        path = root / ("benchmark_tools/results/threadripper_shared_attempt_%d.json" % row["job_id"])
        # The final panel pins terminal reviews, not every intermediate summary.
        summary = json.loads(checked(record(path)).read_text())
        validate_attempt(row, summary, panel["primary_scopes"])
        for category, ref in summary["reviews"].items():
            require(ref == evidence[ref["path"]], "Review differs from admitted panel")
            review = json.loads(checked(ref).read_text())
            require(review["category"] == category and review["decision"] == "passed"
                    and review["index"] == index and review["job_id"] == row["job_id"]
                    and review["plan_sha256"] == PLAN[1], "Invalid retained terminal review")
        require(set(summary["reviews"]) == {"runtime", "resources", "environment", "outputs_or_failure"},
                "Missing terminal review category")
        native_files = {ref["path"]: ref for ref in summary["native_outputs"]["checked_files"]}
        metrics = json.loads(checked(native_files[run["configuration"]["metrics"]]).read_text())
        expansion = prep["candidate_arms"]["p1_c1"]["expansion"] if cell.endswith("r1") else None
        validate_metrics(metrics, run, expansion)
        output = Path(run["configuration"]["output"])
        path = output / ("orthohmm_phylogeny/orthohmm_root_hogs.tsv" if expansion else "orthohmm_orthogroups.txt")
        current = partition(checked(native_files[str(path)]), "root_hogs" if expansion else "named_groups")
        candidate_comparison = None
        if expansion:
            native_candidate = Path(metrics["metadata"]["phylogeny_candidate_profile"]["candidate_checkpoint"])
            candidate_comparison = compare(candidate, partition(checked(record(native_candidate)),
                                           "space_separated_groups"), references)
        points.append({**row, "cell": cell, "boundary": run["boundary"],
                       "native_python": metrics["python"], "partition": compare(old[cell], current, references),
                       "candidate_partition": candidate_comparison,
                       "retained_terminal_reviews_reused": True})
    for ref in list(inputs.values()):
        checked(ref)
    return {"schema": "factorial_native_resource_linkage_v1", "status": "native_configuration_cost_linkage_checked",
            "inputs": sorted(inputs.values(), key=lambda ref: ref["path"]),
            "points": points, "summaries": summaries(points), "resource_scopes": SCOPES,
            "frozen_resources": plan["resources"], "original_core_commit": prep["core_commit"],
            "unmatched_cells": [{"dataset": dataset, "cell": cell} for dataset in ("OrthoBench", "Corrected QfO")
                                for cell in CELLS if dataset != "OrthoBench" or cell not in SELECTED.values()],
            "original_factorial_stage_costs_modified": False, "accuracy_recomputed": False,
            "native_inference_repeated": False, "publication_ready": False,
            "limitations": [
                "All six points are retained shared-host native CLI observations. Contention distortion is unknown and potentially method dependent; do not infer isolated tool speed or causal HMM/phylogeny overhead.",
                "Native launch-to-exit includes inference and output/metrics writing; excludes preparation, harness hashing, conversion and scoring. CPU includes wrapper work; lifetime cgroup memory peak includes launcher work.",
                "Frozen scientific source and original twelve-proteome bytes are checked, but the patched private runtime differs from the original factorial deployment. These are not costs of the original cached-stage executions.",
                "Phylogenetic candidate partitions differ despite identical prescribed settings and counts; two final partitions differ in nonreference groups. Retain every repeat, not only matching outputs. These costs describe the prescribed native configurations, not exact reproduction of every original factorial intermediate. Reference-touching final group equality is a structural check, not newly recomputed accuracy or proof of exact whole-output reproducibility.",
                "Native pair files, inferred trees and alignment bytes are not compared. No pairwise-orthology equivalence or cause of output differences is established.",
                "Reuse pinned terminal environment/resource/runtime reviews; do not repeat their raw forensic audits. New candidate files are checked as current readbacks, not independently attested historical write events.",
                "Fourteen other factorial configurations have no associated full-native costs in this report. Do not impute them, combine old/new memory scopes, or pool observations across deployments.",
            ]}


def render(report):
    lines = ["# Native Costs Linked To OrthoBench Configurations", "",
             "All repeats retained; shared-host observations with unknown, potentially method-dependent contention.",
             "Different private runtime from original factorial; native inference/output writing only, not whole workflow.", "",
             "| Cell | Repeats | Exact partitions | Wall median [range], s | CPU median [range], s | Lifetime peak median [range], GiB |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in report["summaries"]:
        def interval(key, divisor=1):
            v = row["resources"][key]
            return "%.3f [%.3f, %.3f]" % tuple(v[k] / divisor for k in ("median", "minimum", "maximum"))
        lines.append("| %s | 3 | %d/3 | %s | %s | %s |" % (row["cell"], row["exact_partition_repeats"],
                     interval("wall_seconds"), interval("cpu_seconds"), interval("peak_memory_bytes", 2**30)))
    lines += ["", "## Partition Readback", "", "| Run | Cell | Partition equal | Genes in changed groups | Reference genes in changed groups |",
              "| ---: | --- | --- | ---: | ---: |"]
    for point in report["points"]:
        p = point["partition"]
        lines.append("| %d | %s | %s | %d | %d |" % (point["index"], point["cell"], p["partition_equal"],
                     p["genes_in_changed_groups"], len(p["reference_genes_in_changed_groups"])))
    lines += ["", "## Limits", "", *["- " + item for item in report["limitations"]]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output_directory.exists(), "Refusing existing output directory")
    report = collect(args.root)
    report["source"] = record(Path(__file__))
    args.output_directory.mkdir(parents=True)
    (args.output_directory / "linkage.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (args.output_directory / "linkage.md").write_text(render(report))
    print(json.dumps({"status": report["status"], "points": len(report["points"]),
                      "configurations": len(report["summaries"]), "unmatched_cells": len(report["unmatched_cells"])}))


if __name__ == "__main__":
    main()
