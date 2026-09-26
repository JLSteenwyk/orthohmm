"""Strict partition/score readback for the installed OrthoBench reproduction."""

import argparse
import csv
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.verify_ygob_validation import require_completed_job

BASELINE = "benchmark_tools/results/orthobench_factorial_results_20260916.json"
BASELINE_SHA = "6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141"
PLAN_SHA = "5fd8dc70c337951706fe9246c2bf08d1191196939da066c5b74bba0bd3dd2b12"


def read_root_hogs(path, universe):
    groups, seen, labels = [], set(), set()
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["root_hog", "source_family", "genes"]:
            raise ValueError("Unexpected root-HOG columns")
        for row in reader:
            members = row["genes"].split(",")
            if (not row["root_hog"] or not row["source_family"] or row["root_hog"] in labels
                    or not all(members) or len(set(members)) != len(members)
                    or seen.intersection(members) or not set(members) <= universe):
                raise ValueError("Invalid, duplicate or foreign group membership")
            seen.update(members)
            labels.add(row["root_hog"])
            groups.append(frozenset(members))
    if seen != universe:
        raise ValueError("Partition does not cover every input gene exactly once")
    return groups


def compare_partitions(original, current):
    a, b = set(original), set(current)
    return dict(original_groups=len(a), current_groups=len(b), identical_groups=len(a & b),
                original_only_groups=len(a - b), current_only_groups=len(b - a),
                genes_in_changed_groups=len(set().union(*(a - b))) if a - b else 0,
                label_invariant_equal=a == b)


def baseline_inputs(repo):
    source = record(repo / BASELINE)
    if source["sha256"] != BASELINE_SHA:
        raise ValueError("Changed baseline score evidence")
    baseline = json.loads((repo / BASELINE).read_text())
    references, uncertain = {}, {}
    checked = [source, baseline["predictions"]["p1_c1_r1"], *baseline["references"]]
    for item in checked:
        check(item)
    for item in baseline["references"]:
        path = Path(item["path"])
        target = uncertain if path.parent.name == "low_certainty_assignments" else references
        target[path.name] = set(path.read_text().splitlines())
    if len(references) != 70:
        raise ValueError("Expected all 70 RefOGs")
    return baseline, references, uncertain, checked


def verify_baseline(repo, input_dir):
    baseline, references, uncertain, checked = baseline_inputs(repo)
    inputs = sorted(input_dir.glob("*.fa"))
    universe = fasta_ids(inputs)
    if len(inputs) != 12 or len(universe) != 251378:
        raise ValueError("Wrong full input universe")
    groups = read_root_hogs(Path(baseline["predictions"]["p1_c1_r1"]["path"]), universe)
    score = score_partition(groups, references, uncertain)
    if score != baseline["scores"]["p1_c1_r1"]:
        raise ValueError("Strict baseline readback differs from retained score")
    return dict(status="historical_partition_strict_readback_verified", checked_records=checked,
                fasta_records=[record(p) for p in inputs], genes=len(universe), groups=len(groups), score=score,
                source=record(__file__), scorer=record(Path(__file__).with_name("score_orthobench_partition.py")),
                installed_run_admitted=False)


def audit(repo, directory, job):
    plan_path = directory / "plan.json"
    plan_record = record(plan_path)
    if plan_record["sha256"] != PLAN_SHA:
        raise ValueError("Changed installed-run plan")
    plan = json.loads(plan_path.read_text())
    execution_path = directory / "execution.json"
    execution = json.loads(execution_path.read_text())
    if (execution["status"] != "native_completed_pending_independent_scientific_readback"
            or execution["returncode"] != 0 or execution["job_id"] != str(job)
            or execution["plan"] != plan_record):
        raise ValueError("Native execution has not completed successfully")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, job)
    for item in plan["checked_records"]:
        check(item)
    baseline, refs, uncertain, checked = baseline_inputs(repo)
    universe = fasta_ids(sorted((directory / "input").glob("*.fa")))
    previous = read_root_hogs(Path(plan["baseline_partition"]["path"]), universe)
    phylo = directory / "inference/orthohmm_phylogeny"
    groups = read_root_hogs(phylo / "orthohmm_root_hogs.tsv", universe)
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    manifest = json.loads((phylo / "provenance_manifest.json").read_text())
    if (summary["checkpoint_hits"] != 0 or summary["remapped_checkpoint_hits"] != 0
            or summary["species_tree_checkpoint_hit"] is not False or summary["species_tree_families"] <= 0
            or summary["reconciled_families"] <= 0 or summary["root_hogs"] != len(groups)
            or manifest["species_tree_mode"] != "infer" or len(manifest["species_tree_taxa"]) != 12):
        raise ValueError("Expected fresh full inferred phylogeny")
    trees = sorted((phylo / "gene_trees").glob("*.reconciled.nwk"))
    if len(trees) != summary["reconciled_families"] or any(p.stat().st_size == 0 for p in trees):
        raise ValueError("Missing reconciled tree files")
    score = score_partition(groups, refs, uncertain)
    old = baseline["scores"]["p1_c1_r1"]
    partition_comparison = compare_partitions(previous, groups)
    result = dict(status="installed_partition_and_scores_read_back", plan=plan_record,
        execution=record(execution_path), scheduler=scheduler, accounting=accounting,
        genes=len(universe), score=score, baseline_score=old,
        score_differences={k:score[k]-old[k] for k in ("f_score", "precision", "recall")},
        family_records_equal=score["refog_records"] == old["refog_records"],
        partition_comparison=partition_comparison, summary=summary,
        checked_records=[*plan["checked_records"], *checked],
        outputs=[record(phylo / name) for name in ("orthohmm_root_hogs.tsv", "reconciliation_summary.json",
                                                  "provenance_manifest.json", "species_tree.rooted.nwk")],
        source=record(__file__), scorer=record(Path(__file__).with_name("score_orthobench_partition.py")),
        full_phylogeny_validation_complete=False, historical_scores_replaced=False, publication_ready=False,
        limitations=["Strict coverage and score readback; gene-tree contents and pair table require separate validation.",
                     "No score or partition difference is silently promoted or treated as failure-free equivalence.",
                     "Shared-host installation reproduction, not independent validation or controlled timing."])
    for item in result["checked_records"]:
        check(item)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    if args.baseline_only:
        result = verify_baseline(args.repo.resolve(), args.directory.resolve() / "input")
    elif args.job:
        result = audit(args.repo.resolve(), args.directory.resolve(), args.job)
    else:
        parser.error("Require --job or --baseline-only")
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
