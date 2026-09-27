"""Audit fresh canonical phylogeny and compare both retained OrthoBench runs."""

import argparse
import json
from pathlib import Path

from benchmark_tools import admit_canonical_ob_phylogeny as admission
from benchmark_tools import audit_phylogeny_structure as structure
from benchmark_tools import audit_phylogeny_sequences as sequences
from benchmark_tools import audit_phylogeny_events as events
from benchmark_tools import audit_phylogeny_hierarchy as hierarchy
from benchmark_tools.audit_installed_orthobench import (
    baseline_inputs, compare_partitions, read_root_hogs, verify_input_inventory,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.score_orthobench_partition import score_partition

SCORER_SHA = "0e24c302a9d12ce8c82a4111ff0321c87303469b3c866c317dc2629ad66f9746"
FRESH_RECEIPT_SHA = "79fc6ce6a169d486235c811c7d6eec9019127958c07e56d8ff982675cd4149e0"


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def comparison(original, current, original_score, current_score):
    old = {r["refog"]: r for r in original_score["refog_records"]}
    new = {r["refog"]: r for r in current_score["refog_records"]}
    if len(old) != len(original_score["refog_records"]) or len(new) != len(current_score["refog_records"]) or old.keys() != new.keys():
        raise ValueError("Reference-family inventories differ or contain duplicates")
    return dict(partitions=compare_partitions(original, current),
        score_differences_percentage_points={k:current_score[k]-original_score[k]
            for k in ("f_score", "precision", "recall")},
        changed_reference_families=[dict(refog=k, original=old[k], current=new[k])
            for k in sorted(old) if old[k] != new[k]],
        score_objects_equal=original_score == current_score)


def readback(repo, directory, job, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = admission.audit(directory, job)
    plan = json.loads((directory / "plan.json").read_text())
    scorer = record(repo / "benchmark_tools/score_orthobench_partition.py")
    if scorer["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    baseline, refs, uncertain, baseline_records = baseline_inputs(repo)
    fresh_path = repo / "benchmark_tools/results/installed_orthobench_readback_20260926.json"
    fresh_record = record(fresh_path)
    if fresh_record["sha256"] != FRESH_RECEIPT_SHA:
        raise ValueError("Changed fresh-installed evidence")
    fresh = json.loads(fresh_path.read_text())
    for item in fresh["readback_records"]:
        check(item)
    fresh_scores_record = next(r for r in fresh["readback_records"] if r["path"].endswith("/readback_scores.json"))
    fresh_scores = json.loads(Path(fresh_scores_record["path"]).read_text())
    previous_output = next(r for r in fresh_scores["outputs"] if r["path"].endswith("/orthohmm_root_hogs.tsv"))
    check(previous_output)
    output.mkdir(parents=True)
    save(output / "admission.json", admitted)
    phylo = Path(plan["inference"]) / "orthohmm_phylogeny"
    reports = {"structure": structure.audit(phylo, Path(plan["inputs"]))}
    save(output / "structure.json", reports["structure"])
    reports["sequences"] = sequences.audit(phylo, output / "structure.json")
    save(output / "sequences.json", reports["sequences"])
    reports["events"] = events.audit(phylo, output / "structure.json", Path(plan["constraints"]["path"]))
    save(output / "events.json", reports["events"])
    reports["hierarchy"] = hierarchy.audit(phylo, output / "events.json")
    save(output / "hierarchy.json", reports["hierarchy"])
    universe, _ = verify_input_inventory(Path(plan["inputs"]), dict(
        checked_records=plan["checked_records"], expected_species=12, expected_genes=251378))
    paths = dict(historical=Path(baseline["predictions"]["p1_c1_r1"]["path"]),
                 fresh_installed=Path(previous_output["path"]), canonical=phylo / "orthohmm_root_hogs.tsv")
    groups = {name:read_root_hogs(path, universe) for name, path in paths.items()}
    scores = {name:score_partition(partition, refs, uncertain) for name, partition in groups.items()}
    if scores["historical"] != baseline["scores"]["p1_c1_r1"] or scores["fresh_installed"] != fresh_scores["score"]:
        raise ValueError("Retained scores did not reproduce exactly")
    checked = [*baseline_records, fresh_record, *fresh["readback_records"], previous_output,
               *plan["checked_records"], *[record(p) for p in paths.values()]]
    for report in reports.values():
        checked.extend(report["checked_records"])
    for item in checked:
        check(item)
    if admission.audit(directory, job) != admitted:
        raise ValueError("Native admission changed during readback")
    result = dict(status="canonical_phylogeny_scientific_readback_complete", job_id=job,
        plan=admitted["plan"], genes=len(universe), scores=scores,
        comparisons={name:comparison(groups[name], groups["canonical"], scores[name], scores["canonical"])
            for name in ("historical", "fresh_installed")},
        partitions={name:record(path) for name, path in paths.items()},
        readbacks=[record(output / (name + ".json")) for name in ("admission", *reports)],
        sources=[record(Path(module.__file__)) for module in (admission, structure, sequences, events, hierarchy)],
        scorer=scorer, source=record(__file__), summary=admitted["summary"],
        historical_scores_replaced=False, production_defaults_changed=False, publication_ready=False,
        limitations=["Development-exposed partial-stage reproducibility experiment, not independent validation.",
                     "Conditional rule/content audits do not establish biological truth or optimal trees.",
                     "Shared-host resource measurements cannot support controlled speed comparisons."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    readback(args.repo.resolve(), args.directory.resolve(), args.job, args.output.resolve())
