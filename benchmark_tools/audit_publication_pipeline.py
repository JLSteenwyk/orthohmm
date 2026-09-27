"""Read back the explicit canonical full-pipeline experiment without scoring."""

import argparse
import json
from pathlib import Path

from benchmark_tools import audit_phylogeny_structure as structure
from benchmark_tools import audit_phylogeny_sequences as sequences
from benchmark_tools import audit_phylogeny_events as events
from benchmark_tools import audit_phylogeny_hierarchy as hierarchy
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def constraint_argument(path):
    value = json.loads(path.read_text())
    if not isinstance(value, list):
        raise ValueError("Constraint trace must be a list")
    return path if value else None


def audit(directory, output):
    if output.exists():
        raise FileExistsError(output)
    if (directory / "failure.json").exists():
        raise ValueError("Failed attempt cannot be admitted")
    started = json.loads((directory / "started.json").read_text())
    complete = json.loads((directory / "complete.json").read_text())
    if (complete["status"] != "native_complete_pending_scientific_readback"
            or complete["started"] != record(directory / "started.json")
            or started["attempts"] != 1 or started["checkpoint_reuse"] is not False
            or started["production_default_changed"] is not False
            or len(complete["policy_application"]) != 1
            or complete["policy_application"][0]["policy"] != "canonical_directed_pair_v1"
            or complete["policy_application"][0]["scores_modified"] is not False):
        raise ValueError("Wrong execution scope or policy")
    records = [started["source"], started["policy"], *started["inputs"],
               *started["scientific_sources"], *started["tools"]]
    for item in records:
        check(item)
    inputs = Path(started["arguments"][0])
    if [record(p) for p in sorted(inputs.iterdir()) if p.is_file()] != started["inputs"]:
        raise ValueError("Changed input inventory")
    phylo = directory / "inference/orthohmm_phylogeny"
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    if (summary["checkpoint_hits"] != 0 or summary["remapped_checkpoint_hits"] != 0
            or summary["species_tree_checkpoint_hit"] is not False
            or summary["species_tree_families"] < 1 or summary["reconciled_families"] < 1):
        raise ValueError("No fresh inferred species/gene phylogeny evidence")
    output.mkdir(parents=True)
    save(output / "structure.json", structure.audit(phylo, inputs))
    save(output / "sequences.json", sequences.audit(phylo, output / "structure.json"))
    constraints = directory / "inference/orthohmm_working_res/phylogeny_candidate_merges.json"
    records.append(record(constraints))
    save(output / "events.json", events.audit(phylo, output / "structure.json", constraint_argument(constraints)))
    save(output / "hierarchy.json", hierarchy.audit(phylo, output / "events.json"))
    for item in records:
        check(item)
    result = dict(status="canonical_full_pipeline_scientific_readback_verified",
        started=record(directory / "started.json"), complete=record(directory / "complete.json"),
        source=record(__file__), summary=summary,
        constraints=record(constraints),
        reports={name: record(output / (name + ".json")) for name in
                 ("structure", "sequences", "events", "hierarchy")},
        accuracy_evaluated=False, publication_ready=False,
        limitations=["Same-host execution; no cross-platform reproducibility claim.",
                     "Canonical ordering remains an explicit experimental policy."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.directory.resolve(), args.output.resolve()), indent=2))
