"""Admit fresh QfO phylogeny and compare all retained native predictions."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import subprocess

from benchmark_tools import audit_phylogeny_structure as structure
from benchmark_tools import audit_phylogeny_sequences as sequences
from benchmark_tools import audit_phylogeny_events as events
from benchmark_tools import audit_phylogeny_hierarchy as hierarchy
from benchmark_tools.audit_installed_orthobench import read_root_hogs, compare_partitions
from benchmark_tools.run_qfo_fresh_phylogeny import config_arguments, record, save, validate
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.qfo_input_inventory import verify_qfo_inputs

JOB = 22329
PLAN_SHA = "57247bf76d26c1ff54cf13aa4c9d0a9622f0b8a2902c45ef222842513888d030"
SUBMISSION_SHA = "3d54dc754b5d4c467f1c2fa049cf8969c0be9b209290c76c992f6be24a75d5b0"


def verify_summary(summary, native_summary):
    if summary != dict(native_summary, schema_version=1):
        raise ValueError("Completion summary differs")


def verify_native(plan, identity, started, complete, started_record):
    pinned = {r["path"]: r for r in plan["checked_records"]}
    source = str(Path(plan["repo"]) / "benchmark_tools/run_qfo_fresh_phylogeny.py")
    if (started["plan"] != identity or started["source"] != pinned[source]
            or started["job_id"] != str(JOB) or started["executable"] != plan["python"]
            or started["environment"] != plan["environment"]
            or started["config"] != config_arguments(plan) or started["attempts"] != 1
            or started["checkpoint_reuse"] is not False
            or complete["status"] != "fresh_phylogeny_complete_pending_readback"
            or complete["plan"] != identity or complete["started"] != started_record
            or complete["accuracy_evaluated"] is not False):
        raise ValueError("Native provenance or scope differs")
    summary = complete["summary"]
    if (summary["checkpoint_hits"] != 0 or summary["remapped_checkpoint_hits"] != 0
            or summary["species_tree_checkpoint_hit"] is not False
            or summary["reconciled_families"] < 1 or summary["species_tree_families"] < 1):
        raise ValueError("Fresh trees not established")


def admit(repo, directory):
    accounting = subprocess.check_output(["sacct", "-j", str(JOB), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, JOB)
    if scheduler["AllocCPUS"] != "32":
        raise ValueError("Wrong CPU allocation")
    path = directory / "plan.json"
    plan = validate(path, PLAN_SHA)
    if plan["repo"] != str(repo) or plan["directory"] != str(directory):
        raise ValueError("Wrong run locations")
    frozen = repo / "benchmark_tools/results/qfo_fresh_phylogeny_submission_22329.json"
    if record(frozen)["sha256"] != SUBMISSION_SHA:
        raise ValueError("Changed submission")
    submission = json.loads((directory / "submission.json").read_text())
    if (submission != json.loads(frozen.read_text()) or submission["job_id"] != str(JOB)
            or submission["plan"] != record(path)):
        raise ValueError("Wrong submitted attempt")
    native = directory / "native"
    if (native / "failure.json").exists():
        raise ValueError("Retained failure cannot be admitted")
    started = json.loads((native / "started.json").read_text())
    complete = json.loads((native / "complete.json").read_text())
    verify_native(plan, record(path), started, complete, record(native / "started.json"))
    outputs = [record(p) for p in sorted((native / "inference").rglob("*")) if p.is_file()]
    if complete["outputs"] != outputs:
        raise ValueError("Changed or missing native outputs")
    summary = json.loads((native / "inference/orthohmm_phylogeny/reconciliation_summary.json").read_text())
    verify_summary(summary, complete["summary"])
    inputs = verify_qfo_inputs(plan)
    return dict(plan=record(path), scheduler=scheduler, outputs=outputs, inputs=inputs,
        records=[record(p) for p in (frozen, directory / "submission.json", native / "started.json",
            native / "complete.json", directory / "time.txt", directory / f"native-{JOB}.log")])


def pair_rows(path):
    previous = None
    for row in structure.rows(path, ["gene_a", "species_a", "gene_b", "species_b", "confidence"]):
        key = (row["gene_a"], row["gene_b"])
        if not key[0] < key[1] or (previous is not None and key <= previous):
            raise ValueError("Pairs must be unique and sorted")
        if row["confidence"] not in {"high", "medium", "low"}:
            raise ValueError("Invalid pair confidence")
        previous = key
        yield key, (row["species_a"], row["species_b"], row["confidence"])


def compare_pairs(left, right):
    left, right = iter(left), iter(right)
    a, b = next(left, None), next(right, None)
    counts = Counter(shared=0, left_only=0, right_only=0, annotations_changed=0)
    while a is not None or b is not None:
        if b is None or (a is not None and a[0] < b[0]):
            counts["left_only"] += 1
            a = next(left, None)
        elif a is None or b[0] < a[0]:
            counts["right_only"] += 1
            b = next(right, None)
        else:
            counts["shared"] += 1
            counts["annotations_changed"] += a[1] != b[1]
            a, b = next(left, None), next(right, None)
    return dict(counts, pair_sets_equal=not(counts["left_only"] or counts["right_only"]),
                confidence_tables_equal=not(counts["left_only"] or counts["right_only"] or counts["annotations_changed"]))


def family_partitions(path):
    result = defaultdict(set)
    for row in structure.rows(path, ["root_hog", "source_family", "genes"]):
        result[row["source_family"]].add(frozenset(row["genes"].split(",")))
    return result


def scientific_readback(phylo, inputs, constraints, output):
    output.mkdir(exist_ok=False)
    summary = structure.audit(phylo, inputs)
    save(output / "structure.json", summary)
    save(output / "sequences.json", sequences.audit(phylo, output / "structure.json"))
    save(output / "events.json", events.audit(phylo, output / "structure.json", constraints))
    save(output / "hierarchy.json", hierarchy.audit(phylo, output / "events.json"))
    return summary


def readback(repo, directory, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = admit(repo, directory)
    plan = json.loads((directory / "plan.json").read_text())
    output.mkdir(parents=True)
    save(output / "admission.json", admitted)
    constraints = Path(plan["constraints"]["path"])
    paths = dict(historical=Path(plan["historical"]) / "orthohmm_phylogeny",
                 fresh=directory / "native/inference/orthohmm_phylogeny")
    checked = []
    for label, phylo in paths.items():
        report = scientific_readback(phylo, Path(plan["input_directory"]), constraints, output / label)
        if report["genes"] != plan["expected_genes"] or report["species"] != plan["expected_species"]:
            raise ValueError("Wrong complete input universe")
        checked.extend(report["checked_records"])
    universe = structure.fasta_ids([Path(r["path"]) for r in plan["fastas"]])
    root_paths = {label: p / "orthohmm_root_hogs.tsv" for label, p in paths.items()}
    roots = {label: read_root_hogs(p, universe) for label, p in root_paths.items()}
    families = {label: family_partitions(p) for label, p in root_paths.items()}
    changed_families = [name for name in sorted(families["historical"].keys() | families["fresh"].keys())
                       if families["historical"].get(name) != families["fresh"].get(name)]
    pairs = compare_pairs(*(pair_rows(p / "orthohmm_pairwise_orthologs_confidence.tsv") for p in paths.values()))
    manifests = {label: json.loads((p / "provenance_manifest.json").read_text()) for label, p in paths.items()}
    summaries = {label: json.loads((p / "reconciliation_summary.json").read_text()) for label, p in paths.items()}
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed scientific artifact")
    if admit(repo, directory) != admitted:
        raise ValueError("Changed native admission during readback")
    result = dict(status="fresh_qfo_phylogeny_scientific_readback_complete", job_id=JOB,
        plan=record(directory / "plan.json"), admission=record(output / "admission.json"),
        partition=compare_partitions(roots["historical"], roots["fresh"]), changed_families=changed_families,
        pairs=pairs, summaries=summaries,
        species_tree_bytes_equal=manifests["historical"]["species_tree_sha256"] == manifests["fresh"]["species_tree_sha256"],
        membership={k: m["membership_reconciliation"] for k, m in manifests.items()},
        reports=[record(p) for p in sorted(output.glob("*/*.json"))], source=record(__file__),
        accuracy_evaluated=False, historical_scores_replaced=False, publication_ready=False,
        limitations=["Downstream cached-candidate comparison, not fresh search or independent generalization",
                     "Tree byte differences, if present, do not by themselves establish topology differences"])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(readback(args.repo.resolve(), args.directory.resolve(), args.output.resolve()), indent=2))
