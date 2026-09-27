"""Bind fresh canonical phylogeny artifacts before independent scientific audits."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_installed_orthobench import verify_input_inventory
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.verify_ygob_validation import require_completed_job

FULL_PLAN_SHA = "2ded17c15cf3510e02f1683973edf211dc8b7f7f86037dcf83e121a948dfafdc"


def verify_artifacts(directory, expected_sha, species, genes):
    def read(name):
        return json.loads((directory / name).read_text())

    plan_record = record(directory / "plan.json")
    if plan_record["sha256"] != expected_sha:
        raise ValueError("Changed plan")
    plan = read("plan.json")
    if (Path(plan["output"]) != directory or plan["attempts"] != 1
            or plan["checkpoint_reuse"] is not False or plan["scoring"] is not False
            or (directory / "failure.json").exists()):
        raise ValueError("Wrong scope or failed native attempt")
    for item in [plan["source"], *plan["checked_records"]]:
        check(item)
    universe, inputs = verify_input_inventory(Path(plan["inputs"]), dict(
        checked_records=plan["checked_records"], expected_species=species, expected_genes=genes))
    started, completed, replay = (read(n) for n in
        ("started.json", "native_complete.json", "replay.json"))
    if (started["plan"] != plan_record or started["runtime"] != plan["runtime"]
            or completed["plan"] != plan_record
            or completed["status"] != "fresh_canonical_phylogeny_complete_pending_readback"
            or completed["accuracy_evaluated"] is not False
            or completed["replay"] != record(directory / "replay.json")
            or replay["status"] != "complete"):
        raise ValueError("Native completion/runtime binding differs")
    expected_command = [plan["source"]["path"], "--run", plan_record["path"],
                        "--plan-sha256", expected_sha]
    if started["command"] != expected_command or replay["command"] != [plan["python"], *expected_command]:
        raise ValueError("Native command differs")
    expected_args = ["--fasta-directory", plan["inputs"], "--candidate-clusters", plan["candidate"]["path"],
        "--membership-constraints", plan["constraints"]["path"], "--output-directory", plan["inference"],
        "--json", str(directory / "replay.json"), "--cpu", str(plan["cpu"]),
        "--species-tree-mode", "infer", "--species-tree-rooting", "min_variance",
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
        "--aligner", plan["aligner"], "--tree-builder", plan["tree_builder"]]
    if started["arguments"] != expected_args:
        raise ValueError("Replay arguments differ")
    actual = replay["input"]
    if (actual["candidate_clusters"] != plan["candidate"]
            or actual["membership_constraints"] != plan["constraints"]
            or actual["fasta_directory"] != plan["inputs"]
            or sorted(actual["files"]) != sorted(Path(r["path"]).name for r in inputs)):
        raise ValueError("Replay input binding differs")
    phylo = Path(plan["inference"]) / "orthohmm_phylogeny"
    outputs = replay["outputs"]
    for key, name in (("manifest", "provenance_manifest.json"),
                      ("root_hogs", "orthohmm_root_hogs.tsv"),
                      ("summary", "reconciliation_summary.json")):
        if outputs[key] != record(phylo / name):
            raise ValueError("Replay output binding differs")
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    manifest = json.loads((phylo / "provenance_manifest.json").read_text())
    for key, value in summary.items():
        if key != "schema_version" and replay["counts"].get(key) != value:
            raise ValueError("Replay counts differ")
    if (summary["checkpoint_hits"] != 0 or summary["remapped_checkpoint_hits"] != 0
            or summary["species_tree_checkpoint_hit"] is not False
            or summary["species_tree_families"] <= 0 or summary["reconciled_families"] <= 0
            or manifest["input_cluster_sha256"] != plan["candidate"]["sha256"]
            or manifest["results"] != summary or len(manifest["species_tree_taxa"]) != species):
        raise ValueError("Not a fresh inferred phylogeny")
    fixed = dict(species_tree_mode="infer", species_tree_rooting="min_variance",
                 root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    if any(manifest.get(k) != v or replay["parameters"].get(k) != v for k, v in fixed.items()):
        raise ValueError("Frozen phylogeny rules differ")
    if (replay["parameters"]["checkpoint_source"] is not None
            or replay["parameters"]["species_tree"] is not None
            or replay["parameters"]["explicit_unconstrained_ablation"] is not False
            or manifest["membership_reconciliation"]["policy"] != "high_confidence_pair"):
        raise ValueError("Unexpected reuse or unconstrained inference")
    return dict(status="native_artifacts_bound_pending_scientific_readback", plan=plan_record,
        genes=len(universe), species=species, input_inventory=inputs, summary=summary,
        evidence=[record(directory / n) for n in ("started.json", "native_complete.json", "replay.json")],
        outputs=outputs, source=record(__file__), scientific_scores_admitted=False)


def audit(directory, job):
    submission = json.loads((directory / "submission.json").read_text())
    if submission["job_id"] != str(job) or submission["plan"] != record(directory / "plan.json"):
        raise ValueError("Scheduler submission differs")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, job)
    result = verify_artifacts(directory, FULL_PLAN_SHA, 12, 251378)
    result.update(scheduler=scheduler, accounting=accounting, submission=record(directory / "submission.json"))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.directory.resolve(), args.job)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
