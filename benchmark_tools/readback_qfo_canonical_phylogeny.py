"""Validate the canonical QfO arm and all three frozen native contrasts."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.compare_qfo_phylogeny_arms import compare
from benchmark_tools.readback_qfo_fresh_phylogeny import scientific_readback, verify_summary
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.run_qfo_canonical_phylogeny import config_arguments, record, save, validate
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.qfo_input_inventory import verify_qfo_inputs


def verify_execution(plan, identity, job, started, started_record, complete, cache, cache_record, copies_record):
    pinned = {r["path"]: r for r in plan["checked_records"]}
    source = str(Path(plan["repo"]) / "benchmark_tools/run_qfo_canonical_phylogeny.py")
    if (started["plan"] != identity or started["job_id"] != str(job)
            or started["source"] != pinned[source] or started["executable"] != plan["python"]
            or started["environment"] != plan["environment"] or started["config"] != config_arguments(plan)
            or started["reuse_admission"] != plan["reuse_admission"] or started["attempts"] != 1
            or complete["status"] != "canonical_phylogeny_complete_pending_readback"
            or complete["started"] != started_record or complete["plan"] != identity
            or complete["cache"] != cache_record or complete["copies"] != copies_record
            or complete["accuracy_evaluated"] is not False):
        raise ValueError("Wrong canonical provenance")
    summary = complete["summary"]
    if (summary["checkpoint_hits"] != len(cache["included_families"])
            or summary["remapped_checkpoint_hits"] != 0 or summary["species_tree_checkpoint_hit"] is not False
            or cache["species_tree_cache_copied"] is not False or cache["reconciliation_outputs_copied"] is not False
            or summary["species_tree_families"] < 1):
        raise ValueError("Wrong canonical reuse scope")


def verify_cache(plan, gate, cache, copies, target):
    from benchmark_tools.qfo_phylogeny_cache import cache_paths
    prefix = Path(plan["aligner"]).parent.parent
    tools = [r for r in plan["checked_records"] if Path(r["path"]).is_relative_to(prefix)
             or r["path"] == plan["tree_builder"]]
    if (cache["source_directory"] != gate["source_directory"] or cache["tool_records"] != plan["tool_records"]
            or sorted(tools, key=lambda r: r["path"]) != sorted(plan["tool_records"], key=lambda r: r["path"])):
        raise ValueError("Cache source or external tool identities differ")
    included = cache["included_families"]
    if len(included) != len(set(included)):
        raise ValueError("Duplicate cache families")
    expected = {relative for family in included for relative in cache_paths(family)}
    relative = [r["relative"] for r in cache["files"]]
    if len(relative) != len(set(relative)) or set(relative) != expected:
        raise ValueError("Unexpected cache inventory")
    admitted = {r["path"]: r for r in gate["admitted_outputs"]}
    expected_copies = []
    for row in cache["files"]:
        source = row["source"]
        if (source["path"] != str(Path(gate["source_directory"]) / row["relative"])
                or admitted.get(source["path"]) != source or record(source["path"]) != source):
            raise ValueError("Unadmitted or changed cache source")
        expected_copies.append(dict(source, path=str(target / row["relative"])))
    # Copies record the pre-execution checkpoint bytes; reconciliation rewrites them.
    if copies != expected_copies:
        raise ValueError("Cache copy receipt differs")


def admit(repo, directory, job, plan_sha, submission_sha):
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    if scheduler["AllocCPUS"] != "32":
        raise ValueError("Wrong canonical CPU allocation")
    plan_path = directory / "plan.json"
    plan = validate(plan_path, plan_sha)
    if plan["repo"] != str(repo) or plan["directory"] != str(directory):
        raise ValueError("Wrong canonical locations")
    submission_path = directory / "submission.json"
    if record(submission_path)["sha256"] != submission_sha:
        raise ValueError("Changed frozen canonical submission")
    submission = json.loads(submission_path.read_text())
    if submission["job_id"] != str(job) or submission["plan"] != record(plan_path):
        raise ValueError("Wrong canonical submission binding")
    native = directory / "native"
    if (native / "failure.json").exists():
        raise ValueError("Cannot admit failed canonical execution")
    receipt_names = ("started.json", "complete.json", "cache_manifest.json", "cache_copies.json")
    started, complete, cache, copies = [json.loads((native / n).read_text()) for n in receipt_names]
    verify_execution(plan, record(plan_path), job, started, record(native / "started.json"), complete,
                     cache, record(native / "cache_manifest.json"), record(native / "cache_copies.json"))
    gate = json.loads(Path(plan["reuse_admission"]["path"]).read_text())
    if gate["status"] != "fresh_qfo_raw_tree_reuse_admitted" or gate["native_job"] != 22329 or gate["readback_job"] != 22330:
        raise ValueError("Wrong retained-arm reuse admission")
    verify_cache(plan, gate, cache, copies, native / "inference/orthohmm_phylogeny")
    outputs = [record(p) for p in sorted((native / "inference").rglob("*")) if p.is_file()]
    if outputs != complete["outputs"]:
        raise ValueError("Changed canonical output inventory")
    verify_summary(json.loads((native / "inference/orthohmm_phylogeny/reconciliation_summary.json").read_text()),
                   complete["summary"])
    inputs = verify_qfo_inputs(plan)
    return dict(plan=record(plan_path), scheduler=scheduler, outputs=outputs, inputs=inputs,
        records=[record(native / n) for n in receipt_names] + [record(submission_path),
            record(directory / "time.txt"), record(directory / f"native-{job}.log")])


def readback(repo, directory, job, plan_sha, submission_sha, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = admit(repo, directory, job, plan_sha, submission_sha)
    plan = json.loads((directory / "plan.json").read_text())
    gate = json.loads(Path(plan["reuse_admission"]["path"]).read_text())
    output.mkdir(parents=True)
    save(output / "admission.json", admitted)
    canonical = directory / "native/inference/orthohmm_phylogeny"
    report = scientific_readback(canonical, Path(plan["input_directory"]),
        Path(plan["constraints"]["path"]), output / "canonical")
    if report["genes"] != plan["expected_genes"] or report["species"] != plan["expected_species"]:
        raise ValueError("Wrong canonical universe")
    paths = dict(historical=Path(plan["historical"]) / "orthohmm_phylogeny",
                 fresh_retained=Path(gate["source_directory"]), canonical=canonical)
    contrasts = compare(paths, fasta_ids([Path(r["path"]) for r in plan["fastas"]]))
    if admit(repo, directory, job, plan_sha, submission_sha) != admitted:
        raise ValueError("Canonical admission changed during scientific readback")
    result = dict(status="canonical_qfo_three_way_readback_complete", job_id=job,
        plan=record(directory / "plan.json"), admission=record(output / "admission.json"),
        reports=[record(p) for p in sorted((output / "canonical").glob("*.json"))],
        comparison=contrasts, source=record(__file__), accuracy_evaluated=False,
        historical_scores_replaced=False, publication_ready=False,
        limitation="Downstream cached-search experiment, not fresh full-pipeline validation or independent accuracy")
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--plan-sha", required=True)
    parser.add_argument("--submission-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(readback(args.repo.resolve(), args.directory.resolve(), args.job,
        args.plan_sha, args.submission_sha, args.output.resolve()), indent=2))
