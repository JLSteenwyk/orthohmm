"""Gate canonical QfO cache use on the completed native and scientific audits."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.readback_qfo_fresh_phylogeny import admit as admit_native
from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.verify_ygob_validation import require_completed_job

READBACK_JOB = 22332
READBACK_PLAN_SHA = "d42e759f2d42442f250a440a1138be30660b2d842c2477bb4d44024e54edcdca"
SUBMISSION_SHA = "d06a1671d93418ef7ee8620bcd44b478d6b19aef3ea1efdcd4be811c0843d24d"
ARCHIVED_ORCHESTRATION = (
    "admit_qfo_raw_tree_reuse.py", "run_qfo_canonical_phylogeny.py",
    "readback_qfo_canonical_phylogeny.py",
)
REPORTS = dict(structure="phylogeny_structure_verified", sequences="phylogeny_sequence_content_verified",
               events="frozen_phylogeny_event_pair_semantics_verified",
               hierarchy="phylogeny_selection_and_hierarchy_verified")


def audited_records(repo, directory, records):
    """Preserve audited orchestration bytes while binding new orchestration separately."""
    replacements = {str(repo / "benchmark_tools" / name):
                    directory / "readback_v2_source_archive" / name
                    for name in ARCHIVED_ORCHESTRATION}
    if not replacements.keys() <= {r["path"] for r in records}:
        raise ValueError("Missing audited orchestration identity")
    verified, archived = [], []
    for item in records:
        path = replacements.get(item["path"], Path(item["path"]))
        actual = record(path)
        if actual["bytes"] != item["bytes"] or actual["sha256"] != item["sha256"]:
            raise ValueError("Changed audited dependency: " + item["path"])
        verified.append(actual)
        if item["path"] in replacements:
            archived.append(dict(original=item, archived=actual))
    return verified, archived


def collect_records(value):
    """Collect nested scientific artifact identities, rejecting contradictions."""
    found = {}

    def visit(item):
        if isinstance(item, dict):
            if {"path", "bytes", "sha256"} <= item.keys():
                normalized = {k: item[k] for k in ("path", "bytes", "sha256")}
                if item["path"] in found and found[item["path"]] != normalized:
                    raise ValueError("Conflicting artifact identities")
                found[item["path"]] = normalized
            else:
                for child in item.values():
                    visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)

    visit(value)
    return list(found.values())


def verify_completion(execution, result, plan_record, result_record, admission_record):
    if (execution != dict(status="readback_complete", job_id=str(READBACK_JOB),
                          plan=plan_record, result=result_record)
            or result["status"] != "fresh_qfo_phylogeny_scientific_readback_complete"
            or result["job_id"] != 22329 or result["admission"] != admission_record
            or result["accuracy_evaluated"] is not False
            or result["historical_scores_replaced"] is not False):
        raise ValueError("Unbound or incomplete readback")


def admit(repo, directory):
    accounting = subprocess.check_output(["sacct", "-j", str(READBACK_JOB), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, READBACK_JOB)
    if scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong readback CPU allocation")
    path = directory / "readback_v2_plan.json"
    if record(path)["sha256"] != READBACK_PLAN_SHA:
        raise ValueError("Changed readback plan")
    plan = json.loads(path.read_text())
    if (plan["repo"] != str(repo) or plan["directory"] != str(directory)
            or plan["output"] != str(directory / "readback_v2") or plan["native_job_id"] != 22329):
        raise ValueError("Wrong readback locations or native job")
    frozen = repo / "benchmark_tools/results/qfo_fresh_phylogeny_readback_v2_submission_20260927.json"
    if record(frozen)["sha256"] != SUBMISSION_SHA:
        raise ValueError("Changed frozen readback submission")
    submission_path = directory / "readback_v2_submission.json"
    submission = json.loads(submission_path.read_text())
    if (submission != json.loads(frozen.read_text()) or submission["job_id"] != str(READBACK_JOB)
            or submission["plan"] != record(path)):
        raise ValueError("Wrong readback submission")
    native = admit_native(repo, directory)
    output = directory / "readback_v2"
    result_path, admission_path = output / "result.json", output / "admission.json"
    result = json.loads(result_path.read_text())
    execution_path = directory / "readback_v2_execution.json"
    verify_completion(json.loads(execution_path.read_text()), result, record(path),
                      record(result_path), record(admission_path))
    if native != json.loads(admission_path.read_text()) or result["plan"] != native["plan"]:
        raise ValueError("Native admission changed since scientific readback")
    expected_reports = [output / arm / (name + ".json") for arm in ("historical", "fresh") for name in REPORTS]
    if sorted(result["reports"], key=lambda r: r["path"]) != [record(p) for p in sorted(expected_reports)]:
        raise ValueError("Missing, extra or altered scientific report")
    reports = []
    for report_path in expected_reports:
        report = json.loads(report_path.read_text())
        if report["status"] != REPORTS[report_path.stem]:
            raise ValueError("Unverified scientific report")
        if report_path.stem == "structure" and (report["genes"] != 984137 or report["species"] != 78):
            raise ValueError("Wrong scientific universe")
        reports.append(report)
    native_plan = json.loads((directory / "plan.json").read_text())
    verified_sources, archived = audited_records(repo, directory, plan["checked_records"])
    checked = collect_records([verified_sources, native_plan["checked_records"], native,
        reports, result, *[record(p) for p in (path, frozen, submission_path, result_path, admission_path,
                                             execution_path, Path(__file__))]])
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed reuse dependency: " + item["path"])
    return dict(status="fresh_qfo_raw_tree_reuse_admitted", native_job=22329,
        readback_job=READBACK_JOB, scheduler=scheduler, checked_records=checked,
        audited_orchestration_archive=archived,
        source_directory=str(directory / "native/inference/orthohmm_phylogeny"),
        admitted_outputs=native["outputs"], native_plan=record(directory / "plan.json"),
        readback_result=record(result_path), historical_equivalence_required=False,
        scope="Eligibility only; each copied family still requires matching tree inputs and tools",
        accuracy_evaluated=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = admit(args.repo.resolve(), args.directory.resolve())
    save(args.output, result)
    print(json.dumps(dict(status=result["status"], output=record(args.output)), indent=2))
