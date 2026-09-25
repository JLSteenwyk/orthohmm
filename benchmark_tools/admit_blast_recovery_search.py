"""Admit recovered search evidence without releasing downstream inference."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_blast_recovery_candidate import audit as audit_candidate
from benchmark_tools.audit_orthomcl_database import audit as audit_database
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_blast_recovery_merge import prepare, diagnostic_signature
from benchmark_tools.run_qfo_corrected_blast import verify
from benchmark_tools.verify_blast_recovery_panel import unique_records

MERGE_COMMIT = "a449ff580aca58e8d1fdc653e7615e007093771a"
MERGE_JOB = "22150"
REPLACEMENT_MERGE_JOB = "22162"
REPLACEMENT_MERGE_COMMIT = "f6ab36db0cbed151ac2b46581345cd187d514e79"


def completed_merge(accounting, replacement=False):
    job = REPLACEMENT_MERGE_JOB if replacement else MERGE_JOB
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == job]
    if len(rows) != 1 or tuple(rows[0][k] for k in (
            "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem")) != (
            "COMPLETED", "0:0", "bizon", "2", "64G"):
        raise ValueError("Require completed frozen recovery merge allocation")
    return rows[0]


def validate_merge(status, source, directory, replacement=False):
    job = REPLACEMENT_MERGE_JOB if replacement else MERGE_JOB
    if (status["status"] != "merged_candidate_pending_full_table_admission"
            or status["job_id"] != job or status["source"] != source
            or any(status[k] is not False for k in (
                "search_admitted", "reuse_authorized", "publication_ready"))):
        raise ValueError("Wrong merge execution identity or status")
    candidate = status["candidate"]
    if (candidate["status"] != "merged_candidate_requires_full_admission"
            or candidate["path"] != str(directory / "table/all.blast.candidate")
            or status["selected_log"]["path"] != str(directory / "selected.blast.log")
            or any(candidate[k] is not False for k in (
                "search_admitted", "reuse_authorized", "publication_ready"))):
        raise ValueError("Wrong merge candidate or diagnostic output")
    return {key: candidate[key] for key in ("path", "bytes", "sha256")}


def admit(root, output, replacement=False, reviewed_native_o=False):
    if reviewed_native_o and not replacement:
        raise ValueError("Reviewed native representation requires replacement recovery")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    job = REPLACEMENT_MERGE_JOB if replacement else MERGE_JOB
    commit = REPLACEMENT_MERGE_COMMIT if replacement else MERGE_COMMIT
    accounting = subprocess.check_output(["sacct", "-j", job, "--parsable2",
        "--format=JobID,State,ExitCode,NodeList,AllocCPUS,ReqMem,Elapsed"], text=True)
    scheduler = completed_merge(accounting, replacement)
    executor = root / ("benchmarks/work/blast_replacement_merge_v1_20260925" if replacement
                       else "benchmarks/work/blast_recovery_merge_v1_20260923")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed merge executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    frozen_source = record(executor / "benchmark_tools/run_blast_recovery_merge.py")
    local_source = record(Path(__file__).with_name("run_blast_recovery_merge.py"))
    if local_source["sha256"] != frozen_source["sha256"]:
        raise ValueError("Merge prerequisite revalidator differs from executed source")
    directory = root / ("benchmarks/results/qfo_blast_replacement_merge_v1" if replacement
                        else "benchmarks/results/qfo_blast_recovery_merge_v1")
    status_path = directory / "status.json"
    status_record = record(status_path)
    status = json.loads(status_path.read_text())
    candidate_record = validate_merge(status, frozen_source, directory, replacement)
    checked = unique_records([record(__file__), status_record, frozen_source, local_source,
        candidate_record, status["selected_log"], *status["checked_inputs"]])
    for item in checked:
        check(item)
    output.mkdir(exist_ok=False)
    report = dict(status="validating_recovered_search", scheduler=scheduler, accounting=accounting,
        checked_records=checked, search_admitted=False, accuracy_admitted=False,
        publication_ready=False, downstream_execution_authorized=False)
    save_status(output / "report.json", report)
    try:
        # Re-run the same scientific prerequisites without copying or merging.
        fresh = output / "fresh_prerequisites"
        fresh.mkdir()
        prepared = prepare(root, fresh, replacement)
        expected_diagnostics = {g: diagnostic_signature(d) for g, d in prepared["selected_diagnostics"].items()}
        if {g: diagnostic_signature(d) for g, d in status["diagnostics"].items()} != expected_diagnostics:
            raise ValueError("Merged diagnostic selection differs from fresh prerequisites")
        if status["replay_coverage"] != prepared["panel"]["coverage"]:
            raise ValueError("Merged replay coverage differs from fresh prerequisites")
        plan_path = root / "benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json"
        runtime_path = root / "benchmark_tools/results/qfo_corrected_legacy_blast_runtime_20260918.json"
        plan = verify(plan_path, runtime_path)
        fasta = Path(plan["output_root"]) / "work/all.fa"
        database = audit_database(fasta, runtime_path, output / "database")
        representation = None
        if reviewed_native_o:
            from benchmark_tools.reviewed_legacy_database import verify as verify_representation
            representation = verify_representation(root, database)
            report["database_representation"] = representation
            checked = unique_records([*checked, *representation["checked_records"]])
        elif (database["status"] != "database_exact_sequence_parity_verified"
                or database["content"]["input_sequences"] != 984137):
            raise ValueError("Recovered search database lacks exact corrected-input parity")
        table = audit_candidate(Path(candidate_record["path"]), fasta,
            Path(status["selected_log"]["path"]), Path(prepared["prefix_audit"]["blocks"]["path"]),
            [Path(a["report"]["path"]) for a in prepared["panel"]["admissions"]], output / "table.json")
        if (table["content"]["input_proteins"] != 984137
                or table["content"]["hsp_rows"] != status["candidate"]["rows"]
                or table["query_blocks"] != status["candidate"]["query_blocks"]):
            raise ValueError("Candidate totals disagree with merged output")
        actual_diagnostics = {d["gene"]: diagnostic_signature(d) for d in table["content"]["diagnostics"]}
        if actual_diagnostics != expected_diagnostics:
            raise ValueError("Final diagnostic log differs from freshly selected evidence")
        if verify(plan_path, runtime_path) != plan:
            raise ValueError("Prepared search provenance changed during admission")
        checked = unique_records([*checked, *prepared["records"], *database["checked_records"],
            *database["outputs"], *table["checked_records"], record(plan_path), record(runtime_path),
            record(output / "database/report.json"), record(output / "table.json")])
        for item in checked:
            check(item)
        report.update(status="recovered_orthomcl_search_evidence_verified", search_admitted=True,
            checked_records=checked, query_coverage=table["content"], database_content=database["content"],
            candidate=candidate_record, selected_log=status["selected_log"],
            limitations=["Recovered search, not a successful uninterrupted original BLAST run.",
                "Conditional prefix reuse retains historical durability and source-build uncertainty.",
                "Logged failed queries remain explicit; incoming hits do not repair outgoing failures.",
                "BPO, clustering, failure-impact analysis and scoring require separate admission.",
                "No downstream job release or matched-resource timing claim is authorized."])
        if representation is not None:
            report["status"] = "recovered_search_native_representation_verified"
            report["limitations"].extend(representation["limitations"])
    except BaseException as error:
        report.update(status="recovery_search_admission_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replacement", action="store_true")
    parser.add_argument("--reviewed-native-o", action="store_true",
                        help="Require the exact pinned seven-deletion native representation review")
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.absolute(), args.replacement, args.reviewed_native_o)
