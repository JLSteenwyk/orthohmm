"""Check complete admitted replay coverage; never merge or authorize prefix reuse."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess

from Bio import SeqIO
from benchmark_tools.admit_blast_recovery_batch import completed_task, coverage, EXECUTOR_COMMIT, interrupted_attempt
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

BATCH_SHA = "5e59e7deedb71feb6453057f0b379bda0e06a3556651b13eb7275b26155f45e0"
PARTITION_SHA = "c37381bef1cdb406849f62accde679480a4abe76193cc0fdc65a90e16c8d6a34"
ADMIT_SHA = "b8081b449ad17adc0ba16dc964249e5d068a303063e0b7c921a00326b0c4062d"
ADMIT_COMMIT = "7306c5540e579854800c7a944b570ad478124485"
REPLACEMENT_SHA = "5b5d4abb0aaf08d88e56d1df5b74a306445fe522652d8b03282fc43752767c98"
REPLACEMENT_COMMIT = "55e20ecc4eed86766b42f9c5a6ff1c8cf214f15f"


def validate_batch_provenance(report, accounting, index, source, replacement=False, recovery=None):
    replaced = replacement and index == 14
    validator_id = "22161" if replaced else f"22105_{index}"
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    validators = [r for r in rows if r["JobID"] == validator_id]
    if len(validators) != 1 or tuple(validators[0][k] for k in (
            "State", "ExitCode", "NodeList", "AllocCPUS")) != ("COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require uniquely completed matching admission job")
    sources = [r for r in report["records"] if Path(r["path"]).name == "admit_blast_recovery_batch.py"]
    if sources != [source]:
        raise ValueError("Batch was not checked by the frozen admission source")
    if report["scheduler"] != completed_task(accounting, index, replaced):
        raise ValueError("Native scheduler identity differs from batch admission")
    if replaced:
        if recovery is None or report.get("replacement") != recovery:
            raise ValueError("Replacement failed-attempt evidence differs")
        if any(r not in report["records"] for r in recovery["records"]):
            raise ValueError("Replacement original records missing")
    elif "replacement" in report:
        raise ValueError("Unexpected replacement provenance")
    return validators[0]


def unique_records(records):
    by_path = {}
    for item in records:
        if item["path"] in by_path and by_path[item["path"]] != item:
            raise ValueError("Conflicting provenance for one path")
        by_path[item["path"]] = item
    return list(by_path.values())


def combine(reports, expected_batches, replay_ids):
    if len(reports) != len(expected_batches) or not reports:
        raise ValueError("Incomplete recovery panel")
    seen, diagnostics, failures, no_hits, rows, summaries = [], {}, [], [], 0, []
    for index, (report, genes) in enumerate(zip(reports, expected_batches)):
        if (type(report["index"]) is not int or report["index"] != index
                or report["status"] != "recovery_batch_execution_and_rows_verified"
                or report["batch_admitted"] is not True
                or any(report[k] is not False for k in ("search_admitted", "reuse_authorized", "publication_ready"))
                or report["executor_commit"] != EXECUTOR_COMMIT
                or not genes or report["query_ids"] != genes or len(set(genes)) != len(genes)):
            raise ValueError("Wrong batch admission or query order")
        block_ids = [b["query"] for b in report["query_blocks"]]
        block_set = set(block_ids)
        if len(block_set) != len(block_ids) or block_ids != [g for g in genes if g in block_set]:
            raise ValueError("Invalid query block order")
        observed = coverage(genes, report["query_blocks"], report["diagnostics"])
        if observed != report["coverage"] or observed["failed_queries_with_outgoing_hits"]:
            raise ValueError("Changed coverage or failed-query partial output")
        if set(diagnostics) & set(report["diagnostics"]):
            raise ValueError("Diagnostic query repeated across batches")
        seen.extend(genes)
        diagnostics.update(report["diagnostics"])
        failures.extend(observed["failed_queries"])
        no_hits.extend(observed["no_hits_without_logged_failure"])
        rows += sum(b["rows"] for b in report["query_blocks"])
        summaries.append(dict(index=index, coverage=observed))
    if len(set(seen)) != len(seen) or seen != replay_ids:
        raise ValueError("Replay universe duplicated, reordered or incomplete")
    return dict(replay_queries=len(seen), hsp_rows=rows, failed_queries=sorted(failures),
                no_hits_without_logged_failure=sorted(no_hits), diagnostics=diagnostics,
                batches=summaries)


def verify(root, output, replacement=False):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    batch_path = root / "benchmark_tools/results/qfo_blast_recovery_batches_20260923.json"
    partition_path = root / "benchmark_tools/results/qfo_blast_recovery_partition_20260923.json"
    batches = read_frozen(batch_path, BATCH_SHA)
    partition = read_frozen(partition_path, PARTITION_SHA)
    if batches["total_queries"] != 98913 or len(batches["batches"]) != 20 or partition["replay_queries"] != 98913:
        raise ValueError("Wrong frozen replay panel")
    executor = root / "benchmarks/work/blast_recovery_admission_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMIT_COMMIT:
        raise ValueError("Wrong admission executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    admit_source = record(executor / "benchmark_tools/admit_blast_recovery_batch.py")
    if admit_source["sha256"] != ADMIT_SHA or record(Path(__file__).with_name("admit_blast_recovery_batch.py"))["sha256"] != REPLACEMENT_SHA:
        raise ValueError("Changed admission or coverage implementation")
    replacement_source = None
    if replacement:
        replacement_executor = root / "benchmarks/work/blast_replacement_admission_v1_20260925"
        if subprocess.check_output(["git", "-C", str(replacement_executor), "rev-parse", "HEAD"], text=True).strip() != REPLACEMENT_COMMIT:
            raise ValueError("Wrong replacement admission executor")
        if subprocess.check_output(["git", "-C", str(replacement_executor), "status", "--porcelain"], text=True).strip():
            raise ValueError("Dirty replacement admission executor")
        replacement_source = record(replacement_executor / "benchmark_tools/admit_blast_recovery_batch.py")
        if replacement_source["sha256"] != REPLACEMENT_SHA:
            raise ValueError("Changed replacement admission source")
    accounting = subprocess.check_output(["sacct", "-j", "22103,22105,22160,22161" if replacement else "22103,22105", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,NodeList,AllocCPUS,Elapsed"], text=True)
    recovery = interrupted_attempt(root, accounting) if replacement else None
    reports, expected, records, admissions = [], [], [record(__file__), record(batch_path), record(partition_path),
        admit_source, *batches["inputs"], *partition["outputs"]], []
    for index, batch in enumerate(batches["batches"]):
        path = root / f"benchmarks/work/qfo_blast_recovery_admission_{index}_20260923.json"
        item = record(path)
        report = json.loads(path.read_text())
        source = replacement_source if replacement and index == 14 else admit_source
        validator = validate_batch_provenance(report, accounting, index, source, replacement, recovery)
        if batch["index"] != index or batch["replay_ordinal_start"] != sum(len(g) for g in expected):
            raise ValueError("Wrong frozen batch ordinal")
        check(batch["input"])
        genes = [r.id for r in SeqIO.parse(batch["input"]["path"], "fasta")]
        if (len(genes) != batch["queries"] or genes[0] != batch["first_query"] or genes[-1] != batch["last_query"]
                or batch["replay_ordinal_end_exclusive"] != batch["replay_ordinal_start"] + len(genes)):
            raise ValueError("Frozen FASTA and batch manifest disagree")
        expected.append(genes)
        reports.append(report)
        records.extend([item, batch["input"], *report["records"]])
        admissions.append(dict(index=index, report=item, scheduler=validator))
    records = unique_records(records)
    for item in records:
        check(item)
    replay = [r for r in partition["outputs"] if Path(r["path"]).name == "replay.fa"]
    if len(replay) != 1:
        raise ValueError("Require one frozen replay FASTA")
    replay_ids = [r.id for r in SeqIO.parse(replay[0]["path"], "fasta")]
    result = combine(reports, expected, replay_ids)
    if result["replay_queries"] != 98913:
        raise ValueError("Wrong total replay coverage")
    for item in records:
        check(item)
    report = dict(status="complete_replay_panel_checked_no_prefix_reuse", coverage=result,
        admissions=admissions, accounting=accounting, records=records,
        replay_panel_admitted=True, search_admitted=False, reuse_authorized=False, publication_ready=False,
        limitations=["Retained-prefix evidence review and complete-table merge/validation remain required.",
                     "This aggregates independently admitted batches; it does not rerun native search or its row parser.",
                     "Failed and no-hit queries remain distinct; incoming subject hits are not outgoing search success."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replacement", action="store_true")
    args = parser.parse_args()
    verify(args.root.resolve(), args.output.absolute(), args.replacement)
