"""Describe completed replay batches without admitting a whole recovered search."""

import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_blast_recovery_batch import completed_task, coverage, EXECUTOR_COMMIT
from benchmark_tools.audit_orthomcl_blast import parse_diagnostics
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def summarize(paths, accounting):
    if not paths:
        raise ValueError("Require explicit completed batch admissions")
    checked, batches, failed, seen, indices = [], [], [], set(), set()
    for path in paths:
        item = record(path)
        report = json.loads(path.read_text())
        index = report["index"]
        if (report["status"] != "recovery_batch_execution_and_rows_verified"
                or report["batch_admitted"] is not True
                or report["executor_commit"] != EXECUTOR_COMMIT or index in indices):
            raise ValueError("Wrong or repeated admitted batch")
        scheduler = completed_task(accounting, index)
        if scheduler != report["scheduler"]:
            raise ValueError("Batch accounting differs from admission")
        indices.add(index)
        selected = []
        for name in (f"queries_{index:02d}.fa", "blast.log"):
            matches = [r for r in report["records"] if Path(r["path"]).name == name]
            if len(matches) != 1:
                raise ValueError("Missing or ambiguous batch FASTA/log record")
            check(matches[0])
            selected.append(matches[0])
        entries = list(SeqIO.parse(selected[0]["path"], "fasta"))
        genes = [entry.id for entry in entries]
        if (not genes or len(set(genes)) != len(genes) or seen.intersection(genes)
                or genes != report["query_ids"] or any(not len(entry) for entry in entries)):
            raise ValueError("Changed, empty or overlapping batch queries")
        seen.update(genes)
        diagnostics = parse_diagnostics(Path(selected[1]["path"]))
        if diagnostics != report["diagnostics"]:
            raise ValueError("Diagnostic log differs from admission")
        observed = coverage(genes, report["query_blocks"], diagnostics)
        if observed != report["coverage"] or observed["failed_queries_with_outgoing_hits"]:
            raise ValueError("Changed coverage or failed-query partial output")
        sequences = {entry.id: str(entry.seq).upper() for entry in entries}
        for gene in observed["failed_queries"]:
            sequence = sequences[gene]
            failed.append(dict(gene=gene, batch=index, length=len(sequence),
                residue_counts=dict(Counter(sequence)),
                categories=sorted({m["category"] for m in diagnostics[gene]["messages"]}),
                has_outgoing_hits=False))
        batches.append(dict(index=index, scheduler=scheduler, input_queries=len(genes),
            queries_with_hits=observed["queries_with_hits"], failed_queries=len(observed["failed_queries"]),
            no_hits_without_logged_failure=len(observed["no_hits_without_logged_failure"])))
        checked.extend([item, *selected])
    for item in checked:
        check(item)
    lengths = [row["length"] for row in failed]
    return dict(status="partial_replay_failure_description_not_search_admission",
        batches=sorted(batches, key=lambda row: row["index"]), input_queries=len(seen),
        queries_with_hits=sum(row["queries_with_hits"] for row in batches), failed_queries=len(failed),
        no_hits_without_logged_failure=sum(row["no_hits_without_logged_failure"] for row in batches),
        failed_query_lengths=dict(min=min(lengths), median=statistics.median(lengths), max=max(lengths)) if lengths else None,
        failed_queries_by_category=dict(Counter(c for row in failed for c in row["categories"])),
        records=sorted(failed, key=lambda row: row["gene"]), checked_records=checked,
        search_admitted=False, accuracy_admitted=False, publication_ready=False,
        limitations=[
            "Only explicitly selected admitted replay batches; not the original prefix or the whole dataset.",
            "Replay queries are selected by interrupted-search recovery, not a random sample; do not extrapolate failure rates.",
            "Sequence lengths and logged categories describe failures, not experimentally established causes.",
            "Categories overlap: setup failure can accompany a more specific diagnostic.",
            "No-hit without a logged failure is not classified as failure.",
            "Query-block inventories are reused from admission; HSP bytes and database parity are not revalidated here.",
            "Incoming hits, final-group membership, reference exposure and counterfactual accuracy require later analysis."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    accounting = subprocess.check_output(["sacct", "-j", "22103", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,NodeList,AllocCPUS,Elapsed"], text=True)
    result = summarize([p.resolve() for p in args.batch], accounting)
    result.update(source=record(__file__), helpers=[record(Path(__file__).with_name(name)) for name in (
        "admit_blast_recovery_batch.py", "audit_orthomcl_blast.py", "prepare_ob_candidate_neighborhood.py")])
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
