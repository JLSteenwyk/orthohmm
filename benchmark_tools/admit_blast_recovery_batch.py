"""Validate a completed recovery batch, retaining failures and indexed query blocks."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

from Bio import SeqIO
from benchmark_tools.audit_orthomcl_blast import parse_diagnostics
from benchmark_tools.audit_orthomcl_search_table import audit_table
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

EXECUTOR_COMMIT = "1a73dc619bf50d5cc37b59f7b7820df3f4787c5d"
ARRAY_JOB = "22103"
REPLACEMENT_ARRAY = "22160"
INTERRUPTED_HASHES = {
    "hits.blast.partial": "6e22696fd4050d15933401c1d28b4e72c9eb585083c51ea08f5795646b51dabf",
    "preflight.json": "65c3ef26d3c98ca45a11a0be8e9c0f0f84f5ba36f9a4055fd55a0313aa858dc9",
    "status.json": "65c3ef26d3c98ca45a11a0be8e9c0f0f84f5ba36f9a4055fd55a0313aa858dc9",
    "blast.log": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "blast.time.txt": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
}


def execution_array(index, replacement=False):
    if type(index) is not int or index not in range(20):
        raise ValueError("Unknown recovery batch")
    if type(replacement) is not bool or (replacement and index != 14):
        raise ValueError("Replacement authorized only for batch 14")
    return REPLACEMENT_ARRAY if replacement else ARRAY_JOB


def interrupted_attempt(root, accounting):
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == "22103_14"]
    if len(rows) != 1 or tuple(rows[0][k] for k in
            ("JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "22143", "TIMEOUT", "0:0", "bizon", "180"):
        raise ValueError("Interrupted original task identity changed")
    directory = root / "benchmarks/results/qfo_blast_recovery_v1/batch_14_interrupted_22103_14_20260925"
    if {p.name for p in directory.iterdir()} != set(INTERRUPTED_HASHES):
        raise ValueError("Interrupted attempt inventory changed")
    records = [record(directory / name) for name in INTERRUPTED_HASHES]
    if any(r["sha256"] != INTERRUPTED_HASHES[Path(r["path"]).name] for r in records):
        raise ValueError("Interrupted attempt bytes changed")
    return dict(original_scheduler=rows[0], records=records,
                replacement_task="22160_14", partial_rows_reused=False)


def completed_task(accounting, index, replacement=False):
    array = execution_array(index, replacement)
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == f"{array}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != ("COMPLETED", "0:0", "bizon", "180"):
        raise ValueError("Require uniquely completed matching recovery task")
    return rows[0]


def index_blocks(blast, genes):
    if not genes or len(set(genes)) != len(genes):
        raise ValueError("Empty or duplicate batch queries")
    ordinals = {gene: i for i, gene in enumerate(genes)}
    blocks, current, start, end, rows = [], None, 0, 0, 0
    digest = hashlib.sha256()
    previous = -1

    def append():
        blocks.append(dict(query=current, batch_query_ordinal=ordinals[current], start=start, end=end,
                           rows=rows, sha256=digest.hexdigest()))

    with blast.open("rb") as stream:
        for line in stream:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if not line.endswith(b"\n") or len(fields) != 12 or b"\x00" in line:
                raise ValueError("Incomplete or malformed recovery HSP")
            query = fields[0].decode("ascii")
            if query not in ordinals:
                raise ValueError("Recovery output query outside its batch")
            if query != current:
                if ordinals[query] <= previous:
                    raise ValueError("Repeated or out-of-order recovery query block")
                if current is not None:
                    append()
                current, start, rows, digest = query, end, 0, hashlib.sha256()
                previous = ordinals[query]
            digest.update(line)
            end += len(line)
            rows += 1
    if current is not None:
        append()
    return blocks


def coverage(genes, blocks, diagnostics):
    universe = set(genes)
    hits = {b["query"] for b in blocks}
    if not set(diagnostics) <= universe or not hits <= universe:
        raise ValueError("Query diagnostics or blocks outside batch")
    failed = {g for g, r in diagnostics.items() if r["query_failed"]}
    return dict(input_queries=len(genes), queries_with_hits=len(hits),
        queries_without_hits=sorted(universe-hits), failed_queries=sorted(failed),
        failed_queries_with_outgoing_hits=sorted(failed & hits),
        no_hits_without_logged_failure=sorted(universe-hits-failed))


def admit(root, index, output, replacement=False):
    root, output = root.resolve(), output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    array = execution_array(index, replacement)
    accounting = subprocess.check_output(["sacct", "-j", f"{ARRAY_JOB},{array}" if replacement else array, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,NodeList,AllocCPUS,Elapsed"], text=True)
    scheduler = completed_task(accounting, index, replacement)
    recovery = interrupted_attempt(root, accounting) if replacement else None
    executor = root / "benchmarks/work/blast_recovery_executor_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Recovery executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    # Re-derive the frozen execution contract without launching native search.
    code = "import json,sys; from pathlib import Path; from benchmark_tools.run_blast_recovery_batch import prepare; print(json.dumps(prepare(Path(sys.argv[1]), int(sys.argv[2]))))"
    expected = json.loads(subprocess.check_output(["/home/bizon/anaconda3/bin/python", "-B", "-c", code,
        str(root), str(index)], cwd=executor, text=True))
    directory = Path(expected["directory"])
    status_path, preflight_path = directory / "status.json", directory / "preflight.json"
    status, initial = [json.loads(p.read_text()) for p in (status_path, preflight_path)]
    if (status["status"] != "native_batch_completed_pending_admission" or status["exit_code"] != 0
            or status["preflight"] != expected or initial["preflight"] != expected or initial["status"] != "starting"
            or status["index"] != index or status["array_job_id"] != array
            or status["job_id"] != scheduler["JobIDRaw"] or status["node"] != "bizon"
            or status["finished_epoch"] < status["started_epoch"]
            or any(status[k] is not False for k in ("search_admitted", "reuse_authorized", "publication_ready"))):
        raise ValueError("Terminal native status contradicts execution contract")
    if any(initial[k] != status[k] for k in ("index", "job_id", "array_job_id", "node", "started_epoch")):
        raise ValueError("Preflight execution identity changed")
    blast, log, timing = [directory / name for name in ("hits.blast", "blast.log", "blast.time.txt")]
    if status["outputs"] != [record(p) for p in (log, timing, blast)]:
        raise ValueError("Native output identities changed")
    actual = {p.name for p in directory.iterdir()}
    if actual != {"preflight.json", "status.json", "hits.blast", "blast.log", "blast.time.txt"}:
        raise ValueError("Unexpected native batch artifact inventory")
    checked = [record(status_path), record(preflight_path), record(__file__),
        *expected["checked_records"], *status["outputs"],
        *[record(Path(__file__).with_name(name)) for name in (
            "audit_orthomcl_blast.py", "audit_orthomcl_search_table.py", "convert_orthomcl_blast.py")]]
    if recovery:
        checked.extend(recovery["records"])
    for item in checked:
        check(item)
    genes = [r.id for r in SeqIO.parse(expected["batch"]["input"]["path"], "fasta")]
    batch = expected["batch"]
    if len(genes) != batch["queries"] or genes[0] != batch["first_query"] or genes[-1] != batch["last_query"]:
        raise ValueError("Batch query inventory changed")
    blocks = index_blocks(blast, genes)
    diagnostics = parse_diagnostics(log)
    full_fasta = Path(expected["command"][expected["command"].index("-d")+1])
    structural = audit_table(blast, full_fasta, log) if blast.stat().st_size else None
    if structural and (structural["hsp_rows_above_1e_minus_5"] or structural["hsp_rows"] != sum(b["rows"] for b in blocks)):
        raise ValueError("Invalid cutoff or row accounting")
    if structural:
        # Whole-database absent-query lists are not batch completion evidence.
        structural.pop("query_ids_without_hits")
    observed = coverage(genes, blocks, diagnostics)
    if observed["failed_queries_with_outgoing_hits"]:
        raise ValueError("Logged failed query has partial outgoing hits; review required")
    for item in checked:
        check(item)
    report = dict(status="recovery_batch_execution_and_rows_verified", index=index, scheduler=scheduler,
        accounting=accounting, executor_commit=EXECUTOR_COMMIT, records=checked,
        query_ids=genes, query_blocks=blocks, coverage=observed,
        full_database_structural_audit=structural, diagnostics=diagnostics,
        batch_admitted=True, search_admitted=False, reuse_authorized=False, publication_ready=False,
        limitations=["Single completed batch only; not full recovery coverage, prefix reuse, or whole-search admission.",
            "Logged failed queries remain explicit despite native exit zero; no-hit queries are not automatically failures.",
            "Native execution contract is re-derived using frozen preparer code; row indexing and coverage are independently checked."])
    if recovery:
        report["replacement"] = recovery
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(20), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replacement", action="store_true")
    args = parser.parse_args()
    admit(args.root, args.index, args.output, args.replacement)
