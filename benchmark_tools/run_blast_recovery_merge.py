"""Gate and write a recovery candidate; never admit it to OrthoMCL inference."""

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

from Bio import SeqIO
from benchmark_tools.audit_orthomcl_blast import DIAGNOSTIC, parse_diagnostics
from benchmark_tools.merge_blast_recovery_blocks import ordered_blocks, merge
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_recovery_batch import save_status, sync_directory
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_blast_recovery_panel import verify, unique_records

PINS = {
    "QFO_BLAST_PREFIX_REUSE_REVIEW_20260923.md": "14e8eb1a7a02bc77e3712c741da6c61a23611722852906a3b716b1b7fe94a24f",
    "qfo_blast_prefix_recheck_22148.json": "315913c456280b55a5c591b7e4c138b5af58fc3c62a783fadde056473362d0d6",
    "qfo_blast_recovery_partition_20260923.json": "c37381bef1cdb406849f62accde679480a4abe76193cc0fdc65a90e16c8d6a34",
    "qfo_blast_replay_comparison_22055.json": "9d49960fd22aaa1b538dcae728b06d0a67802d8d3d1f4296b31c459a98f702e1",
}
HELPERS = {
    "verify_blast_recovery_panel.py": "8f40db306ccbc567dc1ed5de1a6b9828fb99a26df57991df3db9b3bf16a99465",
    "merge_blast_recovery_blocks.py": "3b9814d13dedfea369e1deb800ce40e86c87f3f83a8b2016a1b64abc48cbb3e1",
}
ARCHIVES = {
    "ncbi.tar.gz": "3f311eb066a49c73eef36cdb8305c97c18ed35f7f9bc2f8904e18d4298e5a85d",
    "blast-2.2.13-x64-linux.tar.gz": "e284b4f95adf52267b21d3cbdffb9987c8ef038c557e0b3527e54fa97796bb0c",
    "MD5SUM.txt": "8f6c4abdf3396fef4b0e11d58ecfb0450f2954b68296d40efc707cc2fca09948",
    "ncbi/demo/blastall.c": "3344597098b9ce248be14cab885646a9a9604c7e6021f143b4ede356e172ff60",
}


def require_completed(accounting, replacement=False):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    expected = {"22148": "2", **{f"{job}_{i}": cpu for job, cpu in
                (("22103", "180"), ("22105", "2")) for i in range(20)}}
    if replacement:
        del expected["22103_14"], expected["22105_14"]
        expected.update({"22160_14": "180", "22161": "2"})
    for job, cpu in expected.items():
        selected = [r for r in rows if r["JobID"] == job]
        if len(selected) != 1 or tuple(selected[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
                "COMPLETED", "0:0", "bizon", cpu):
            raise ValueError("Incomplete or failed prerequisite: " + job)


def diagnostic_signature(item):
    return (item["query_failed"], Counter((m["level"], m["category"], m["message"]) for m in item["messages"]))


def select_diagnostics(original, replay, replay_ids):
    replay_ids = set(replay_ids)
    if not set(replay) <= replay_ids:
        raise ValueError("Replay diagnostics outside replay universe")
    for gene, item in original.items():
        if gene in replay_ids:
            if gene not in replay or diagnostic_signature(item) != diagnostic_signature(replay[gene]):
                raise ValueError("Repeated original diagnostic changed")
        elif item["query_failed"]:
            raise ValueError("Original failed query cannot use prefix output")
    return {**{g: r for g, r in original.items() if g not in replay_ids}, **replay}


def raw_diagnostics(path):
    result = {}
    with path.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            match = DIAGNOSTIC.fullmatch(line.rstrip("\r\n"))
            if match is None or not line.endswith("\n"):
                raise ValueError("Invalid raw diagnostic line")
            result.setdefault(match.group(2), []).append(line)
    return result


def prepare(root, output, replacement=False):
    accounting = subprocess.check_output(["sacct", "-j", "22148,22103,22105,22160,22161" if replacement else "22148,22103,22105", "--parsable2",
        "--format=JobID,State,ExitCode,NodeList,AllocCPUS,Elapsed"], text=True)
    require_completed(accounting, replacement)
    results = root / "benchmark_tools/results"
    records = [record(__file__)]
    for name, digest in {**PINS, **HELPERS}.items():
        path = results / name if name in PINS else Path(__file__).with_name(name)
        item = record(path)
        if item["sha256"] != digest:
            raise ValueError("Changed merge gate: " + name)
        records.append(item)
    recheck = read_frozen(results / "qfo_blast_prefix_recheck_22148.json", PINS["qfo_blast_prefix_recheck_22148.json"])
    if recheck["status"] != "prefix_block_bytes_rechecked_not_admitted" or recheck["scheduler_job_id"] != "22148":
        raise ValueError("Wrong prefix byte recheck")
    executor = root / "benchmarks/work/blast_prefix_recheck_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != "d61769462e67151482f2e817dc9926b26cbb33b1":
        raise ValueError("Changed prefix recheck executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if recheck["source"] != record(executor / "benchmark_tools/check_blast_prefix_bytes.py"):
        raise ValueError("Wrong prefix recheck source identity")
    partition = read_frozen(results / "qfo_blast_recovery_partition_20260923.json", PINS["qfo_blast_recovery_partition_20260923.json"])
    comparison = read_frozen(results / "qfo_blast_replay_comparison_22055.json", PINS["qfo_blast_replay_comparison_22055.json"])
    if comparison["status"] != "diagnostic_replay_comparison_complete" or comparison["all_diagnostics_compatible"] is not True:
        raise ValueError("Diagnostic replay incompatible")
    records.extend([recheck["source"], *recheck["checked_inputs"], *partition["inputs"],
                    *partition["outputs"], *comparison["records"]])
    archive_dir = root / "benchmarks/work/legacy_blast_source_20260923"
    for name, digest in ARCHIVES.items():
        item = record(archive_dir / name)
        if item["sha256"] != digest:
            raise ValueError("Changed reviewed source archive")
        records.append(item)
    binary = root.parents[1] / "SOFTWARE/blast-2.2.13/bin/blastall"
    binary_record = record(binary)
    with tarfile.open(archive_dir / "blast-2.2.13-x64-linux.tar.gz", "r:gz") as archive:
        content = archive.extractfile("blast-2.2.13/bin/blastall").read()
    if (hashlib.sha256(content).hexdigest() != binary_record["sha256"]
            or binary_record["sha256"] != "34b91fedd2e7858a3478e83758e3f5339c09f96f543d9338488ea4b84628fb75"):
        raise ValueError("Installed binary differs from reviewed release")
    records.append(binary_record)
    records = unique_records(records)
    for item in records:
        check(item)
    # Re-run the complete panel gate, rather than trusting a user-supplied summary.
    panel = verify(root, output / "replay_panel.json", replacement)
    reports = [json.loads(Path(row["report"]["path"]).read_text()) for row in panel["admissions"]]
    replay_ids = [g for report in reports for g in report["query_ids"]]
    prefix_audit = read_frozen(Path(partition["inputs"][0]["path"]), partition["inputs"][0]["sha256"])
    fasta = next(Path(r["path"]) for r in prefix_audit["inputs"] if Path(r["path"]).name == "all.fa")
    genes = [r.id for r in SeqIO.parse(fasta, "fasta")]
    query_table = next(Path(r["path"]) for r in partition["outputs"] if Path(r["path"]).name == "queries.tsv")
    replay_from_partition = []
    with query_table.open() as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        count = 0
        for count, row in enumerate(rows, 1):
            if count > len(genes) or row["query"] != genes[count-1] or int(row["input_ordinal_0based"]) != count-1:
                raise ValueError("Original query partition/order differs")
            if row["disposition"] in {"replay_absent", "replay_final_incomplete"}:
                replay_from_partition.append(row["query"])
            elif row["disposition"] != "candidate_prefix_not_admitted":
                raise ValueError("Unknown partition disposition")
    if count != 984137 or len(set(genes)) != count or replay_from_partition != replay_ids:
        raise ValueError("Incomplete or reordered query universe")
    original_log = next(Path(r["path"]) for r in prefix_audit["inputs"] if Path(r["path"]).name == "blast.log")
    selected = select_diagnostics(parse_diagnostics(original_log), panel["coverage"]["diagnostics"], replay_ids)
    replay_set = set(replay_ids)
    raw = {g: lines for g, lines in raw_diagnostics(original_log).items() if g not in replay_set}
    sources = {"prefix": Path(recheck["partial_path"])}
    new_blocks = []
    for index, report in enumerate(reports):
        hits = [Path(r["path"]) for r in report["records"] if Path(r["path"]).name == "hits.blast"]
        if len(hits) != 1:
            raise ValueError("Ambiguous admitted batch table")
        key = f"batch_{index:02d}"
        sources[key] = hits[0]
        new_blocks.extend(dict(block, path=key) for block in report["query_blocks"])
        diagnostic_rows = raw_diagnostics(hits[0].with_name("blast.log"))
        if set(raw) & set(diagnostic_rows):
            raise ValueError("Duplicated selected diagnostics")
        raw.update(diagnostic_rows)
    if set(raw) != set(selected) or not set(raw) <= set(genes):
        raise ValueError("Raw diagnostic selection differs")
    records = unique_records([*records, *panel["records"], record(output / "replay_panel.json")])
    return dict(genes=genes, replay_ids=replay_ids, blocks=new_blocks, sources=sources,
        prefix_audit=prefix_audit, recheck=recheck, panel=panel, records=records,
        selected_diagnostics=selected, raw_diagnostics=raw, accounting=accounting)


def run(root, output, replacement=False):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOB_NODELIST") != "bizon":
        raise ValueError("Require a scheduled merge job on bizon")
    output.mkdir(exist_ok=False)
    sync_directory(output.parent)
    status = dict(status="checking_merge_prerequisites", started_epoch=time.time(),
        job_id=os.environ["SLURM_JOB_ID"], source=record(__file__),
        search_admitted=False, reuse_authorized=False, publication_ready=False)
    save_status(output / "status.json", status)
    try:
        prepared = prepare(root, output, replacement)
        save_status(output / "preflight.json", dict(checked_inputs=prepared["records"],
            accounting=prepared["accounting"], job_id=status["job_id"],
            total_queries=len(prepared["genes"]), replay_queries=len(prepared["replay_ids"]),
            sources={key: str(path) for key, path in prepared["sources"].items()},
            search_admitted=False, publication_ready=False))
        recheck = prepared["recheck"]["result"]
        with Path(prepared["prefix_audit"]["blocks"]["path"]).open() as stream:
            prefix = (dict(block, path="prefix") for line in stream
                      if not (block := json.loads(line))["final_observed_query"])
            plan = ordered_blocks(prepared["genes"], prefix, prepared["blocks"],
                                  prepared["replay_ids"], recheck["retained_end"])
            candidate = merge(plan, prepared["sources"], output / "table")
        if candidate["rows"] != recheck["retained_rows"] + prepared["panel"]["coverage"]["hsp_rows"]:
            raise ValueError("Merged HSP count differs from admitted inventories")
        log = output / "selected.blast.log"
        with log.open("x") as stream:
            for gene in prepared["genes"]:
                stream.writelines(prepared["raw_diagnostics"].get(gene, []))
            stream.flush()
            os.fsync(stream.fileno())
        actual = parse_diagnostics(log)
        if {g: diagnostic_signature(r) for g, r in actual.items()} != {
                g: diagnostic_signature(r) for g, r in prepared["selected_diagnostics"].items()}:
            raise ValueError("Selected diagnostic log differs")
        for item in prepared["records"]:
            check(item)
        status.update(status="merged_candidate_pending_full_table_admission", candidate=candidate,
            selected_log=record(log), diagnostics=actual, checked_inputs=prepared["records"],
            accounting=prepared["accounting"], replay_coverage=prepared["panel"]["coverage"],
            prefix_reuse_review=record(root / "benchmark_tools/results/QFO_BLAST_PREFIX_REUSE_REVIEW_20260923.md"),
            limitations=["No final numerical/alignment audit or whole-search admission yet.",
                         "Conditional prefix reuse retains historical durability/source-build uncertainty.",
                         "Shared-host recovery cost, not successful uninterrupted or controlled timing."])
    except BaseException as error:
        status.update(status="merge_failed_preserved", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        status["finished_epoch"] = time.time()
        save_status(output / "status.json", status)
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replacement", action="store_true")
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute(), args.replacement)
