"""Audit every recovered HSP and query block; not whole-search admission."""

import argparse
import json
from pathlib import Path

from Bio import SeqIO
from benchmark_tools.admit_blast_recovery_batch import index_blocks
from benchmark_tools.audit_orthomcl_search_table import audit_table
from benchmark_tools.check_blast_recovery_dispositions import reconcile
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def expected_blocks(prefix_path, batches):
    retained, boundary, end = [], None, 0
    with prefix_path.open() as stream:
        for line in stream:
            block = json.loads(line)
            if boundary is not None:
                raise ValueError("Inventory continues after incomplete boundary query")
            if block["start"] != end or block["end"] <= end:
                raise ValueError("Noncontiguous original prefix inventory")
            end = block["end"]
            if block["final_observed_query"] is True:
                boundary = block["query"]
            elif block["final_observed_query"] is False:
                retained.append(block)
            else:
                raise ValueError("Invalid prefix boundary flag")
    if boundary is None:
        raise ValueError("Missing excluded boundary query")
    replay_ids = [g for batch in batches for g in batch["query_ids"]]
    if replay_ids.count(boundary) != 1:
        raise ValueError("Excluded boundary query must be replayed exactly once")
    blocks = retained + [block for batch in batches for block in batch["query_blocks"]]
    expected = {}
    for block in blocks:
        if block["query"] in expected or type(block["rows"]) is not int or block["rows"] <= 0:
            raise ValueError("Repeated block or invalid row count")
        expected[block["query"]] = (block["rows"], block["sha256"])
    return retained, expected


def audit(blast, fasta, log, prefix, batch_paths, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    helpers = [Path(__file__).with_name(name) for name in (
        "admit_blast_recovery_batch.py", "audit_orthomcl_search_table.py",
        "audit_orthomcl_blast.py", "convert_orthomcl_blast.py",
        "check_blast_recovery_dispositions.py", "prepare_ob_candidate_neighborhood.py")]
    checked = [record(p) for p in [blast, fasta, log, prefix, *batch_paths, Path(__file__), *helpers]]
    batches = [json.loads(path.read_text()) for path in batch_paths]
    if not batches or any(b["status"] != "recovery_batch_execution_and_rows_verified"
                          or b["batch_admitted"] is not True for b in batches):
        raise ValueError("Require batch validation reports")
    genes = [entry.id for entry in SeqIO.parse(fasta, "fasta")]
    retained, expected = expected_blocks(prefix, batches)
    # This pass independently checks order, complete lines and every block hash;
    # the separate table auditor then checks all numeric/alignment semantics.
    observed = index_blocks(blast, genes)
    if {b["query"]: (b["rows"], b["sha256"]) for b in observed} != expected:
        raise ValueError("Candidate query bytes differ from selected inventories")
    content = audit_table(blast, fasta, log)
    dispositions = reconcile(genes, [b["query"] for b in retained],
        sum(b["rows"] for b in retained), batches, content)
    for item in checked:
        check(item)
    result = dict(status="recovery_candidate_content_verified_not_admitted",
        checked_records=checked, content=content, dispositions=dispositions,
        query_blocks=len(observed), search_admitted=False, reuse_authorized=False,
        publication_ready=False, downstream_execution_authorized=False,
        limitations=[
            "Supplied inventory provenance and native scheduler execution require separate admission.",
            "This is not formatted-database sequence parity or authorization for historical prefix reuse.",
            "No final BPO conversion, clustering, failure-impact analysis or accuracy assessment is admitted."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("blast", "fasta", "log", "prefix", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--batch", type=Path, action="append", required=True)
    args = parser.parse_args()
    audit(args.blast.resolve(), args.fasta.resolve(), args.log.resolve(),
          args.prefix.resolve(), [p.resolve() for p in args.batch], args.output.absolute())
