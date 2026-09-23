"""Prepare an exhaustive proposed query partition, without authorizing BLAST reuse."""

import argparse
import csv
import json
from pathlib import Path
import re

from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

AUDIT_SHA = "a0d024ff33544b64eafcb7ebe88ba3948ad0cffe8aabcb1e4ef2ce300cc71bc9"


def partition(ids, blocks, boundary):
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError("Empty or duplicate query universe")
    blocks = iter(blocks)
    current = next(blocks, None)
    end, count, final = 0, 0, None
    for ordinal, gene in enumerate(ids):
        row = dict(query=gene, input_ordinal_0based=ordinal, disposition="replay_absent",
                   retained_start=None, retained_end=None)
        if current is not None:
            index = current["input_ordinal_0based"]
            if type(index) is not int or index < ordinal:
                raise ValueError("Nonincreasing query block order")
            if index == ordinal:
                if (final is not None or current["query"] != gene or current["reuse_authorized"] is not False
                        or type(current["rows"]) is not int or current["rows"] < 1
                        or type(current["start"]) is not int or type(current["end"]) is not int
                        or current["start"] != end or current["end"] <= end
                        or type(current["final_observed_query"]) is not bool
                        or re.fullmatch("[0-9a-f]{64}", current["sha256"]) is None):
                    raise ValueError("Invalid query block inventory")
                end, count = current["end"], count + 1
                if current["final_observed_query"]:
                    final = current
                    row["disposition"] = "replay_final_incomplete"
                else:
                    row.update(disposition="candidate_prefix_not_admitted",
                               retained_start=current["start"], retained_end=current["end"])
                current = next(blocks, None)
        yield row
    if (current is not None or final is None or count != boundary["block_count"]
            or final["query"] != boundary["last_query"] or final["start"] != boundary["last_query_start"]
            or end != boundary["prefix_end"]):
        raise ValueError("Incomplete block coverage or inconsistent final boundary")


def prepare(audit_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    audit = read_frozen(audit_path, AUDIT_SHA)
    if (audit["status"] != "observed_prefix_rows_validated_not_admitted"
            or audit["reuse_authorized"] is not False or audit["search_admitted"] is not False):
        raise ValueError("Unexpected prefix-audit scope")
    source = next(r for r in audit["inputs"] if Path(r["path"]).name == "all.fa")
    records = [record(audit_path), source, audit["blocks"], record(__file__)]
    for item in records:
        check(item)
    output.mkdir(parents=True)
    index = SeqIO.index(source["path"], "fasta")
    counts = dict(candidate_prefix_not_admitted=0, replay_absent=0, replay_final_incomplete=0)
    try:
        ids = list(index)
        if len(ids) != 984137:
            raise ValueError("Wrong corrected query universe")
        with Path(audit["blocks"]["path"]).open() as blocks, (output / "queries.tsv").open("x") as table, (output / "replay.fa").open("xb") as replay:
            writer = csv.DictWriter(table, delimiter="\t", lineterminator="\n",
                fieldnames=["query", "input_ordinal_0based", "disposition", "retained_start", "retained_end"])
            writer.writeheader()
            for row in partition(ids, (json.loads(line) for line in blocks), audit["boundary"]):
                counts[row["disposition"]] += 1
                writer.writerow({k: "NA" if v is None else v for k, v in row.items()})
                if row["disposition"] != "candidate_prefix_not_admitted":
                    replay.write(index.get_raw(row["query"]))
    finally:
        index.close()
    if counts != dict(candidate_prefix_not_admitted=885224, replay_absent=98912, replay_final_incomplete=1):
        raise ValueError("Frozen partition counts differ")
    for item in records:
        check(item)
    result = dict(status="proposed_blast_query_partition_prepared", inputs=records,
        total_queries=sum(counts.values()), counts=counts, replay_queries=98913,
        proposed_retained_prefix_end=audit["boundary"]["last_query_start"],
        outputs=[record(p) for p in sorted(output.iterdir())],
        search_admitted=False, reuse_authorized=False, replay_executed=False,
        limitations=["Query partition only; does not authorize reuse or verify native replay equivalence.",
            "Absent rows include unprocessed, no-hit and failed queries; absence is not a failure label.",
            "Every replay query must search the unchanged full database, not the replay subset.",
            "Original partial BLAST bytes are not copied, merged, changed or rehashed by this preparer."])
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.audit.resolve(), args.output.absolute())
