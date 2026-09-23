"""Prepare frozen full-database recovery query batches, without running BLAST."""

import argparse
import json
from pathlib import Path

from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

PARTITION_SHA = "c37381bef1cdb406849f62accde679480a4abe76193cc0fdc65a90e16c8d6a34"
COMPARISON_SHA = "9d49960fd22aaa1b538dcae728b06d0a67802d8d3d1f4296b31c459a98f702e1"


def batch_members(ids):
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Empty or duplicate recovery query inventory")
    return [ids[start:start+5000] for start in range(0, len(ids), 5000)]


def prepare(root, output, protocol_sha256):
    root, output = root.resolve(), output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    partition_path = results / "qfo_blast_recovery_partition_20260923.json"
    comparison_path = results / "qfo_blast_replay_comparison_22055.json"
    protocol = record(results / "QFO_BLAST_RECOVERY_PROTOCOL_20260923.md")
    if protocol["sha256"] != protocol_sha256:
        raise ValueError("Recovery protocol changed")
    partition = read_frozen(partition_path, PARTITION_SHA)
    comparison = read_frozen(comparison_path, COMPARISON_SHA)
    if (comparison["all_diagnostics_compatible"] is not True or comparison["reuse_authorized"] is not False
            or partition["replay_queries"] != 98913 or partition["replay_executed"] is not False):
        raise ValueError("Unexpected recovery evidence")
    source = next(r for r in partition["outputs"] if Path(r["path"]).name == "replay.fa")
    checked = [record(partition_path), record(comparison_path), protocol, source, record(__file__)]
    for item in checked:
        check(item)
    index = SeqIO.index(source["path"], "fasta")
    try:
        ids = list(index)
        batches = batch_members(ids)
        if len(ids) != 98913 or len(batches) != 20 or len(batches[-1]) != 3913:
            raise ValueError("Wrong frozen query-batch panel")
        output.mkdir(parents=True)
        rows = []
        for number, members in enumerate(batches):
            path = output / f"queries_{number:02d}.fa"
            with path.open("xb") as stream:
                for gene in members:
                    stream.write(index.get_raw(gene))
            rows.append(dict(index=number, queries=len(members), first_query=members[0], last_query=members[-1],
                input=record(path), replay_ordinal_start=number*5000,
                replay_ordinal_end_exclusive=number*5000+len(members)))
    finally:
        index.close()
    for item in checked:
        check(item)
    result = dict(status="blast_recovery_batches_prepared_not_executed", inputs=checked, batches=rows,
        total_queries=len(ids), batch_size=5000, replay_executed=False, search_admitted=False,
        reuse_authorized=False, publication_ready=False,
        limitations=["Query inputs only; no database reformatted, search launched or prefix admitted.",
            "Every batch must search the unchanged complete database."])
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root, args.output, args.protocol_sha256)
