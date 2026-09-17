"""Compare one declared ARM smoke dataset with its retained x86 predictions."""

import argparse
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.simulation_method_outputs import orthohmm_pairs


def partition(path, format, owners):
    groups = read_predictions(path, format)
    if set(membership(groups)) != set(owners):
        raise ValueError("Partition does not cover exact input universe")
    return {tuple(sorted(genes)) for genes in groups.values()}


def pair_set(path, owners):
    rows = list(orthohmm_pairs(path, owners))
    if any(a >= b for a, b in rows) or rows != sorted(set(rows)):
        raise ValueError("Pairs are not canonical, sorted and unique")
    return set(rows)


def compare(left, right):
    return {"x86_count": len(left), "arm_count": len(right), "equal": left == right,
            "x86_only": len(left - right), "arm_only": len(right - left)}


def audit(inputs, reference, arm):
    owners, input_records = {}, []
    paths = sorted(inputs.glob("*.fasta"))
    for path in paths:
        input_records.append(record(path))
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in owners or not seq.seq:
                raise ValueError("Invalid input genes")
            owners[seq.id] = path.stem
    if len(paths) != 8 or len(owners) != 645:
        raise ValueError("Wrong declared missing20_20261101 fixture")
    transferred = {}
    for line in (arm / "input.before.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        name = Path(path).name
        if name in transferred:
            raise ValueError("Duplicate transferred input")
        transferred[name] = digest
    if transferred != {Path(r["path"]).name: r["sha256"] for r in input_records}:
        raise ValueError("ARM input bytes differ from reference")
    files = [record(arm / "input.before.sha256")]
    results = {}
    for method, short in (("orthohmm_high_sensitivity", "high"), ("orthohmm_satellite_v2", "satellite")):
        ref_metrics, arm_metrics = reference / (method + ".json"), arm / (short + ".json")
        metrics = [json.loads(p.read_text()) for p in (ref_metrics, arm_metrics)]
        for m in metrics:
            if m["status"] != "complete" or m["counts"]["genes"] != 645 or m["counts"]["species"] != 8:
                raise ValueError("Incomplete native metrics")
        left, right = reference / method, arm / short
        a, b = left / "orthohmm_orthogroups.txt", right / "orthohmm_orthogroups.txt"
        groups = [partition(p, "named_groups", owners) for p in (a, b)]
        if [len(g) for g in groups] != [m["counts"]["orthogroups"] for m in metrics]:
            raise ValueError("Metric group count differs")
        result = {"orthogroups": compare(*groups)}
        files.extend(record(p) for p in (ref_metrics, arm_metrics, a, b))
        if short == "satellite":
            for name, filename, reader, count_key in (
                ("root_hogs", "orthohmm_root_hogs.tsv", lambda p: partition(p, "root_hogs", owners), "phylogeny_root_hogs"),
                ("ortholog_pairs", "orthohmm_pairwise_orthologs.tsv", lambda p: pair_set(p, owners), "phylogeny_ortholog_pairs"),
            ):
                paths = [directory / "orthohmm_phylogeny" / filename for directory in (left, right)]
                sets = [reader(p) for p in paths]
                if [len(s) for s in sets] != [m["counts"][count_key] for m in metrics]:
                    raise ValueError("Metric phylogeny count differs")
                result[name] = compare(*sets)
                files.extend(record(p) for p in paths)
        results[method] = result
    return {"status": "fixture_outputs_audited", "dataset": "missing20_20261101",
            "input_genes": len(owners), "inputs": input_records, "source": record(Path(__file__)),
            "files": files, "comparisons": results,
            "limitations": ["One development-exposed simulation fixture, not general equivalence or independent accuracy.",
                            "Not a timing comparison or full environment/command admission.",
                            "Gene/species tree equivalence and intermediate search scores are not established here."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("inputs", "reference", "arm", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.inputs.resolve(), args.reference.resolve(), args.arm.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
