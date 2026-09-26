"""Freeze label-free numeric inputs for the prespecified matched-search graph arms."""

import argparse
import csv
import json
import math
from pathlib import Path

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.score_search_sensitivity import read_hits


RESULT_SHA = "d58b1603e563eda27d90cc113950fb28091c94d57233a60ec27b32e297ad729c"


def convert(paths, genes, arm):
    if arm not in {"hmm", "diamond"}:
        raise ValueError("Unknown search arm")
    names = sorted(genes)
    indices = {g: i for i, g in enumerate(names)}
    species = {s: i for i, s in enumerate(sorted({v[0] for v in genes.values()}))}
    hits = {}
    for path, target in paths:
        read_hits(path, genes, arm, target)
        with path.open() as stream:
            reader = csv.reader(stream, delimiter="\t")
            if arm == "hmm":
                next(reader)
            for row in reader:
                if arm == "hmm":
                    q, t, value = row[2], row[3], float(row[4])
                else:
                    if float(row[6]) > 1e-40:
                        continue
                    q, t = row[:2]
                    value = float(row[4]) / math.sqrt(genes[q][1] * genes[t][1])
                key = indices[q], indices[t]
                if key in hits:
                    raise ValueError("Duplicate directed hit across files")
                if not math.isfinite(value) or value <= 0:
                    raise ValueError("Invalid graph score")
                hits[key] = value
    ordered = sorted(hits)
    return dict(gene_names=names, gene_to_species=[species[genes[g][0]] for g in names],
                hit_queries=[q for q, t in ordered], hit_targets=[t for q, t in ordered],
                hit_scores=[hits[pair] for pair in ordered])


def prepare(report_path, protocol, output):
    report_record = record(report_path)
    if report_record["sha256"] != RESULT_SHA:
        raise ValueError("Unrecognized matching result")
    report = json.loads(report_path.read_text())
    if not report["reporting_match_gate_passed"] or report["selected_cutoff"] != 1e-40:
        raise ValueError("Matching gate/cutoff mismatch")
    check(report["plan"])
    check(report["submission"])
    submission = json.loads(Path(report["submission"]["path"]).read_text())
    plan = json.loads(Path(report["plan"]["path"]).read_text())
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    checked = [report_record, record(protocol), report["plan"], report["submission"]]
    for index, dataset in enumerate(plan["datasets"]):
        if dataset["split"] != "reporting":
            continue
        cell = Path(submission["output"]) / f"cell_{index}"
        receipt_record = next(r for r in report["receipts"] if Path(r["path"]).parent == cell)
        check(receipt_record)
        receipt = json.loads(Path(receipt_record["path"]).read_text())
        checked.append(receipt_record)
        genes = {}
        for item in dataset["inputs"]:
            check(item)
            checked.append(item)
            for seq in SeqIO.parse(item["path"], "fasta"):
                if seq.id in genes or not len(seq):
                    raise ValueError("Repeated/empty input gene")
                genes[seq.id] = (Path(item["path"]).name, len(seq))
        for item in receipt["outputs"]:
            check(item)
        for arm in ("hmm", "diamond"):
            paths = ([(cell / "hmm/hits.tsv", None)] if arm == "hmm" else
                     [(cell / "diamond" / (Path(r["path"]).stem + ".tsv"), Path(r["path"]).name)
                      for r in dataset["inputs"]])
            source_records = [record(p) for p, _ in paths]
            checked.extend(source_records)
            numeric = convert(paths, genes, arm)
            destination = output / f"cell_{index}_{arm}.json"
            with destination.open("x") as stream:
                json.dump(numeric, stream, separators=(",", ":"))
                stream.write("\n")
            rows.append(dict(dataset_index=index, condition=dataset["condition"], seed=dataset["seed"],
                             arm=arm, numeric=record(destination), sources=source_records,
                             genes=len(genes), hits=len(numeric["hit_scores"])))
    if len(rows) != 70 or {r["seed"] for r in rows} != set(range(20261106, 20261111)):
        raise ValueError("Expected two arms for all 35 reporting datasets")
    for item in checked:
        check(item)
    result = dict(status="numeric_inputs_prepared_no_graph_inference", source=record(__file__),
                  protocol=record(protocol), search_result=report_record, checked_records=checked,
                  cells=rows, graph_settings=dict(cpm_resolution=.1, leiden_seed=4,
                      include_isolates=True, refinement=True, profile_expansion=False,
                      candidate_expansion=False, phylogeny=False),
                  normalization=dict(hmm="native score unchanged", diamond="raw / sqrt(query_length * target_length)"),
                  independent_validation=False, scientific_defaults_changed=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.report.resolve(), args.protocol.resolve(), args.output.absolute())
