"""Localize retained VGNC misses without assigning causes to absent search hits."""

import argparse
from collections import Counter
import csv
import gzip
import json
import math
from pathlib import Path
import resource
import sqlite3
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.audit_accuracy_checkpoint import check_arrays
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.audit_qfo_vgnc_mapping import reference_data, mapped_reference, validate_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import argument
from benchmark_tools.score_ygob_groups import membership, read_predictions
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

PROTOCOL = "hillclimb_qfo_failure_localization_20261012.json"


def normalized_index(names):
    index = {}
    for i, gene in enumerate(names):
        accession = _strip_to_uniprot(gene)
        if accession in index:
            raise ValueError("Ambiguous normalized accession")
        index[accession] = i
    return index


def accession_truth(database, truth, labels):
    metadata = {}
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        for protein, accession in db.execute("SELECT prot_nr, uniprot_id FROM proteomes ORDER BY rowid"):
            if protein in labels:
                metadata[protein] = accession
    if set(metadata) != set(labels) or len(set(metadata.values())) != len(labels):
        raise ValueError("Incomplete or non-bijective selected reference mapping")
    result = {}
    for (a, b), family in truth.items():
        pair = tuple(sorted((metadata[a], metadata[b])))
        if pair in result:
            raise ValueError("Reference normalization collapses asserted pairs")
        result[pair] = family
    return result


def watched_hits(names, species, queries, targets, scores, wanted, chunk_size=1000000):
    summary = check_arrays(names, species, queries, targets, scores, chunk_size)
    if summary["nonpositive_scores"]:
        raise ValueError("Require positive normalized initial hit scores")
    watched = np.zeros(len(names), dtype=bool)
    watched[list({i for pair in wanted for i in pair})] = True
    found = {}
    for start in range(0, len(scores), chunk_size):
        q, t, s = (values[start:start + chunk_size] for values in (queries, targets, scores))
        selected = watched[q] & watched[t] & (q != t)
        for a, b, score in zip(q[selected], t[selected], s[selected]):
            a, b = int(a), int(b)
            if tuple(sorted((a, b))) not in wanted:
                continue
            if (a, b) in found:
                raise ValueError("Duplicate ordered watched hit")
            found[a, b] = float(score)
    return found, summary


def trace_rows(truth, predicted, names, index, species, candidates, roots, hits):
    rows = []
    for (a, b), family in sorted(truth.items()):
        i, j = index[a], index[b]
        left, right = names[i], names[j]
        if species[i] == species[j]:
            raise ValueError("VGNC reference pair is not cross-species")
        same_candidate = candidates[left] == candidates[right]
        retained = (a, b) in predicted
        if retained and not same_candidate:
            raise ValueError("Native predicted pair crosses candidate families")
        stage = ("predicted" if retained else "within_candidate_not_predicted" if same_candidate
                 else "different_candidate_families")
        forward, reverse = hits.get((i, j)), hits.get((j, i))
        rows.append({"left": a, "right": b, "reference_label": family,
            "category": "TP" if retained else "FN", "stage": stage,
            "initial_hit_directions": int(forward is not None) + int(reverse is not None),
            "forward_normalized_score": forward, "reverse_normalized_score": reverse,
            "same_root_hog": roots[left] == roots[right]})
    return rows


def native_truth_pairs(path, names, wanted, expected_count):
    index = normalized_index(names)
    selected, previous, count = set(), None, 0
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["gene_a", "species_a", "gene_b", "species_b"]:
            raise ValueError("Unexpected native pair header")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed native pair row")
            a, b = row["gene_a"], row["gene_b"]
            if a >= b or (previous is not None and (a, b) <= previous):
                raise ValueError("Noncanonical or duplicate native pair")
            if row["species_a"] == row["species_b"]:
                raise ValueError("Native pair is not cross-species")
            previous, count = (a, b), count + 1
            a, b = _strip_to_uniprot(a), _strip_to_uniprot(b)
            if a not in index or b not in index:
                raise ValueError("Native gene absent from initial checkpoint")
            pair = tuple(sorted((a, b)))
            if pair in wanted:
                selected.add(pair)
    if count != expected_count:
        raise ValueError("Native pair count differs from admitted baseline")
    return selected


def prepare(repo, comparison, output):
    if output.exists():
        raise FileExistsError(output)
    report = json.loads(comparison.read_text())
    method = next(m for m in report["methods"] if m["key"] == "membership:control_direct")
    cell = next(c for c in method["cells"] if c["dataset"] == "QfO")
    if cell["status"] != "independent_content_readback_verified":
        raise ValueError("Need independently verified baseline")
    factorial = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    original = json.loads(factorial.read_text())
    checkpoint = original["numeric_checkpoint"]["manifest"]
    check(checkpoint)
    manifest = json.loads(Path(checkpoint["path"]).read_text())
    candidates = original["candidate_arms"]["p1_c1"]["candidate_partition"]
    command = cell["resources"]["resources"]["command"]
    if argument(command, "--candidate-clusters") != candidates["path"] or candidates not in cell["checked"]:
        raise ValueError("Candidate ancestry differs from verified baseline")
    folder = Path(cell["products"]["orthohmm_root_hogs.tsv"]["path"]).parents[2]
    reference_inventory = repo / "benchmark_tools/results/vgnc_reference_inventory_20260917.json"
    reference = json.loads(reference_inventory.read_text())["reference"]
    inputs = {"comparison": record(comparison), "factorial": record(factorial),
        "reference_inventory": record(reference_inventory), "reference": reference,
        "checkpoint_manifest": checkpoint, "candidates": candidates,
        "roots": cell["products"]["orthohmm_root_hogs.tsv"],
        "pairs": cell["products"]["orthohmm_pairwise_orthologs.tsv"],
        "raw": record(folder / "qfo_scoring/results/VGNC/VGNC_ohmm-hillclimb-r1-0_raw.txt.gz"),
        "database": record(folder / "qfo_scoring/other/ohmm_hillclimb_r1_0.db"),
        "protocol": record(repo / "benchmark_tools" / PROTOCOL)}
    for name, value in manifest["files"].items():
        inputs["checkpoint_" + name] = {"path": str(Path(checkpoint["path"]).parent / name), **value}
    for value in inputs.values():
        check(value)
    sources = [record(Path(module.__file__)) for module in
        (sys.modules[__name__], sys.modules["orthohmm.accuracy"],
         sys.modules["benchmark_tools.audit_accuracy_checkpoint"],
         sys.modules["benchmark_tools.audit_historical_profile_ablation"],
         sys.modules["benchmark_tools.audit_qfo_vgnc_mapping"],
         sys.modules["benchmark_tools.prepare_ob_candidate_neighborhood"],
         sys.modules["benchmark_tools.run_membership_hillclimb"],
         sys.modules["benchmark_tools.score_ygob_groups"], sys.modules["qfo_benchmark.og_to_pairwise"])]
    plan = {"status": "prepared_unrun", "inputs": inputs, "sources": sources,
        "source_revision": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "expected_universe": original["numeric_checkpoint"]["summary"],
        "expected_pairs": cell["coverage"]["pairs"], "expected_metrics": cell["details"]["VGNC F"],
        "expected_f1": cell["claimed_scores"]["VGNC F"]}
    with output.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return plan


def run(plan_path, output):
    if output.exists():
        raise FileExistsError(output)
    plan_record = record(plan_path)
    plan = json.loads(plan_path.read_text())
    checked = [plan_record, *plan["inputs"].values(), *plan["sources"]]
    for value in checked:
        check(value)
    started, cpu = time.monotonic(), time.process_time()
    inputs = {key: Path(value["path"]) for key, value in plan["inputs"].items()}
    names, species, queries, targets, scores = load_accuracy_checkpoint(inputs["checkpoint_manifest"].parent)
    index = normalized_index(names)
    truth, labels = reference_data(inputs["reference"])
    asserted, annotations, mapping_digest, aliases = mapped_reference(inputs["database"], truth, labels)
    raw = validate_raw(inputs["raw"], asserted, annotations)
    tp, fp, fn = (raw["counts"][label] for label in ("TP", "FP", "FN"))
    for metric, value, expected in (("precision", tp / (tp + fp), plan["expected_metrics"]["precision"]),
                                    ("recall", tp / (tp + fn), plan["expected_metrics"]["recall"]),
                                    ("f1", 2 * tp / (2 * tp + fp + fn), plan["expected_f1"])):
        if not math.isfinite(value) or not math.isclose(value, expected, rel_tol=0, abs_tol=5e-8):
            raise ValueError("Raw VGNC metric differs from verified baseline: " + metric)
    missing = set(annotations) - index.keys()
    if missing:
        raise ValueError(f"Reference accession missing from checkpoint: {len(missing)}")
    wanted = {tuple(sorted((index[a], index[b]))) for a, b in asserted}
    hits, summary = watched_hits(names, species, queries, targets, scores, wanted)
    if summary != plan["expected_universe"]:
        raise ValueError("Initial checkpoint universe differs from admitted baseline")
    universe = set(names)
    candidates = {g: i for i, group in enumerate(read_partition(inputs["candidates"], universe)) for g in group}
    roots = membership(read_predictions(inputs["roots"], "root_hogs"))
    if set(roots) != universe:
        raise ValueError("Root-HOG gene universe differs from checkpoint")
    predicted = native_truth_pairs(inputs["pairs"], names, asserted, plan["expected_pairs"])
    # Native pair output, not RootHOG co-membership, determines true positives.
    observed = {"TP": set(), "FN": set()}
    with gzip.open(inputs["raw"], "rt") as stream:
        for line in stream:
            a, b, category, *_ = line.rstrip("\n").split("\t")
            if category in observed:
                observed[category].add(tuple(sorted((a, b))))
    if predicted != observed["TP"] or asserted - predicted != observed["FN"]:
        raise ValueError("Native truth pairs differ from raw VGNC TP/FN")
    # Map each asserted row's own label, retaining multi-label reference proteins.
    mapped_truth = accession_truth(inputs["database"], truth, labels)
    if set(mapped_truth) != asserted:
        raise ValueError("Selected reference mappings disagree")
    rows = trace_rows(mapped_truth, predicted, names, index, species, candidates, roots, hits)
    for value in checked:
        check(value)
    output.mkdir(parents=True, exist_ok=False)
    table = output / "reference_pair_stages.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    result = {"status": "native_reference_pair_stage_trace_verified", "checked": checked,
        "source_revision": plan["source_revision"], "checkpoint": summary,
        "reference_pairs": len(rows), "raw_categories": raw, "mapping_sha256": mapping_digest,
        "mapping_extra_rows": aliases, "table": record(table),
        "summary": dict(Counter(f'{r["category"]}/{r["stage"]}/hits_{r["initial_hit_directions"]}/root_{r["same_root_hog"]}' for r in rows)),
        "resources": {"wall_s": time.monotonic() - started, "cpu_s": time.process_time() - cpu,
                      "process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "defaults_promoted": False, "goal_completed": False,
        "limitations": ["Development-exposed diagnostic, not a new search or independent validation.",
            "Absent significant initial hits do not identify prefilter, cap, ranking, or significance rejection.",
            "Direct hit presence is not native RBNH edge presence or a causal explanation of grouping.",
            "Within-candidate omission combines reconciliation and membership mechanisms; RootHOGs are not pair cliques.",
            "VGNC rows may share proteins and labels; no independent-family uncertainty is asserted.",
            "Duplicate ordered hits are rejected for watched pairs only; this is not a global duplicate-hit audit.",
            "Diagnostic resource receipts are shared-host analysis timings, not inference timings."]}
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    preparation = commands.add_parser("prepare")
    preparation.add_argument("--repo", required=True, type=Path)
    preparation.add_argument("--comparison", required=True, type=Path)
    preparation.add_argument("--output", required=True, type=Path)
    execution = commands.add_parser("run")
    execution.add_argument("--plan", required=True, type=Path)
    execution.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.operation == "prepare":
        prepare(args.repo.resolve(), args.comparison.resolve(), args.output.resolve())
    else:
        result = run(args.plan.resolve(), args.output.resolve())
        print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
