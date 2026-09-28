"""Descriptive fixed-screen strata of retained YGOB pillar sufficient statistics."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_ygob_groups import METRICS, statistics
from benchmark_tools.verify_ygob_overlap_screen import verify
from benchmark_tools.probe_dgx_step_separation import save

METHODS = ("orthohmm_satellite_v2", "orthohmm_high_sensitivity",
           "orthofinder_full", "orthofinder_sequence_only")
PINS = {
    "ygob_frozen_results_20260916.json": "3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef",
    "ygob_homology_screen_20260916.json": "c142f6ecba369c8192eea131de0827eafabd498399ad25c80b7136b33d34a74f",
    "ygob_overlap_admission_20260916.json": "d2a8aa5d2901c15e50a162aae8c53c54df261753d5a4c3999e89e641477f7054",
}


def validate_records(records):
    signature = {}
    for row in records:
        name, genes = row["pillar"], row["genes"]
        if not isinstance(name, str) or not name or name in signature or type(genes) is not int or genes < 1:
            raise ValueError("Invalid or duplicate pillar")
        signature[name] = genes
        for key in ("tp", "fp", "fn"):
            value = row[key]
            scale = 2 if key == "fp" else 1
            if (type(value) not in (int, float) or not math.isfinite(value)
                    or value < 0 or value * scale != math.floor(value * scale)):
                raise ValueError("Invalid pillar counts")
        if row["tp"] + row["fn"] != genes * (genes - 1) // 2:
            raise ValueError("Truth pair count differs")
        covered = row["covered_genes"]
        if type(covered) is not int or not 0 <= covered <= genes or type(row["exact"]) is not bool:
            raise ValueError("Invalid coverage or exact flag")
        if row["tp"] > covered * (covered - 1) // 2:
            raise ValueError("True pairs exceed covered genes")
        if row["exact"] and (covered != genes or row["fp"] or row["fn"]):
            raise ValueError("Exact pillar has errors or missing genes")
    if not signature:
        raise ValueError("Empty pillar universe")
    return signature


def aggregate(records):
    counts = {k: sum(r[k] for r in records) for k in ("tp", "fp", "fn")}
    tp, fp, fn = (counts[k] for k in ("tp", "fp", "fn"))
    genes = sum(r["genes"] for r in records)
    covered = sum(r["covered_genes"] for r in records)
    return dict(counts=counts, metrics=dict(zip(METRICS, statistics([tp, fp, fn]).tolist())),
        defined=dict(zip(METRICS, [2*tp+fp+fn > 0, tp+fp > 0, tp+fn > 0])),
        reference_groups=len(records), reference_genes=genes, truth_pairs=tp+fn,
        zero_truth_pair_pillars=sum(r["genes"] == 1 for r in records),
        covered_reference_genes=covered, reference_gene_coverage=covered/genes if genes else 0.,
        represented_reference_groups=sum(r["covered_genes"] > 0 for r in records),
        exact_reference_groups=sum(r["exact"] for r in records))


def summarize(scores, positive_ids):
    if set(scores) != set(METHODS) or not isinstance(positive_ids, list):
        raise ValueError("Require all four methods and explicit screen IDs")
    if any(not isinstance(p, str) or not p for p in positive_ids) or len(set(positive_ids)) != len(positive_ids):
        raise ValueError("Duplicate or invalid screen labels")
    positive = set(positive_ids)
    signature = None
    output = {"screen_positive": {}, "screen_negative": {}}
    for method in METHODS:
        score = scores[method]
        records = score["records"]
        current = validate_records(records)
        if signature is not None and current != signature:
            raise ValueError("Methods have different reference signatures")
        signature = current
        if not positive <= set(current):
            raise ValueError("Screen includes unknown pillars")
        full = aggregate(records)
        for key in ("counts", "reference_groups", "reference_genes", "covered_reference_genes",
                    "reference_gene_coverage", "represented_reference_groups", "exact_reference_groups"):
            if full[key] != score[key]:
                raise ValueError("Pillar statistics do not recover full totals: " + key)
        if any(not math.isclose(full["metrics"][k], score["metrics"][k], rel_tol=0, abs_tol=1e-12) for k in METRICS):
            raise ValueError("Pillar metrics do not recover full score")
        for label, selected in (("screen_positive", positive), ("screen_negative", set(current)-positive)):
            output[label][method] = aggregate([r for r in records if r["pillar"] in selected])
        left, right = (output[k][method] for k in output)
        for key in ("tp", "fp", "fn"):
            if left["counts"][key] + right["counts"][key] != full["counts"][key]:
                raise ValueError("Stratum counts do not add up")
        for key in ("reference_groups", "reference_genes", "covered_reference_genes", "exact_reference_groups"):
            if left[key] + right[key] != full[key]:
                raise ValueError("Stratum totals do not add up")
    for methods in output.values():
        baseline = methods["orthofinder_full"]
        for score in methods.values():
            score["difference_vs_orthofinder_percentage_points"] = {
                key: 100*(score["metrics"][key]-baseline["metrics"][key])
                if score["defined"][key] and baseline["defined"][key] else None for key in METRICS}
    return output


def run(root, output):
    base = root / "benchmark_tools/results"
    inputs = {name: read_frozen(base / name, digest) for name, digest in PINS.items()}
    summary, screen, admission = (inputs[name] for name in PINS)
    full_ref = summary["full_results"]
    full = read_frozen(Path(full_ref["path"]), full_ref["sha256"])
    if full_ref["sha256"] != "5104e000d0c3bb8687103646be8011fbc97f8b05f7eb0a925bbfdd1e76200e83":
        raise ValueError("Full statistics differ from protocol")
    if {m: {k: v for k, v in s.items() if k != "records"} for m, s in full["scores"].items()} != summary["scores"]:
        raise ValueError("Full results differ from admitted summary")
    overlap = verify(root)
    if overlap != admission:
        raise ValueError("Retained overlap admission does not reproduce")
    result = summarize(full["scores"], screen["reference_pillars_with_hit_ids"])
    refs = [record(base / name) for name in PINS] + [full_ref]
    for ref in refs:
        check(ref)
    save(output, dict(status="ygob_overlap_strata_descriptive", strata=result, inputs=refs,
        protocol=record(base / "YGOB_OVERLAP_STRATA_PROTOCOL_20260928.md"),
        protocol_commit="e51ffc54", source=record(__file__),
        helpers=[record(Path(__file__).with_name(n)) for n in ("score_ygob_groups.py", "verify_ygob_overlap_screen.py")],
        overlap_reverified=overlap, publication_ready=False, independent_confirmation=False,
        confidence_intervals_computed=False,
        limitations=["Secondary descriptive partition of existing evaluated data, not independent-family confirmation.",
            "No-hit screen status does not establish absence of remote homology or shared annotation ancestry.",
            "Original allocated cross-pillar false positives are retained; these are not subset-rescored predictions.",
            "Undefined ratios follow the frozen zero convention and are explicitly flagged.",
            "No causal attribution, significance, new method tuning or general superiority claim."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    run(args.root.resolve(), args.output)
