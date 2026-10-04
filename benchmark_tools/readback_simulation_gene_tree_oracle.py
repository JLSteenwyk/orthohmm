"""Independent stdlib count/decomposition readback; does not rerun reconciliation."""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


ADMISSION_SHA = "17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820"
SCORES_SHA = "1f965d18f572864922e1cdaffb9e4f0443d4ba91ab8b7418dbc99482a162674b"
ARMS = ("inferred", "generating_root", "generating_rerooted")
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control"}


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def load(path, digest):
    ref = record(path)
    if ref["sha256"] != digest:
        raise ValueError("Changed readback input: " + str(path))
    return json.loads(Path(path).read_text()), ref


def validate_counts(score):
    for name in ("tp", "fp", "fn", "predicted_pairs", "eligible_true_pairs"):
        if type(score[name]) is not int or score[name] < 0:
            raise ValueError("Invalid nonnegative integer count")
    tp, fp, fn = [score[k] for k in ("tp", "fp", "fn")]
    if score["predicted_pairs"] != tp + fp or score["eligible_true_pairs"] != tp + fn:
        raise ValueError("Pair-count identity mismatch")
    for name, numerator, denominator in (("f1", 2 * tp, 2 * tp + fp + fn),
        ("precision", tp, tp + fp), ("recall", tp, tp + fn)):
        value = numerator / denominator if denominator else 0.0
        if not math.isclose(score[name], value, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Metric differs from counts")


def retained_baselines(scores):
    rows = [r for r in scores["records"]
            if r["method"] == "orthohmm_satellite_v2" and r["arm"] == "generating"]
    expected = {(c, s) for c in CONDITIONS for s in range(20261101, 20261111)}
    if (len(rows) != 70 or {(r["condition"], r["seed"]) for r in rows} != expected
            or any(r["status"] != "complete" for r in rows)):
        raise ValueError("Incomplete retained baseline inventory or completion status")
    return {(r["condition"], r["seed"]): r["score"] for r in rows}


def decompose(cell, groups, ancestors, true_pairs):
    seen, index = set(), {}
    if len(cell["candidates"]) != len(groups):
        raise ValueError("Candidate inventory mismatch")
    by_name = {r["family"]: r for r in cell["candidates"]}
    if len(by_name) != len(groups):
        raise ValueError("Duplicate candidate result")
    for family, genes in groups.items():
        if not genes or seen & genes or not genes <= ancestors.keys():
            raise ValueError("Invalid complete candidate partition")
        seen.update(genes)
        index.update({g: family for g in genes})
        row = by_name[family]
        ancestry = sorted({ancestors[g] for g in genes})
        if row["genes"] != len(genes) or row["ancestral_families"] != ancestry:
            raise ValueError("Candidate size or ancestry mismatch")
        if row["status"] == "oracle_eligible" and len(ancestry) != 1:
            raise ValueError("Mixed candidate received oracle intervention")
        if row["status"] not in {"oracle_eligible", "mixed_ancestry_ineligible", "unambiguous_bypass"}:
            raise ValueError("Unknown candidate status")
        if row["status"] != "oracle_eligible" and any(row["arms"][arm] != row["arms"]["inferred"] for arm in ARMS):
            raise ValueError("Ineligible or bypass predictions changed")
        local_true = sum(set(pair) <= genes for pair in true_pairs)
        for arm in ARMS:
            validate_counts(row["arms"][arm])
            if row["arms"][arm]["eligible_true_pairs"] != local_true:
                raise ValueError("Local truth-pair count mismatch")
    if seen != set(ancestors):
        raise ValueError("Incomplete candidate universe")
    unavailable = sum(index[a] != index[b] for a, b in true_pairs)
    result = {"true_pairs_across_candidates": unavailable, "candidate_counts": dict(Counter(
        r["status"] for r in cell["candidates"])), "residual_by_arm": {}}
    for arm in ARMS:
        score = cell["arms"][arm]
        validate_counts(score)
        summed = {k: sum(r["arms"][arm][k] for r in cell["candidates"]) for k in ("tp", "fp", "fn")}
        if (summed["tp"] != score["tp"] or summed["fp"] != score["fp"]
                or summed["fn"] + unavailable != score["fn"]):
            raise ValueError("Whole-dataset/candidate count decomposition mismatch")
        result["residual_by_arm"][arm] = {
            status: {k: sum(r["arms"][arm][k] for r in cell["candidates"] if r["status"] == status)
                     for k in ("tp", "fp", "fn")}
            for status in ("oracle_eligible", "mixed_ancestry_ineligible", "unambiguous_bypass")}
    return result


def run(repo, path, digest):
    detailed, detailed_ref = load(path, digest)
    admission, admission_ref = load(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json", ADMISSION_SHA)
    scores, scores_ref = load(repo / "benchmark_tools/results/simulation_tree_robustness_summary_20260917.json", SCORES_SHA)
    inputs = {}
    for ref in detailed["inputs"]:
        if ref["path"] in inputs or record(ref["path"]) != ref:
            raise ValueError("Duplicate or changed retained input")
        inputs[ref["path"]] = ref
    old = retained_baselines(scores)
    native = {(r["condition"], r["seed"]): r for r in admission["records"]
              if r["method"] == "orthohmm_satellite_v2" and r["variant"] == "generating"}
    expected = {(c, s) for c in CONDITIONS for s in range(20261101, 20261111)}
    cells = detailed["cells"]
    if len(cells) != 70 or {(r["condition"], r["seed"]) for r in cells} != expected or set(old) != expected or set(native) != expected:
        raise ValueError("Incomplete fixed panel or retained baseline")
    checked_cells, topology = [], Counter()
    for cell in cells:
        key = (cell["condition"], cell["seed"])
        if cell["status"] != "baseline_reproduced_oracle_scored" or cell["arms"]["inferred"] != old[key]:
            raise ValueError("Baseline differs from independently retained original score")
        truth_ref = native[key]["truth"]
        truth_path = truth_ref.get("absolute_path") or truth_ref["path"]
        if inputs[str(Path(truth_path).resolve())]["sha256"] != truth_ref["sha256"]:
            raise ValueError("Truth identity differs from native admission")
        truth = json.loads(Path(truth_path).read_text())
        ancestors = {g: f for f, genes in truth["families"].items() for g in genes}
        if len(ancestors) != sum(map(len, truth["families"].values())):
            raise ValueError("Repeated ancestral-family membership")
        pairs = {tuple(sorted(p)) for p in truth["ortholog_pairs"]}
        if len(pairs) != len(truth["ortholog_pairs"]):
            raise ValueError("Duplicate truth pairs")
        directory = Path(native[key]["native_report"]["path"]).parent / "orthohmm_satellite_v2"
        candidate_path = directory / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
        if str(candidate_path) not in inputs:
            raise ValueError("Unbound candidate input")
        groups = {}
        for i, line in enumerate(candidate_path.read_text().splitlines()):
            genes = line.split()
            if len(genes) != len(set(genes)):
                raise ValueError("Repeated candidate gene")
            groups[f"Family{i:07d}"] = set(genes)
        decomposition = decompose(cell, groups, ancestors, pairs)
        for candidate in cell["candidates"]:
            if candidate["status"] == "oracle_eligible":
                topology["eligible"] += 1
                topology["rooted_disagreement"] += candidate["rooted_clade_distance"] > 0
                topology["unrooted_disagreement"] += candidate["unrooted_split_distance"] > 0
                topology["root_only_disagreement"] += candidate["rooted_clade_distance"] > 0 and candidate["unrooted_split_distance"] == 0
        checked_cells.append({"label": cell["label"], "condition": cell["condition"], "seed": cell["seed"],
                              "baseline_matches_retained_score": True, "arms": cell["arms"], **decomposition})
    summaries = []
    for condition in sorted(CONDITIONS):
        rows = [r for r in checked_cells if r["condition"] == condition]
        means = {a: {m: sum(r["arms"][a][m] for r in rows) / len(rows)
                      for m in ("f1", "precision", "recall")} for a in ARMS}
        prior = next(r for r in detailed["summary"] if r["condition"] == condition)
        if prior["mean_metrics"] != means or prior["scored_cells"] != len(rows):
            raise ValueError("Finite-panel mean mismatch")
        residual = {a: {status: {k: sum(r["residual_by_arm"][a][status][k] for r in rows) for k in ("tp", "fp", "fn")}
                      for status in ("oracle_eligible", "mixed_ancestry_ineligible", "unambiguous_bypass")} for a in ARMS}
        summaries.append({"condition": condition, "scored_cells": len(rows), "mean_metrics": means,
            "true_pairs_across_candidates": sum(r["true_pairs_across_candidates"] for r in rows),
            "oracle_total_fn": sum(r["arms"]["generating_root"]["fn"] for r in rows),
            "candidate_counts": dict(sum((Counter(r["candidate_counts"]) for r in rows), Counter())),
            "residual_by_arm": residual})
    return {"status": "independent_count_and_candidate_decomposition_verified", "source": record(__file__),
            "detailed_report": detailed_ref, "native_admission": admission_ref, "retained_scores": scores_ref,
            "input_identities_rechecked": len(inputs), "cells": checked_cells, "summary": summaries,
            "topology": dict(topology), "independent_confirmation": False, "publication_ready": False,
            "limitations": ["Readback independently verifies counts, truth partitions and existing baseline scores, not new oracle pair sets or tree inference.",
                "Descriptive finite-panel means and residual counts, not population intervals or significance.",
                "Cross-candidate true pairs cannot be recovered by this fixed-candidate gene-tree intervention; this does not isolate HMM search from clustering."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.repo.resolve(), args.report, args.report_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
