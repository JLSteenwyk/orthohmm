"""Bounded synthetic profile-weight and threshold diagnostic, not orthology F1."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from benchmark_tools.hmm_generalization_holdout import (
    _score_pairs, classification_metrics, generate_synthetic_dataset, species_sequences,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from orthohmm.search.evalue import batch_estimate_evalues
from orthohmm.search.matrices import get_ka_params
from orthohmm.search.profile_expansion import build_cluster_profiles
from orthohmm.search.sequence_weights import PROFILE_SEQUENCE_WEIGHTING


def training_database(dataset, extra_copies):
    names, sequences, species, clusters = [], [], [], []
    for family, original in sorted(dataset.training_sequences.items()):
        members = list(original) + [original[0]] * extra_copies
        labels = list(range(len(original))) + [0] * extra_copies
        cluster = []
        for index, (sequence, label) in enumerate(zip(members, labels)):
            cluster.append(len(names))
            names.append(f"family_{family}_member_{index}")
            sequences.append(sequence)
            species.append(label)
        clusters.append(cluster)
    return names, sequences, species, clusters


def summarize(query_scores, evalues, thresholds, queries):
    """Evaluate only after profile construction and scoring have finished."""
    predictions, top_correct, margins = [], 0, []
    for index, query in enumerate(queries):
        column = query_scores[:, index]
        eligible = (column >= thresholds) & (evalues[:, index] < 1e-4)
        selected = np.flatnonzero(eligible)
        predictions.append(int(selected[np.argmax(column[selected])]) if len(selected) else None)
        if query.family_id is not None:
            top_correct += int(int(np.argmax(column)) == query.family_id)
            margins.append(float(column[query.family_id] - thresholds[query.family_id]))
    return {
        "top1_positive_family_correct": top_correct,
        "strict_assignment": classification_metrics(queries, predictions),
        "own_family_threshold_margin": {
            "minimum": min(margins), "median": float(np.median(margins)), "maximum": max(margins),
        },
    }


def evaluate(dataset, extra_copies, cpu, directory):
    names, sequences, species, clusters = training_database(dataset, extra_copies)
    n_train = len(names)
    database = species_sequences(
        names + [query.query_id for query in dataset.queries],
        sequences + [query.sequence for query in dataset.queries], "weighting_diagnostic",
    )
    pairs = np.asarray([(family, target) for family in range(len(clusters))
                        for target in range(database.num_sequences)], dtype=np.int32)
    lam, ka_k = get_ka_params("BLOSUM62")
    rows = []
    for mode in PROFILE_SEQUENCE_WEIGHTING:
        started = time.perf_counter()
        profiles = build_cluster_profiles(
            clusters, names, database, "BLOSUM62", cpu,
            sequence_weighting=mode, gene_to_species=species,
        )
        if set(profiles) != set(range(len(clusters))):
            raise ValueError("Synthetic profile inventory differs from training families")
        profile_ids, lengths, flat = _score_pairs(profiles, database, pairs, cpu)
        if not np.array_equal(profile_ids, np.arange(len(clusters))):
            raise ValueError("Noncanonical synthetic profile order")
        scores = flat.reshape(len(clusters), database.num_sequences)
        thresholds = np.array([scores[family, cluster].min() for family, cluster in enumerate(clusters)])
        query_scores = scores[:, n_train:]
        evalues = batch_estimate_evalues(
            query_scores.reshape(-1), np.repeat(lengths, len(dataset.queries)),
            int(database.lengths.sum()), lam, ka_k,
        ).reshape(query_scores.shape)
        if not np.isfinite(scores).all() or np.isnan(evalues).any():
            raise ValueError("Invalid synthetic numerical outputs")
        raw = directory / f"seed_{dataset.seed}_copies_{extra_copies}_{mode}.npz"
        np.savez_compressed(raw, scores=scores, thresholds=thresholds,
                            evalues=evalues, profile_lengths=lengths)
        rows.append({
            "mode": mode, "extra_copies_per_family": extra_copies,
            "seed": dataset.seed, "profiles": len(profiles), "training_sequences": n_train,
            "effective_sequence_counts": [profiles[index].effective_sequence_count for index in range(len(clusters))],
            **summarize(query_scores, evalues, thresholds, dataset.queries),
            "raw": record(raw), "wall_s": time.perf_counter() - started,
        })
    training = directory / f"seed_{dataset.seed}_copies_{extra_copies}_generated.json"
    training.write_text(json.dumps({
        "names": names, "sequences": sequences, "species": species, "clusters": clusters,
        "queries": [vars(query) for query in dataset.queries],
        "query_truth_scope": "Used by summarize only, never profile construction or scoring",
    }, indent=2) + "\n")
    return {"rows": rows, "generated": record(training)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=4)
    args = parser.parse_args(argv)
    if args.output.exists() or not 1 <= args.cpu <= 16:
        raise ValueError("Require a fresh output directory and 1-16 CPUs")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Run the prespecified CPU diagnostic with CUDA_VISIBLE_DEVICES empty")
    source = Path(__file__).resolve().parents[1]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "orthohmm",
                    "benchmark_tools/diagnose_profile_weighting.py",
                    "benchmark_tools/hmm_generalization_holdout.py",
                    "benchmark_tools/hmm_evidence_campaign_20261018.json"], cwd=source, check=True)
    plan = source / "benchmark_tools/hmm_evidence_campaign_20261018.json"
    specification = json.loads(plan.read_text())
    if [arm["label"] for arm in specification["first_experiment"]["arms"]] != list(PROFILE_SEQUENCE_WEIGHTING):
        raise ValueError("Weighting arms differ from the prespecified campaign")
    paths = [plan, Path(__file__).resolve(), source / "benchmark_tools/hmm_generalization_holdout.py",
             source / "benchmark_tools/run_search_recovery.py", source / "benchmark_tools/prepare_ob_candidate_neighborhood.py",
             *sorted((source / "orthohmm").rglob("*.py")),
             *sorted((source / "orthohmm/search/csrc").glob("*.so"))]
    pins = [record(path) for path in paths]
    args.output.mkdir(parents=True)
    before = host_activity()
    started = time.perf_counter()
    results = [evaluate(generate_synthetic_dataset(seed, families=36), copies, args.cpu, args.output)
               for seed in (18102026, 19102026) for copies in (0, 4)]
    for pin in pins:
        check(pin)
    report = {
        "schema_version": 1, "status": "completed_synthetic_diagnostic_not_biological_validation",
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip(),
        "sources": pins, "runtime": runtime_records(), "cpu": args.cpu,
        "command": [sys.executable, "-m", "benchmark_tools.diagnose_profile_weighting", *sys.argv[1:]],
        "scope": "All-profile scoring diagnostic, excluding prefilter, anchors and phylogeny; strict family-assignment F1 is not orthology F1. Both seeds are development diagnostics, not independent biological confirmation.",
        "host_before": before, "host_after": host_activity(),
        "wall_s": time.perf_counter() - started, "results": results,
        "next_action": "Run matched canonical complete-proteome OrthoBench downstream control and three weighted pilots; synthetic signs do not override the prespecified biological pilot.",
    }
    output = args.output / "summary.json"
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(record(output)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
