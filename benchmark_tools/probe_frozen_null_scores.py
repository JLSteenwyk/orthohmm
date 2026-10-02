"""Prespecified synthetic null-score diagnostic; never tune the frozen filter."""

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.stats import binomtest

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

REVISION = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
NATIVE_SHA = "f21d902852d727d6270b65ba900db52e7216cd874ee2d2b89fa72b6f16dcd5c0"
REGIMES = ("blosum_background", "uniform", "half_glutamine")
LENGTHS = (50, 150, 400)
BANDS = (0, 64)
THRESHOLDS = (1.0, 0.1, 0.01, 0.001, 0.0001)
SEEDS = tuple(range(10))
PAIRS_PER_SEED = 1000
ENDPOINTS = len(REGIMES) * len(LENGTHS) * len(BANDS) * len(THRESHOLDS)
PACKAGES = dict(numpy="2.2.6", numba="0.65.0", llvmlite="0.47.0", scipy="1.15.3")


def probabilities(background, regime):
    background = np.asarray(background, dtype=float)
    if (background.shape != (20,) or not np.isfinite(background).all()
            or (background <= 0).any()):
        raise ValueError("Require positive finite 20-residue background")
    background = background / background.sum()
    if regime == "blosum_background":
        return background
    if regime == "uniform":
        return np.full(20, 0.05)
    if regime == "half_glutamine":
        value = background * 0.5
        value[13] += 0.5  # Frozen alphabet ACDEFGHIKLMNPQRSTVWY.
        return value
    raise ValueError("Unknown prespecified composition")


def pairs(seed, regime, length, background, count=PAIRS_PER_SEED):
    if seed not in SEEDS or regime not in REGIMES or length not in LENGTHS or count <= 0:
        raise ValueError("Unknown null panel cell")
    generator = np.random.Generator(np.random.PCG64(
        np.random.SeedSequence([20261002, seed, REGIMES.index(regime), length])))
    weights = probabilities(background, regime)
    query = generator.choice(20, (count, length), p=weights).astype(np.uint8)
    target = generator.choice(20, (count, length), p=weights).astype(np.uint8)
    return query, target


def approximate_e(score, length, lam=0.3176, k=0.134):
    return 1e10 if score <= 0 else k * length * length * math.exp(-lam * score)


def critical_score(length, threshold):
    if type(length) is not int or length <= 0 or not 0 < threshold <= 1:
        raise ValueError("Invalid null search space or threshold")
    score = max(1, math.floor(math.log(0.134 * length * length / threshold) / 0.3176) + 1)
    while approximate_e(score, length) >= threshold:
        score += 1
    while score > 1 and approximate_e(score - 1, length) < threshold:
        score -= 1
    return score


def tail(scores, length, threshold):
    scores = np.asarray(scores)
    if scores.ndim != 1 or not len(scores) or scores.dtype.kind not in "iu" or (scores < 0).any():
        raise ValueError("Require nonnegative integer null scores")
    boundary = critical_score(length, threshold)
    hits, count = int(np.count_nonzero(scores >= boundary)), len(scores)
    result = binomtest(hits, count)
    nominal = result.proportion_ci(confidence_level=0.95, method="exact")
    adjusted = result.proportion_ci(confidence_level=1 - 0.05 / ENDPOINTS, method="exact")
    boundary_e = approximate_e(boundary, length)
    reference = -math.expm1(-boundary_e)
    return dict(threshold=threshold, minimum_integer_score=boundary, trials=count, hits=hits,
                fraction=hits / count, nominal_clopper_pearson=[nominal.low, nominal.high],
                bonferroni_clopper_pearson=[adjusted.low, adjusted.high],
                boundary_approximate_e=boundary_e, poisson_model_tail_reference=reference,
                observed_to_model_tail_ratio=(hits / count) / reference)


def load_frozen(root):
    root = root.resolve()
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip() != REVISION:
        raise ValueError("Wrong frozen scientific checkout")
    if any(name == "orthohmm" or name.startswith("orthohmm.") for name in sys.modules):
        raise ValueError("Require a fresh interpreter without another scientific package")
    files = []
    for name in ("__init__.py", "search/__init__.py", "search/viterbi.py", "search/profile.py",
                 "search/matrices.py", "search/evalue.py", "search/csrc/hmm_viterbi.c"):
        relative = "orthohmm/" + name
        original = subprocess.check_output(["git", "show", REVISION + ":" + relative], cwd=root)
        item = record(root / relative)
        if item["bytes"] != len(original) or item["sha256"] != hashlib.sha256(original).hexdigest():
            raise ValueError("Changed frozen scoring source")
        files.append(item)
    native = record(root / "orthohmm/search/csrc/hmm_viterbi.so")
    if native["sha256"] != NATIVE_SHA:
        raise ValueError("Changed independently retained native kernel")
    files.append(native)
    versions = {name: importlib.metadata.version(name) for name in PACKAGES}
    if versions != PACKAGES:
        raise ValueError("Different null diagnostic runtime")
    sys.path.insert(0, str(root))
    from orthohmm.search import profile, matrices, viterbi, evalue
    for module in (profile, matrices, viterbi, evalue):
        if not Path(module.__file__).resolve().is_relative_to(root):
            raise ValueError("Wrong scientific module origin")
    if not viterbi._load_c_hmm()[1]:
        raise ValueError("Pinned native kernel required; no fallback")
    if matrices.get_ka_params("BLOSUM62") != (0.3176, 0.134):
        raise ValueError("Changed frozen significance constants")
    return (profile, matrices, viterbi, evalue), files, versions


def run(root, protocol, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    protocol_record = record(protocol)
    output.mkdir(parents=True)
    report = dict(status="started", revision=REVISION, protocol=protocol_record,
                  source=record(__file__), python=sys.version, executable=sys.executable,
                  scalar_c_threads=2, jit_threads=2, rows=[], summaries=[],
                  prefilter_executed=False, native_pipeline_rerun=False,
                  empirical_null_diagnostic_executed=True, calibration_established=False,
                  benchmark_scores_or_defaults_changed=False, publication_ready=False)
    try:
        modules, files, versions = load_frozen(root)
        report.update(checked_sources=files, packages=versions)
        profile, matrices, viterbi, evalue = modules
        import numba
        numba.set_num_threads(2)
        matrix = matrices.get_matrix("BLOSUM62")
        background = matrices.get_background_freqs("BLOSUM62")
        report["background"] = background.tolist()
        for regime in REGIMES:
            for length in LENGTHS:
                pooled = {band: [] for band in BANDS}
                for seed in SEEDS:
                    query, target = pairs(seed, regime, length, background)
                    profiles = [profile.build_profile(q, length, matrix, background) for q in query]
                    np.testing.assert_array_equal(profiles[0].transitions, [0, -12, -12, -1, -3, -1, -3])
                    np.testing.assert_array_equal(profiles[0].insert_emissions, np.full(20, -1))
                    count = len(query)
                    offsets = np.arange(count, dtype=np.int64) * length
                    lengths = np.full(count, length, dtype=np.int32)
                    indices = np.column_stack([np.arange(count), np.arange(count)]).astype(np.int32)
                    arguments = (np.concatenate([p.match_emissions for p in profiles]),
                                 profiles[0].insert_emissions, profiles[0].transitions,
                                 offsets, lengths, target.ravel(), offsets, lengths, indices)
                    scored = {}
                    for band in BANDS:
                        scores = viterbi.batch_viterbi_c(*arguments, band, n_threads=2)
                        values = evalue.batch_estimate_evalues(scores, lengths, length, 0.3176, 0.134)
                        for threshold in THRESHOLDS:
                            np.testing.assert_array_equal(values < threshold,
                                                          scores >= critical_score(length, threshold))
                        scored[band] = scores
                        pooled[band].append(scores)
                    if (scored[64] > scored[0]).any():
                        raise ValueError("Restricted band exceeds full-matrix raw score")
                    report["rows"].append(dict(regime=regime, length=length, seed=seed,
                        query_sha256=hashlib.sha256(query.tobytes()).hexdigest(),
                        target_sha256=hashlib.sha256(target.tobytes()).hexdigest(),
                        scores={str(b):scored[b].tolist() for b in BANDS}))
                full, banded = [np.concatenate(pooled[b]) for b in BANDS]
                report["summaries"].append(dict(regime=regime, length=length,
                    probabilities=probabilities(background, regime).tolist(),
                    bands={str(b):[tail(np.concatenate(pooled[b]), length, t) for t in THRESHOLDS]
                           for b in BANDS},
                    band_changed_scores=int(np.count_nonzero(full != banded)),
                    band_lost_gate_hits={str(t):int(np.count_nonzero(
                        (full >= critical_score(length, t)) & (banded < critical_score(length, t))))
                        for t in THRESHOLDS}))
                for item in [protocol_record, report["source"], *files]:
                    check(item)
                (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
                print(regime + " length=" + str(length) + " completed", flush=True)
        report.update(status="prespecified_frozen_null_panel_completed", independent_pairs=90000,
                      native_score_evaluations=180000, tail_endpoints=ENDPOINTS)
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "result.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "protocol", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.protocol, args.output)
