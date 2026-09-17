"""Deterministic native HMM portability probe, separate from benchmark accuracy."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

SEED = 20260917
BANDS = (0, 1, 8, 64, 128)
LENGTHS = (1, 7, 8, 9, 31, 64, 129, 257, 512)
COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fixtures():
    rng = np.random.Generator(np.random.PCG64(SEED))
    queries, targets = [], []
    for length in LENGTHS:
        query = rng.integers(0, 20, length, dtype=np.uint8)
        queries.append(query)
        ambiguous = query.copy()
        ambiguous[::7] = 20
        targets.extend([query.copy(), ambiguous, query[:max(1, length // 2)],
                        np.concatenate([query[:length // 2], np.zeros(11, dtype=np.uint8), query[length // 2:]]),
                        rng.integers(0, 20, length + 3, dtype=np.uint8)])
    # Include ambiguous query residues as well as ambiguous target residues.
    queries[-2] = queries[-2].copy()
    queries[-2][::13] = 20
    pairs = np.array([(q, t) for q in range(len(queries)) for t in range(len(targets))], dtype=np.int32)
    rng.shuffle(pairs)
    return queries, targets, pairs


def flat(sequences):
    lengths = np.array([len(s) for s in sequences], dtype=np.int32)
    offsets = np.concatenate([np.zeros(1, dtype=np.int64), np.cumsum(lengths[:-1], dtype=np.int64)])
    return np.concatenate(sequences), offsets, lengths


def probe(root):
    root = Path(root).resolve()
    if subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Wrong frozen revision")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--", "orthohmm"], check=True)
    sys.path.insert(0, str(root))
    from orthohmm.search import viterbi, profile, matrices, evalue
    for module in (viterbi, profile, matrices, evalue):
        if not Path(module.__file__).resolve().is_relative_to(root):
            raise ValueError("Imported wrong checkout")
    lib, available = viterbi._load_c_hmm()
    if not available:
        raise ValueError("Native HMM required; no silent JIT substitution")
    use_multipair = bool(lib.hmm_have_avx2()) and hasattr(lib, "batch_hmm_viterbi_multipair_avx2_c")
    files = [Path(m.__file__) for m in (viterbi, profile, matrices, evalue)]
    files.append(root / "orthohmm/search/csrc/hmm_viterbi.so")
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    queries, targets, pairs = fixtures()
    matrix, background = matrices.get_matrix("BLOSUM62"), matrices.get_background_freqs("BLOSUM62")
    profiles = [profile.build_profile(q, len(q), matrix, background) for q in queries]
    _, offsets, lengths = flat(queries)
    target_flat, target_offsets, target_lengths = flat(targets)
    args = (np.concatenate([p.match_emissions for p in profiles]), profiles[0].insert_emissions,
            profiles[0].transitions, offsets, lengths, target_flat, target_offsets, target_lengths, pairs)
    fixture = {"queries": [q.tolist() for q in queries], "targets": [t.tolist() for t in targets],
               "pairs": pairs.tolist(), "bands": list(BANDS)}
    report = {"status": "checked", "root": str(root), "commit": COMMIT, "machine": platform.machine(),
              "python": sys.version, "probe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "packages": {p: importlib.metadata.version(p) for p in ("numpy", "numba", "llvmlite")},
              "files": before, "fixture_sha256": digest(fixture), "pairs": len(pairs), "bands": [],
              "backend": "multipair_AVX2" if use_multipair else "scalar_C",
              "limitations": ["Finite synthetic HMM fixtures only, not end-to-end pipeline equivalence or accuracy."]}
    lam, k = matrices.get_ka_params("BLOSUM62")
    for band in BANDS:
        reference = viterbi.batch_viterbi_score(*args, band)
        scalar = viterbi.batch_viterbi_c(*args, band, n_threads=1)
        native = viterbi.batch_viterbi_multipair_c if use_multipair else viterbi.batch_viterbi_c
        scores = native(*args, band, n_threads=4)
        values = evalue.batch_estimate_evalues(scores, lengths[pairs[:, 0]], int(target_lengths.sum()), lam, k)
        normalized = scores / np.sqrt(lengths[pairs[:, 0]].astype(float) * target_lengths[pairs[:, 1]])
        row = {"band": band, "scores": scores.tolist(), "scalar_scores": scalar.tolist(),
               "jit_scores": reference.tolist(), "evalues": values.tolist(), "normalized": normalized.tolist(),
               "decisions": {str(t): (values <= t).tolist() for t in (1e-3, 1e-4, 1e-5)},
               "scalar_mismatches": int(np.count_nonzero(scores != scalar)),
               "jit_mismatches": int(np.count_nonzero(scores != reference))}
        if row["scalar_mismatches"] or row["jit_mismatches"]:
            report["status"] = "mismatch"
        report["bands"].append(row)
    if before != {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}:
        raise ValueError("Runtime changed during probe")
    return report


def compare(left, right):
    for key in ("fixture_sha256", "probe_sha256", "commit", "pairs", "packages"):
        if left[key] != right[key]:
            raise ValueError("Probe provenance differs: " + key)
    if left["status"] != "checked" or right["status"] != "checked":
        raise ValueError("Local backend discrepancy")
    if [r["band"] for r in left["bands"]] != list(BANDS) or [r["band"] for r in right["bands"]] != list(BANDS):
        raise ValueError("Incomplete band checks")
    for a, b in zip(left["bands"], right["bands"]):
        for row in (a, b):
            if row["scalar_mismatches"] or row["jit_mismatches"]:
                raise ValueError("Backend discrepancy counters")
            expected = {str(t): (np.asarray(row["evalues"]) <= t).tolist() for t in (1e-3, 1e-4, 1e-5)}
            if row["decisions"] != expected:
                raise ValueError("Invalid threshold decisions")
        for key in ("scores", "scalar_scores", "jit_scores"):
            if len(a[key]) != left["pairs"] or a[key] != b[key] or a[key] != a["scores"]:
                raise ValueError("Integer scores differ: " + key)
        if a["decisions"] != b["decisions"]:
            raise ValueError("Threshold decisions differ")
        for key in ("evalues", "normalized"):
            if len(a[key]) != left["pairs"] or len(b[key]) != left["pairs"]:
                raise ValueError("Incomplete floating scores")
            if not np.isfinite(a[key]).all() or not np.isfinite(b[key]).all():
                raise ValueError("Nonfinite scores")
            if not np.allclose(a[key], b[key], rtol=1e-12, atol=0):
                raise ValueError("Floating scores differ: " + key)
    return {"status": "synthetic_hmm_scores_match", "pairs_per_band": left["pairs"], "bands": list(BANDS),
            "backends": [left["backend"], right["backend"]], "relative_tolerance": 1e-12,
            "end_to_end_equivalence": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--root", type=Path)
    mode.add_argument("--compare", nargs=2, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = probe(args.root) if args.root else compare(*(json.loads(p.read_text()) for p in args.compare))
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    raise SystemExit(0 if result["status"] != "mismatch" else 1)
