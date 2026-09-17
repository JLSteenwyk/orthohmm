"""Causal diagnostic of per-pair banding; never installs a benchmark binary."""

import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_native_scoring import COMMIT, BANDS, fixtures, flat, digest


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def summarize(report):
    if [row["band"] for row in report["bands"]] != list(BANDS):
        raise ValueError("Incomplete diagnostic band inventory")
    rows = []
    for row in report["bands"]:
        arrays = [np.asarray(row[key]) for key in ("baseline", "scalar", "jit", "variant")]
        if any(a.shape != (report["pairs_per_band"],) or a.dtype.kind not in "iu" for a in arrays):
            raise ValueError("Incomplete integer diagnostic scores")
        baseline, scalar, jit, variant = arrays
        counts = {"baseline_scalar_mismatches": int(np.count_nonzero(baseline != scalar)),
                  "variant_scalar_mismatches": int(np.count_nonzero(variant != scalar)),
                  "variant_jit_mismatches": int(np.count_nonzero(variant != jit))}
        changes = np.flatnonzero(baseline != variant).tolist()
        if any(row[key] != value for key, value in counts.items()) or changes != row["changed_indices"]:
            raise ValueError("Diagnostic counters differ from raw scores")
        rows.append({"band": row["band"], **counts, "changed_indices": changes})
    return {"bands": rows, "pairs_per_band": report["pairs_per_band"],
            "variant_matches_scalar_and_jit": all(not r["variant_scalar_mismatches"] and
                                                  not r["variant_jit_mismatches"] for r in rows),
            "baseline_discrepancy_observed": any(r["baseline_scalar_mismatches"] for r in rows),
            "frozen_runtime_modified": report["frozen_runtime_modified"], "benchmark_admitted": False}


def diagnostic(root, library):
    root, library = root.resolve(), library.resolve()
    if library.is_relative_to(root):
        raise ValueError("Diagnostic library must be outside frozen checkout")
    if subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Wrong frozen revision")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--", "orthohmm"], check=True)
    sys.path.insert(0, str(root))
    from orthohmm.search import viterbi, matrices, profile
    for module in (viterbi, matrices, profile):
        if not Path(module.__file__).resolve().is_relative_to(root):
            raise ValueError("Wrong imported checkout")
    original, available = viterbi._load_c_hmm()
    if not available or not original.hmm_have_avx2():
        raise ValueError("Require frozen AVX2 backend")
    variant = ctypes.CDLL(str(library))
    for name in ("batch_hmm_viterbi_c", "batch_hmm_viterbi_multipair_avx2_c", "hmm_set_num_threads", "hmm_have_avx2"):
        target, source = getattr(variant, name), getattr(original, name)
        target.argtypes, target.restype = source.argtypes, source.restype
    if not variant.hmm_have_avx2():
        raise ValueError("Diagnostic build lacks AVX2")
    files = [Path(__file__), library, library.with_suffix(".c"),
             root / "orthohmm/search/csrc/hmm_viterbi.so",
             root / "orthohmm/search/csrc/hmm_viterbi.c"]
    before = [record(p) for p in files]
    queries, targets, pairs = fixtures()
    matrix, background = matrices.get_matrix("BLOSUM62"), matrices.get_background_freqs("BLOSUM62")
    profiles = [profile.build_profile(q, len(q), matrix, background) for q in queries]
    _, offsets, lengths = flat(queries)
    target_flat, target_offsets, target_lengths = flat(targets)
    args = (np.concatenate([p.match_emissions for p in profiles]), profiles[0].insert_emissions,
            profiles[0].transitions, offsets, lengths, target_flat, target_offsets, target_lengths, pairs)
    rows = []
    try:
        for band in BANDS:
            viterbi._c_hmm_lib = original
            scalar = viterbi.batch_viterbi_c(*args, band, n_threads=1)
            jit = viterbi.batch_viterbi_score(*args, band)
            baseline = viterbi.batch_viterbi_multipair_c(*args, band, n_threads=4)
            viterbi._c_hmm_lib = variant
            rescued = viterbi.batch_viterbi_multipair_c(*args, band, n_threads=4)
            rows.append({"band": band, "baseline_scalar_mismatches": int(np.count_nonzero(baseline != scalar)),
                         "variant_scalar_mismatches": int(np.count_nonzero(rescued != scalar)),
                         "variant_jit_mismatches": int(np.count_nonzero(rescued != jit)),
                         "changed_indices": np.flatnonzero(baseline != rescued).tolist(),
                         "baseline": baseline.tolist(), "scalar": scalar.tolist(),
                         "jit": jit.tolist(), "variant": rescued.tolist()})
    finally:
        viterbi._c_hmm_lib = original
    if before != [record(p) for p in files]:
        raise ValueError("Probe sources or binaries changed")
    fixture = {"queries": [q.tolist() for q in queries], "targets": [t.tolist() for t in targets],
               "pairs": pairs.tolist(), "bands": list(BANDS)}
    return {"status": "diagnostic_complete", "core_commit": COMMIT, "files": before,
            "fixture_sha256": digest(fixture), "pairs_per_band": len(pairs), "bands": rows,
            "frozen_runtime_modified": False, "benchmark_admitted": False,
            "limitations": ["Single isolated patch and finite fixtures; not general backend equivalence.",
                            "No publication baseline, historical scores or runtime manifests were changed.",
                            "Release integration requires broader tests and an explicitly versioned method change."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = diagnostic(args.root, args.library)
    report["summary"] = summarize(report)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
