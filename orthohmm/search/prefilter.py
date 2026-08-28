"""DIAMOND-style two-stage k-mer pre-filtering for candidate generation.

Reduces the search space from O(N_query * N_target) to O(N_query * k)
where k << N_target, by:

Phase A: Count shared k-mer hits per target (fast, coarse filter)
Phase B: Check diagonal consistency of hits (selective, fine filter)

Two implementations:
  1. Numba JIT (always available, fallback)
  2. C/OpenMP (fast path with prefetch hints, OpenMP parallel over queries)

Adapted from ClustKIT's kmer_index.py / kmer_score.c.
"""

import ctypes
import logging
import os
from dataclasses import dataclass
import numpy as np
from numba import njit, prange, int32, int64, int16
from pathlib import Path as _Path

_logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────
# C/OpenMP prefilter extension (lazy-loaded)
# ──────────────────────────────────────────────────────────────────────

_c_pf_lib = None
_c_pf_available = None


def _load_c_prefilter():
    """Load the C/OpenMP prefilter shared library."""
    global _c_pf_lib, _c_pf_available
    if _c_pf_available is not None:
        return _c_pf_lib, _c_pf_available
    try:
        so_path = _Path(__file__).resolve().parent / "csrc" / "kmer_prefilter.so"
        lib = ctypes.cdll.LoadLibrary(str(so_path))
        lib.batch_prefilter_c.restype = None
        lib.batch_prefilter_c.argtypes = [
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,  # q_flat, q_offsets, q_lengths
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,      # nq, k, alpha_size
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,  # kmer_offsets, entries, freqs
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,      # freq_thresh, num_db, min_total
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,      # min_diag, diag_bw, max_cands
            ctypes.c_int32,                                       # phase_a_topk
            ctypes.c_void_p, ctypes.c_void_p,                    # out_targets, out_counts
            ctypes.c_void_p, ctypes.c_void_p,                    # q_flat_orig, t_flat_orig
            ctypes.c_void_p, ctypes.c_void_p,                    # t_offsets_orig, t_lengths_orig
            ctypes.c_void_p,                                      # sub_matrix
        ]
        lib.prefilter_set_num_threads.restype = None
        lib.prefilter_set_num_threads.argtypes = [ctypes.c_int32]
        _c_pf_lib = lib
        _c_pf_available = True
    except (OSError, AttributeError):
        _c_pf_lib = None
        _c_pf_available = False
    return _c_pf_lib, _c_pf_available


def is_c_prefilter_available():
    """Check if C/OpenMP prefilter is available."""
    _, available = _load_c_prefilter()
    return available


# ──────────────────────────────────────────────────────────────────────
# Reduced alphabet for sensitive distant homology detection
# Murphy-10-like (9 groups), same as ClustKIT
# ──────────────────────────────────────────────────────────────────────

REDUCED_ALPHA = np.array(
    [0, 1, 2, 2, 3, 0, 4, 5, 6, 5, 5, 2, 7, 2, 6, 8, 8, 5, 3, 3],
    dtype=np.uint8,
)
REDUCED_ALPHA_SIZE = 9


@njit(cache=True)
def _remap_flat(src, mapping, n):
    """Remap a flat uint8 array through a lookup table."""
    out = np.empty(n, dtype=np.uint8)
    for i in range(n):
        v = src[i]
        if v < 20:
            out[i] = mapping[v]
        else:
            out[i] = 20
    return out


# ──────────────────────────────────────────────────────────────────────
# K-mer index building (CSR format)
# ──────────────────────────────────────────────────────────────────────

@njit(cache=True)
def _count_kmers(flat_seqs, offsets, lengths, k, alpha_size, counts):
    """Pass 1: count k-mer occurrences across all sequences."""
    n = len(lengths)
    for i in range(n):
        start = int64(offsets[i])
        length = int32(lengths[i])
        if length < k:
            continue
        num_kmers = length - k + 1
        for pos in range(num_kmers):
            kmer_val = int64(0)
            valid = True
            for j in range(k):
                r = flat_seqs[start + pos + j]
                if r >= alpha_size:
                    valid = False
                    break
                kmer_val = kmer_val * int64(alpha_size) + int64(r)
            if valid:
                counts[kmer_val] += int64(1)


@njit(cache=True)
def _fill_entries(flat_seqs, offsets, lengths, k, alpha_size,
                  kmer_offsets, entries, cursors):
    """Pass 2: fill CSR entries with packed (seq_id << 32 | position)."""
    n = len(lengths)
    for i in range(n):
        start = int64(offsets[i])
        length = int32(lengths[i])
        if length < k:
            continue
        num_kmers = length - k + 1
        for pos in range(num_kmers):
            kmer_val = int64(0)
            valid = True
            for j in range(k):
                r = flat_seqs[start + pos + j]
                if r >= alpha_size:
                    valid = False
                    break
                kmer_val = kmer_val * int64(alpha_size) + int64(r)
            if valid:
                idx = cursors[kmer_val]
                entries[idx] = (int64(i) << 32) | int64(pos)
                cursors[kmer_val] = idx + int64(1)


def build_kmer_index(flat_sequences, offsets, lengths, k=5, alpha_size=9):
    """Build a CSR-format k-mer inverted index.

    Returns (kmer_offsets, kmer_entries, kmer_freqs).
    """
    num_possible = alpha_size ** k

    counts = np.zeros(num_possible, dtype=np.int64)
    _count_kmers(flat_sequences, offsets, lengths,
                 int32(k), int32(alpha_size), counts)

    kmer_freqs = counts.astype(np.int32)

    kmer_offsets = np.zeros(num_possible + 1, dtype=np.int64)
    np.cumsum(counts, out=kmer_offsets[1:])
    total = int(kmer_offsets[-1])

    kmer_entries = np.empty(total, dtype=np.int64)
    cursors = kmer_offsets[:-1].copy()
    _fill_entries(flat_sequences, offsets, lengths,
                  int32(k), int32(alpha_size),
                  kmer_offsets, kmer_entries, cursors)

    return kmer_offsets, kmer_entries, kmer_freqs


def compute_freq_threshold(kmer_freqs, num_sequences, percentile=99.5):
    """Compute frequency cap for overly common k-mers."""
    nonzero = kmer_freqs[kmer_freqs > 0]
    if len(nonzero) == 0:
        return int32(num_sequences)
    pctl = int(np.percentile(nonzero, percentile))
    floor = max(100, num_sequences // 200)
    return int32(max(pctl, floor))


@dataclass(frozen=True)
class PreparedKmerIndex:
    """Reusable target-side arrays for one prefilter configuration."""

    k: int
    use_reduced_alphabet: bool
    alpha_size: int
    target_flat: np.ndarray
    kmer_offsets: np.ndarray
    kmer_entries: np.ndarray
    kmer_freqs: np.ndarray
    freq_threshold: int


def prepare_kmer_index(target_species, k=5, use_reduced_alphabet=True):
    """Build the target index once for adaptive probing and filtering."""
    if use_reduced_alphabet:
        alpha_size = REDUCED_ALPHA_SIZE
        target_flat = _remap_flat(
            target_species.flat_sequences,
            REDUCED_ALPHA,
            len(target_species.flat_sequences),
        )
    else:
        alpha_size = 20
        target_flat = target_species.flat_sequences

    kmer_offsets, kmer_entries, kmer_freqs = build_kmer_index(
        target_flat,
        target_species.offsets,
        target_species.lengths,
        k=k,
        alpha_size=alpha_size,
    )
    return PreparedKmerIndex(
        k=int(k),
        use_reduced_alphabet=bool(use_reduced_alphabet),
        alpha_size=alpha_size,
        target_flat=target_flat,
        kmer_offsets=kmer_offsets,
        kmer_entries=kmer_entries,
        kmer_freqs=kmer_freqs,
        freq_threshold=int(compute_freq_threshold(
            kmer_freqs, target_species.num_sequences
        )),
    )


# ──────────────────────────────────────────────────────────────────────
# Two-stage query scoring
# ──────────────────────────────────────────────────────────────────────

@njit(cache=True)
def _score_query_two_stage(
    q_seq, q_len, k, alpha_size,
    kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
    num_db, min_total_hits, min_diag_hits, diag_bin_width,
):
    """Two-stage k-mer scoring for one query.

    Phase A: count hits per target, keep top candidates.
    Phase B: diagonal consistency check on Phase A survivors.

    Returns (target_ids, scores) int32 arrays.
    """
    num_kmers = q_len - k + 1
    if num_kmers <= 0:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)

    # Phase A: k-mer counting
    target_counts = np.zeros(num_db, dtype=np.int16)

    for qpos in range(num_kmers):
        kmer_val = int64(0)
        valid = True
        for j in range(k):
            r = q_seq[qpos + j]
            if r >= alpha_size:
                valid = False
                break
            kmer_val = kmer_val * int64(alpha_size) + int64(r)
        if not valid:
            continue
        if kmer_freqs[kmer_val] > freq_thresh:
            continue
        s = kmer_offsets[kmer_val]
        e = kmer_offsets[kmer_val + 1]
        for h in range(s, e):
            tid = int32(kmer_entries[h] >> 32)
            if target_counts[tid] < 32767:
                target_counts[tid] += int16(1)

    # Build survivor mask
    survivor_mask = np.zeros(num_db, dtype=np.bool_)
    n_surv_hits = int64(0)
    for i in range(num_db):
        if target_counts[i] >= min_total_hits:
            survivor_mask[i] = True
            n_surv_hits += int64(target_counts[i])

    if n_surv_hits == 0:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)

    if min_diag_hits <= 1:
        # Skip Phase B: return all Phase A survivors
        num_survivors = int32(0)
        for i in range(num_db):
            if survivor_mask[i]:
                num_survivors += int32(1)
        result_ids = np.empty(num_survivors, dtype=np.int32)
        result_scores = np.empty(num_survivors, dtype=np.int32)
        pos = int32(0)
        for i in range(num_db):
            if survivor_mask[i]:
                result_ids[pos] = int32(i)
                result_scores[pos] = int32(target_counts[i])
                pos += int32(1)
        return result_ids, result_scores

    # Phase B: diagonal consistency
    DIAG_MULT = int64(1000000)
    max_diag_shift = int32(q_len)

    # Allocate with generous padding — the second pass recomputes k-mers
    # and may encounter different rounding, so use 2x safety margin
    surv_keys = np.empty(n_surv_hits + int64(num_kmers) + int64(1000), dtype=np.int64)
    sw = int64(0)
    surv_cap = int64(len(surv_keys))

    for qpos in range(num_kmers):
        kmer_val = int64(0)
        valid = True
        for j in range(k):
            r = q_seq[qpos + j]
            if r >= alpha_size:
                valid = False
                break
            kmer_val = kmer_val * int64(alpha_size) + int64(r)
        if not valid:
            continue
        if kmer_freqs[kmer_val] > freq_thresh:
            continue
        s = kmer_offsets[kmer_val]
        e = kmer_offsets[kmer_val + 1]
        for h in range(s, e):
            entry = kmer_entries[h]
            tid = int32(entry >> 32)
            if survivor_mask[tid]:
                if sw < surv_cap:
                    tpos = int32(entry & 0xFFFFFFFF)
                    diag = int32(tpos) - int32(qpos) + max_diag_shift
                    dbin = int32(diag // diag_bin_width)
                    surv_keys[sw] = int64(tid) * DIAG_MULT + int64(dbin)
                    sw += int64(1)

    surv_keys = surv_keys[:sw]
    surv_keys.sort()

    # Count runs per (tid, dbin), keep best diagonal per target
    max_results = int32(0)
    for i in range(num_db):
        if survivor_mask[i]:
            max_results += int32(1)

    final_ids = np.empty(max_results, dtype=np.int32)
    final_scores = np.empty(max_results, dtype=np.int32)
    num_final = int32(0)
    prev_tid = int32(-1)
    best_count = int32(0)
    i = int64(0)
    n = int64(len(surv_keys))

    while i < n:
        key = surv_keys[i]
        tid = int32(key // DIAG_MULT)
        run = int32(0)
        cur_dbin = int32(key % DIAG_MULT)
        while i < n:
            k2 = surv_keys[i]
            if int32(k2 // DIAG_MULT) != tid or int32(k2 % DIAG_MULT) != cur_dbin:
                break
            run += int32(1)
            i += int64(1)
        if tid != prev_tid:
            if prev_tid >= 0 and best_count >= min_diag_hits:
                final_ids[num_final] = prev_tid
                final_scores[num_final] = best_count
                num_final += int32(1)
            prev_tid = tid
            best_count = run
        else:
            if run > best_count:
                best_count = run

    if prev_tid >= 0 and best_count >= min_diag_hits:
        final_ids[num_final] = prev_tid
        final_scores[num_final] = best_count
        num_final += int32(1)

    return final_ids[:num_final], final_scores[:num_final]


@njit(cache=True)
def batch_prefilter(
    query_flat, query_offsets, query_lengths,
    k, alpha_size,
    kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
    num_db, min_total_hits, min_diag_hits, diag_bin_width,
    max_candidates_per_query,
):
    """Pre-filtering for all queries.

    Returns:
      all_candidate_ids: int32 array (flat, concatenated)
      all_candidate_offsets: int64 array (N_queries + 1), CSR pointers
    """
    N = len(query_lengths)

    # First pass: count candidates per query
    counts = np.empty(N, dtype=np.int32)
    total_candidates = int64(0)

    for qi in range(N):
        q_off = query_offsets[qi]
        q_len = query_lengths[qi]
        q_seq = query_flat[q_off:q_off + q_len]

        cand_ids, cand_scores = _score_query_two_stage(
            q_seq, int32(q_len), int32(k), int32(alpha_size),
            kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
            int32(num_db), int32(min_total_hits),
            int32(min_diag_hits), int32(diag_bin_width),
        )

        n_cands = min(len(cand_ids), max_candidates_per_query)
        counts[qi] = int32(n_cands)
        total_candidates += int64(n_cands)

    # Build CSR offsets
    all_candidate_offsets = np.zeros(N + 1, dtype=np.int64)
    for i in range(N):
        all_candidate_offsets[i + 1] = all_candidate_offsets[i] + int64(counts[i])

    # Second pass: fill candidates
    all_candidate_ids = np.empty(total_candidates, dtype=np.int32)

    for qi in range(N):
        q_off = query_offsets[qi]
        q_len = query_lengths[qi]
        q_seq = query_flat[q_off:q_off + q_len]

        cand_ids, cand_scores = _score_query_two_stage(
            q_seq, int32(q_len), int32(k), int32(alpha_size),
            kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
            int32(num_db), int32(min_total_hits),
            int32(min_diag_hits), int32(diag_bin_width),
        )

        n_cands = min(len(cand_ids), max_candidates_per_query)
        out_off = all_candidate_offsets[qi]
        for j in range(n_cands):
            all_candidate_ids[out_off + j] = cand_ids[j]

    return all_candidate_ids, all_candidate_offsets


def _prefilter_c(
    query_flat, query_offsets, query_lengths,
    k, alpha_size,
    kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
    num_db, min_total_hits, min_diag_hits, diag_bin_width,
    max_candidates_per_query, phase_a_topk=5000, n_threads=4,
    query_flat_original=None, target_flat_original=None,
    target_offsets_original=None, target_lengths_original=None,
    sub_matrix=None,
):
    """C/OpenMP batch prefilter with ungapped extension re-ranking.

    If original sequences and sub_matrix are provided, candidates are
    re-ranked by ungapped extension score (Kadane's on 3 diagonals)
    instead of raw k-mer count. This dramatically improves candidate
    quality for downstream RBH detection.
    """
    lib, available = _load_c_prefilter()
    if not available:
        raise RuntimeError("C prefilter not available")

    lib.prefilter_set_num_threads(n_threads)

    nq = len(query_lengths)
    max_cands = max_candidates_per_query

    out_targets = np.empty(nq * max_cands, dtype=np.int32)
    out_counts = np.empty(nq, dtype=np.int32)

    _qf = np.ascontiguousarray(query_flat, dtype=np.uint8)
    _qo = np.ascontiguousarray(query_offsets, dtype=np.int64)
    _ql = np.ascontiguousarray(query_lengths, dtype=np.int32)
    _ko = np.ascontiguousarray(kmer_offsets, dtype=np.int64)
    _ke = np.ascontiguousarray(kmer_entries, dtype=np.int64)
    _kf = np.ascontiguousarray(kmer_freqs, dtype=np.int32)

    # Extension scoring arrays (NULL if not provided)
    use_ext = (query_flat_original is not None and sub_matrix is not None)
    if use_ext:
        _qfo = np.ascontiguousarray(query_flat_original, dtype=np.uint8)
        _tfo = np.ascontiguousarray(target_flat_original, dtype=np.uint8)
        _too = np.ascontiguousarray(target_offsets_original, dtype=np.int64)
        _tlo = np.ascontiguousarray(target_lengths_original, dtype=np.int32)
        _sm = np.ascontiguousarray(sub_matrix, dtype=np.int8)
        qfo_ptr = _qfo.ctypes.data
        tfo_ptr = _tfo.ctypes.data
        too_ptr = _too.ctypes.data
        tlo_ptr = _tlo.ctypes.data
        sm_ptr = _sm.ctypes.data
    else:
        qfo_ptr = None
        tfo_ptr = None
        too_ptr = None
        tlo_ptr = None
        sm_ptr = None

    lib.batch_prefilter_c(
        _qf.ctypes.data, _qo.ctypes.data, _ql.ctypes.data,
        nq, k, alpha_size,
        _ko.ctypes.data, _ke.ctypes.data, _kf.ctypes.data,
        freq_thresh, num_db, min_total_hits,
        min_diag_hits, diag_bin_width, max_cands, phase_a_topk,
        out_targets.ctypes.data, out_counts.ctypes.data,
        qfo_ptr, tfo_ptr, too_ptr, tlo_ptr, sm_ptr,
    )

    # Convert flat output to CSR format
    all_candidate_offsets = np.zeros(nq + 1, dtype=np.int64)
    np.cumsum(out_counts.astype(np.int64), out=all_candidate_offsets[1:])
    total = int(all_candidate_offsets[-1])

    all_candidate_ids = np.empty(total, dtype=np.int32)
    for qi in range(nq):
        nc = int(out_counts[qi])
        src_off = qi * max_cands
        dst_off = int(all_candidate_offsets[qi])
        all_candidate_ids[dst_off:dst_off + nc] = out_targets[src_off:src_off + nc]

    return all_candidate_ids, all_candidate_offsets


def prefilter_candidates(
    query_species,
    target_species,
    k=5,
    use_reduced_alphabet=True,
    min_total_hits=3,
    min_diag_hits=2,
    diag_bin_width=6,
    max_candidates_per_query=500,
    n_threads=4,
    sub_matrix=None,
    prepared_index=None,
):
    """K-mer pre-filtering with ungapped extension re-ranking.

    Pipeline per query:
      Phase A: k-mer counting → candidate targets
      Phase A.5: Ungapped extension (Kadane's on 3 diags) → re-rank
      Top-K by extension score → high-quality candidates

    Parameters
    ----------
    query_species : SpeciesSequences
    target_species : SpeciesSequences
    k : k-mer size
    use_reduced_alphabet : use Murphy-10 reduced alphabet for k-mer lookup
    min_total_hits : Phase A threshold
    min_diag_hits : Phase B diagonal threshold (1=skip Phase B)
    diag_bin_width : diagonal binning granularity
    max_candidates_per_query : cap on candidates per query
    n_threads : OpenMP thread count (C path only)
    sub_matrix : 20x20 int8 substitution matrix for extension scoring

    Returns
    -------
    all_candidate_ids : int32 array (flat)
    all_candidate_offsets : int64 array (N_queries + 1)
    """
    if prepared_index is None:
        prepared_index = prepare_kmer_index(
            target_species,
            k=k,
            use_reduced_alphabet=use_reduced_alphabet,
        )
    elif (
        prepared_index.k != int(k)
        or prepared_index.use_reduced_alphabet != bool(use_reduced_alphabet)
    ):
        raise ValueError("prepared k-mer index does not match prefilter settings")

    alpha_size = prepared_index.alpha_size
    target_flat = prepared_index.target_flat
    kmer_offsets = prepared_index.kmer_offsets
    kmer_entries = prepared_index.kmer_entries
    kmer_freqs = prepared_index.kmer_freqs
    freq_thresh = int32(prepared_index.freq_threshold)

    if use_reduced_alphabet:
        query_flat = _remap_flat(
            query_species.flat_sequences, REDUCED_ALPHA,
            len(query_species.flat_sequences)
        )
    else:
        query_flat = query_species.flat_sequences

    num_db = target_species.num_sequences

    # Prefer C/OpenMP path
    _, c_available = _load_c_prefilter()
    if c_available:
        return _prefilter_c(
            query_flat, query_species.offsets, query_species.lengths,
            k, alpha_size,
            kmer_offsets, kmer_entries, kmer_freqs, int32(freq_thresh),
            num_db, min_total_hits, min_diag_hits, diag_bin_width,
            max_candidates_per_query, n_threads=n_threads,
            # Pass original sequences for extension re-ranking
            query_flat_original=query_species.flat_sequences,
            target_flat_original=target_species.flat_sequences,
            target_offsets_original=target_species.offsets,
            target_lengths_original=target_species.lengths,
            sub_matrix=sub_matrix,
        )

    # Numba fallback (no extension scoring)
    return batch_prefilter(
        query_flat, query_species.offsets, query_species.lengths,
        int32(k), int32(alpha_size),
        kmer_offsets, kmer_entries, kmer_freqs, freq_thresh,
        int32(num_db), int32(min_total_hits),
        int32(min_diag_hits), int32(diag_bin_width),
        int32(max_candidates_per_query),
    )
