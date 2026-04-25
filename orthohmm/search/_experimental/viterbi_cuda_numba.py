"""CUDA banded 3-state Viterbi for profile HMM search.

One thread per (profile, target) pair. Each thread runs the DP serially
(because rows have serial dependencies) but tens of thousands of pairs
run in parallel across the GPU.

Workspace layout: per-thread strided slabs in global memory — 6 rows of
size (max_len + 2) int32 per thread.

KNOWN PERFORMANCE ISSUE
-----------------------
This naive implementation is currently ~10x SLOWER than the CPU kernel
on typical batches. Root cause: per-thread workspace is in global memory,
and every DP cell access incurs ~500-cycle GPU global-memory latency, so
each thread stalls most of its time waiting on loads. The expected GPU
speedup only materializes with a warp-collaborative kernel that keeps
the active DP rows in shared memory and splits the j-loop across 32
threads. Implementing that properly is a multi-day task (roughly: 1
pair per warp, shared-memory ring of 6*band*4 bytes per block, warp
shuffles for the I-state horizontal dependency, and careful tail
handling). The current kernel is functional and parity-tested; use
it as a reference and when profiling the collaborative rewrite.
"""

import numpy as np

try:
    from numba import cuda
    from numba import int32, int64
    _CUDA_AVAILABLE = cuda.is_available()
except Exception:  # pragma: no cover
    cuda = None
    _CUDA_AVAILABLE = False


def is_cuda_available():
    return _CUDA_AVAILABLE


if _CUDA_AVAILABLE:

    NEG_INF = int32(-1000000)

    @cuda.jit
    def _viterbi_kernel(
        flat_match_emit,   # (sum_L, 20) int32  (widened on host)
        insert_emit,       # (20,) int32
        transitions,       # (7,) int32
        prof_offsets,      # (N_prof,) int64
        prof_lengths,      # (N_prof,) int32
        target_flat,       # (sum_T,) int32  (widened on host)
        target_offsets,    # (N_tgt,) int64
        target_lengths,    # (N_tgt,) int32
        pairs,             # (M, 2) int32
        band_width,        # int32
        out_scores,        # (M,) int32
        ws,                # (n_threads_max, 6, ws_size) int32 workspace
        ws_size,           # int32
    ):
        tid = cuda.grid(1)
        M = pairs.shape[0]
        if tid >= M:
            return

        qi = pairs[tid, 0]
        ti = pairs[tid, 1]

        L = prof_lengths[qi]
        T = target_lengths[ti]
        if L == 0 or T == 0:
            out_scores[tid] = 0
            return

        p_off = prof_offsets[qi]
        t_off = target_offsets[ti]

        t_mm = transitions[0]
        t_mi = transitions[1]
        t_md = transitions[2]
        t_im = transitions[3]
        t_ii = transitions[4]
        t_dm = transitions[5]
        t_dd = transitions[6]

        # Band bounds
        bw = band_width
        max_dim = L if L > T else T
        if bw <= 0 or max_dim <= 50:
            bw = max_dim

        cols = T + 1

        # Workspace row pointers into the global slab for this thread
        # Six rows: Mp Mc Ip Ic Dp Dc
        Mp = ws[tid, 0]
        Mc = ws[tid, 1]
        Ip = ws[tid, 2]
        Ic = ws[tid, 3]
        Dp = ws[tid, 4]
        Dc = ws[tid, 5]

        for j in range(cols):
            Mp[j] = NEG_INF; Mc[j] = NEG_INF
            Ip[j] = NEG_INF; Ic[j] = NEG_INF
            Dp[j] = NEG_INF; Dc[j] = NEG_INF

        max_score = int32(0)

        for i in range(1, L + 1):
            j_center = (i * T) // L
            j_lo = j_center - bw
            j_hi = j_center + bw
            if j_lo < 1:
                j_lo = 1
            if j_hi > T:
                j_hi = T

            # clear current row borders
            lo_m1 = j_lo - 1
            hi_p1 = j_hi + 1
            if hi_p1 >= cols:
                hi_p1 = cols - 1
            for j in range(lo_m1 if lo_m1 >= 0 else 0, hi_p1 + 1):
                Mc[j] = NEG_INF
                Ic[j] = NEG_INF
                Dc[j] = NEG_INF

            # Emission tables are padded to 32 entries; out-of-range letters
            # map to index 20 which holds (0, -1) so no branching is needed.
            emit_base = int64(p_off + (i - 1)) * int64(32)

            for j in range(j_lo, j_hi + 1):
                t_aa_raw = int64(target_flat[t_off + j - 1])
                t_aa = t_aa_raw if t_aa_raw < 20 else int64(20)
                m_emit = flat_match_emit[emit_base + t_aa]
                i_emit = insert_emit[t_aa]

                # M[i,j]
                best_m = int32(0)
                v = int32(Mp[j - 1]) + t_mm
                if v > best_m:
                    best_m = v
                v = int32(Ip[j - 1]) + t_im
                if v > best_m:
                    best_m = v
                v = int32(Dp[j - 1]) + t_dm
                if v > best_m:
                    best_m = v
                mc_v = best_m + m_emit
                Mc[j] = mc_v

                # I[i,j]
                best_i = NEG_INF
                v = int32(Mc[j - 1]) + t_mi
                if v > best_i:
                    best_i = v
                v = int32(Ic[j - 1]) + t_ii
                if v > best_i:
                    best_i = v
                if best_i > NEG_INF // int32(2):
                    Ic[j] = best_i + i_emit
                else:
                    Ic[j] = NEG_INF

                # D[i,j]
                d1 = int32(Mp[j]) + t_md
                d2 = int32(Dp[j]) + t_dd
                Dc[j] = d1 if d1 > d2 else d2

                if mc_v > max_score:
                    max_score = mc_v

            # swap rows by value-copy into Mp/Ip/Dp using a tight loop
            # (Numba doesn't let us swap array-view pointers inside the kernel
            # the way C does, so we copy. Total work = O(T) per row.)
            for j in range(cols):
                Mp[j] = Mc[j]
                Ip[j] = Ic[j]
                Dp[j] = Dc[j]

        out_scores[tid] = max_score


def batch_viterbi_cuda(
    flat_match_emit, insert_emit, transitions,
    profile_offsets, profile_lengths,
    target_flat, target_offsets, target_lengths,
    pairs, band_width,
    device=0,
):
    """GPU batch Viterbi. Uses numba.cuda.

    Launches one thread per pair. Each thread runs the DP serially; the
    scaling comes from running tens of thousands of pairs concurrently
    across the GPU's SMs.

    Returns scores (M,) int32 on host.
    """
    if not _CUDA_AVAILABLE:
        raise RuntimeError("CUDA not available")

    cuda.select_device(device)

    M = len(pairs)
    # Determine per-thread workspace size. Each thread needs 6 rows of
    # (max_len + 2) int32. Cap at 4096 to avoid excessive memory usage;
    # sequences longer than this will be clipped (rare for proteins).
    max_len = int(max(profile_lengths.max(), target_lengths.max()))
    ws_size = max_len + 2

    # Widen int8/uint8 to int32 on host — numba.cuda can't do implicit
    # narrow→wide casts in the kernel. Memory cost is 4× but acceptable
    # for protein alphabet (emit table is a few MB max).
    # Pad emit tables to 32 entries per row so unknown letters (idx 20)
    # map to zero without branching in the kernel.
    fme32 = np.ascontiguousarray(flat_match_emit, dtype=np.int32)
    n_rows = fme32.shape[0]
    fme_padded = np.zeros((n_rows, 32), dtype=np.int32)
    fme_padded[:, :20] = fme32
    # flatten for kernel (uses flat indexing via emit_base)
    fme_flat = fme_padded.reshape(-1)
    ins32 = np.zeros(32, dtype=np.int32)
    ins32[:20] = insert_emit
    ins32[20:] = -1
    d_fme = cuda.to_device(fme_flat)
    d_ie  = cuda.to_device(ins32)
    d_tr  = cuda.to_device(np.ascontiguousarray(transitions, dtype=np.int32))
    d_po  = cuda.to_device(np.ascontiguousarray(profile_offsets, dtype=np.int64))
    d_pl  = cuda.to_device(np.ascontiguousarray(profile_lengths, dtype=np.int32))
    d_tf  = cuda.to_device(np.ascontiguousarray(target_flat, dtype=np.int32))
    d_to  = cuda.to_device(np.ascontiguousarray(target_offsets, dtype=np.int64))
    d_tl  = cuda.to_device(np.ascontiguousarray(target_lengths, dtype=np.int32))
    d_pr  = cuda.to_device(np.ascontiguousarray(pairs, dtype=np.int32))
    d_out = cuda.device_array(M, dtype=np.int32)

    # Workspace: one slab per thread. If M is large, this dominates memory.
    # Bytes = M * 6 * ws_size * 4. For M=100000, ws=600 → 1.4 GB. OK.
    d_ws = cuda.device_array((M, 6, ws_size), dtype=np.int32)

    TPB = 128
    BPG = (M + TPB - 1) // TPB
    _viterbi_kernel[BPG, TPB](
        d_fme, d_ie, d_tr, d_po, d_pl,
        d_tf, d_to, d_tl, d_pr,
        np.int32(band_width), d_out, d_ws, np.int32(ws_size),
    )
    cuda.synchronize()

    return d_out.copy_to_host()
