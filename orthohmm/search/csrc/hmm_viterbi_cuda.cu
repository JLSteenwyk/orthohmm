/*
 * Warp-collaborative CUDA Viterbi kernel for profile HMM search.
 *
 * Design: 1 pair per warp (32 threads). All 32 threads collaborate on
 * the DP for one (profile, target) pair. Shared memory holds the 6 DP
 * row buffers so accesses are 1-cycle instead of ~500-cycle (global).
 *
 * Per row i:
 *   1. All threads parallel-compute M[i,j] and D[i,j] for strided j
 *      (only depend on previous row — no intra-row dependency).
 *   2. Thread 0 serially computes I[i,j] for j = j_lo..j_hi
 *      (I has a same-row j-1 dependency; serial is simplest correct).
 *   3. Row swap via parallel copy (shared-memory pointers can't be
 *      rotated the way CPU pointers can).
 *
 * Shared memory per block: 6 * (T+2) * 4 bytes. For T=1998 fits in 48KB.
 * For longer targets the kernel falls back to a global-memory workspace
 * (not implemented in this version — caller must chunk).
 *
 * Build:
 *   nvcc -O3 -Xcompiler -fPIC -arch=sm_80 -shared
 *        -o hmm_viterbi_cuda.so hmm_viterbi_cuda.cu
 */

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

#define NEG_INF (-1000000)
#define WARP_SIZE 32
#define FULL_MASK 0xFFFFFFFFu


__global__ void hmm_viterbi_warp_kernel(
    const int8_t*  flat_match_emit,
    const int8_t*  insert_emit,
    const int32_t* transitions,
    const int64_t* prof_offsets,
    const int32_t* prof_lengths,
    const uint8_t* target_flat,
    const int64_t* target_offsets,
    const int32_t* target_lengths,
    const int32_t* pairs,
    int32_t        num_pairs,
    int32_t        band_width,
    int32_t        max_T,        /* pitch for shared-memory row stride */
    int32_t*       out_scores
) {
    int32_t pair_idx = blockIdx.x;
    if (pair_idx >= num_pairs) return;
    int32_t tid = threadIdx.x;  /* 0..31 */

    int32_t qi = pairs[pair_idx * 2];
    int32_t ti = pairs[pair_idx * 2 + 1];
    int32_t L = prof_lengths[qi];
    int32_t T = target_lengths[ti];
    if (L == 0 || T == 0) {
        if (tid == 0) out_scores[pair_idx] = 0;
        return;
    }

    int64_t p_off = prof_offsets[qi];
    int64_t t_off = target_offsets[ti];

    /* Shared memory: 6 rows of size (max_T+2). We use max_T as stride
     * so all blocks in a launch have identical layout regardless of
     * the individual pair's T. */
    extern __shared__ int32_t shared_ws[];
    int32_t row_stride = max_T + 2;
    int32_t* Mp = shared_ws + 0 * row_stride;
    int32_t* Mc = shared_ws + 1 * row_stride;
    int32_t* Ip = shared_ws + 2 * row_stride;
    int32_t* Ic = shared_ws + 3 * row_stride;
    int32_t* Dp = shared_ws + 4 * row_stride;
    int32_t* Dc = shared_ws + 5 * row_stride;

    /* Initialize all 6 rows to NEG_INF (parallel) */
    for (int32_t j = tid; j <= T + 1; j += WARP_SIZE) {
        Mp[j] = NEG_INF; Mc[j] = NEG_INF;
        Ip[j] = NEG_INF; Ic[j] = NEG_INF;
        Dp[j] = NEG_INF; Dc[j] = NEG_INF;
    }
    __syncwarp();

    const int32_t t_mm = transitions[0];
    const int32_t t_mi = transitions[1];
    const int32_t t_md = transitions[2];
    const int32_t t_im = transitions[3];
    const int32_t t_ii = transitions[4];
    const int32_t t_dm = transitions[5];
    const int32_t t_dd = transitions[6];

    int32_t max_score = 0;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;

    for (int32_t i = 1; i <= L; i++) {
        int32_t j_center = (int32_t)(((int64_t)i * T) / L);
        int32_t j_lo = j_center - bw;
        int32_t j_hi = j_center + bw;
        if (j_lo < 1) j_lo = 1;
        if (j_hi > T) j_hi = T;

        /* Clear current-row borders in band */
        int32_t clear_start = (j_lo - 1 > 0) ? j_lo - 1 : 0;
        int32_t clear_end = j_hi + 1;
        if (clear_end > T + 1) clear_end = T + 1;
        for (int32_t j = clear_start + tid; j <= clear_end; j += WARP_SIZE) {
            Mc[j] = NEG_INF; Ic[j] = NEG_INF; Dc[j] = NEG_INF;
        }
        __syncwarp();

        /* Parallel: compute M and D for strided j */
        const int8_t* emit_row = flat_match_emit + (p_off + i - 1) * 20;

        for (int32_t j = j_lo + tid; j <= j_hi; j += WARP_SIZE) {
            uint8_t t_aa = target_flat[t_off + j - 1];
            int32_t m_emit = (t_aa < 20) ? (int32_t)emit_row[t_aa] : 0;

            int32_t best_m = 0;
            int32_t v;
            v = Mp[j - 1] + t_mm; if (v > best_m) best_m = v;
            v = Ip[j - 1] + t_im; if (v > best_m) best_m = v;
            v = Dp[j - 1] + t_dm; if (v > best_m) best_m = v;
            int32_t mc_j = best_m + m_emit;
            Mc[j] = mc_j;
            if (mc_j > max_score) max_score = mc_j;

            int32_t d1 = Mp[j] + t_md;
            int32_t d2 = Dp[j] + t_dd;
            Dc[j] = (d1 > d2) ? d1 : d2;
        }
        __syncwarp();

        /* Parallel I-state scan via algebraic reformulation.
         *
         * Recurrence: I[j] = max(M[j-1]+tMI + iE[j], I[j-1]+tII + iE[j])
         * Let a[j] = M[j-1] + tMI + iE[j]   (leaf-value)
         *     c[j] = tII + iE[j]            (continuation cost)
         * Then I[j] = max over k in [1,j] of (a[k] + sum_{m=k+1..j} c[m])
         *          = a[j] IF best recent leaf; else a[k] + C[j] - C[k]
         * where C[j] = prefix sum of c[m] for m=1..j.
         * So I[j] = C[j] + max_prefix(a[k] - C[k]) over k=1..j.
         *
         * We chunk j into warp-sized (32-wide) blocks. Carry-in values
         * from the previous chunk: carry_C (offset) and carry_best
         * (running max prefix in the shifted space).                     */
        {
            int32_t chunk_lo = j_lo;
            int32_t carry_best = NEG_INF;   /* max (a[k]-C[k]) from previous chunks */
            int32_t carry_C = 0;             /* prefix-sum of c up to start of this chunk */
            /* All threads read j_lo..j_hi in chunks of 32 */
            while (chunk_lo <= j_hi) {
                int32_t j_t = chunk_lo + tid;   /* per-thread j in this chunk */
                int32_t in_chunk = (j_t <= j_hi);

                /* Per-thread a[j] and c[j] */
                int32_t aj = NEG_INF, cj = 0;
                if (in_chunk) {
                    uint8_t t_aa = target_flat[t_off + j_t - 1];
                    int32_t i_emit = (t_aa < 20) ? (int32_t)insert_emit[t_aa] : -1;
                    int32_t m_j1 = Mc[j_t - 1];
                    /* a[j] = Mc[j-1] + tMI + iE; if Mc[j-1] is NEG_INF, a = NEG_INF */
                    if (m_j1 > NEG_INF / 2) aj = m_j1 + t_mi + i_emit;
                    cj = t_ii + i_emit;
                }

                /* Inclusive prefix sum of cj over the warp (only the in-chunk
                 * threads' values — out-of-chunk threads contribute 0). */
                int32_t c_prefix = cj;
                #pragma unroll
                for (int32_t off = 1; off < WARP_SIZE; off *= 2) {
                    int32_t up = __shfl_up_sync(FULL_MASK, c_prefix, off);
                    if (tid >= off) c_prefix += up;
                }
                /* absolute prefix including prior chunks */
                int32_t C_j = carry_C + c_prefix;

                /* Compute B = a[j] - C[j]. For out-of-chunk threads use
                 * NEG_INF so they don't affect the max. */
                int32_t B = (in_chunk && aj > NEG_INF / 2) ? (aj - C_j) : NEG_INF;

                /* Inclusive prefix-max of B over the warp */
                int32_t B_max = B;
                #pragma unroll
                for (int32_t off = 1; off < WARP_SIZE; off *= 2) {
                    int32_t up = __shfl_up_sync(FULL_MASK, B_max, off);
                    if (tid >= off && up > B_max) B_max = up;
                }
                /* Combine with carry */
                int32_t best_at_j = (carry_best > B_max) ? carry_best : B_max;

                /* I[j] = C[j] + best_at_j */
                int32_t i_j = (best_at_j > NEG_INF / 2) ? (C_j + best_at_j) : NEG_INF;
                if (in_chunk) Ic[j_t] = i_j;

                /* Update carries for next chunk: broadcast thread 31's values */
                int32_t last_C    = __shfl_sync(FULL_MASK, c_prefix, WARP_SIZE - 1);
                int32_t last_bmax = __shfl_sync(FULL_MASK, B_max,   WARP_SIZE - 1);
                carry_C = carry_C + last_C;
                if (last_bmax > carry_best) carry_best = last_bmax;

                chunk_lo += WARP_SIZE;
            }
        }
        __syncwarp();

        /* Swap rows by parallel copy (can't rotate shared-mem pointers) */
        for (int32_t j = tid; j <= T + 1; j += WARP_SIZE) {
            Mp[j] = Mc[j];
            Ip[j] = Ic[j];
            Dp[j] = Dc[j];
        }
        __syncwarp();
    }

    /* Warp reduce max_score */
    for (int32_t offset = WARP_SIZE / 2; offset > 0; offset /= 2) {
        int32_t other = __shfl_xor_sync(FULL_MASK, max_score, offset);
        if (other > max_score) max_score = other;
    }

    if (tid == 0) {
        out_scores[pair_idx] = max_score;
    }
}


extern "C" {

/* Simple host wrapper that handles D2H/H2D and launch.
 * Caller passes already-on-device pointers (or uses the alloc variant). */
int hmm_viterbi_cuda_launch(
    const int8_t*  d_flat_match_emit,
    const int8_t*  d_insert_emit,
    const int32_t* d_transitions,
    const int64_t* d_prof_offsets,
    const int32_t* d_prof_lengths,
    const uint8_t* d_target_flat,
    const int64_t* d_target_offsets,
    const int32_t* d_target_lengths,
    const int32_t* d_pairs,
    int32_t        num_pairs,
    int32_t        band_width,
    int32_t        max_T,
    int32_t*       d_out_scores
) {
    size_t shared_bytes = (size_t)6 * (max_T + 2) * sizeof(int32_t);
    /* Ampere / Ada: up to 163 KB / 99 KB dynamic shared per block */
    if (shared_bytes > 99 * 1024) return -1;  /* too big — caller must chunk */

    cudaError_t err = cudaFuncSetAttribute(
        hmm_viterbi_warp_kernel,
        cudaFuncAttributeMaxDynamicSharedMemorySize,
        (int)shared_bytes);
    if (err != cudaSuccess) return -2;

    int32_t block_size = WARP_SIZE;
    int32_t grid_size = num_pairs;

    hmm_viterbi_warp_kernel<<<grid_size, block_size, shared_bytes>>>(
        d_flat_match_emit, d_insert_emit, d_transitions,
        d_prof_offsets, d_prof_lengths,
        d_target_flat, d_target_offsets, d_target_lengths,
        d_pairs, num_pairs, band_width, max_T,
        d_out_scores);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA kernel launch error: %s\n", cudaGetErrorString(err));
        return -3;
    }
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA kernel exec error: %s\n", cudaGetErrorString(err));
        return -4;
    }
    return 0;
}


/* All-in-one host-pointer interface: upload, compute, download. */
int hmm_viterbi_cuda_run(
    const int8_t*  h_flat_match_emit,   int64_t fme_n,
    const int8_t*  h_insert_emit,
    const int32_t* h_transitions,
    const int64_t* h_prof_offsets,      int32_t n_profs,
    const int32_t* h_prof_lengths,
    const uint8_t* h_target_flat,       int64_t tf_n,
    const int64_t* h_target_offsets,    int32_t n_tgts,
    const int32_t* h_target_lengths,
    const int32_t* h_pairs,             int32_t num_pairs,
    int32_t        band_width,
    int32_t        max_T,
    int32_t*       h_out_scores
) {
    int8_t*  d_fme = NULL;
    int8_t*  d_ie = NULL;
    int32_t* d_tr = NULL;
    int64_t* d_po = NULL;
    int32_t* d_pl = NULL;
    uint8_t* d_tf = NULL;
    int64_t* d_to = NULL;
    int32_t* d_tl = NULL;
    int32_t* d_pr = NULL;
    int32_t* d_out = NULL;
    int rc = -99;

#define CHECK(x) do { cudaError_t _err = (x); if (_err != cudaSuccess) { \
    fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(_err)); \
    rc = -100; goto cleanup; } } while (0)

    CHECK(cudaMalloc(&d_fme, fme_n));
    CHECK(cudaMalloc(&d_ie, 20));
    CHECK(cudaMalloc(&d_tr, 7 * 4));
    CHECK(cudaMalloc(&d_po, n_profs * 8));
    CHECK(cudaMalloc(&d_pl, n_profs * 4));
    CHECK(cudaMalloc(&d_tf, tf_n));
    CHECK(cudaMalloc(&d_to, n_tgts * 8));
    CHECK(cudaMalloc(&d_tl, n_tgts * 4));
    CHECK(cudaMalloc(&d_pr, num_pairs * 2 * 4));
    CHECK(cudaMalloc(&d_out, num_pairs * 4));

    CHECK(cudaMemcpy(d_fme, h_flat_match_emit, fme_n, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_ie, h_insert_emit, 20, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_tr, h_transitions, 7*4, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_po, h_prof_offsets, n_profs*8, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_pl, h_prof_lengths, n_profs*4, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_tf, h_target_flat, tf_n, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_to, h_target_offsets, n_tgts*8, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_tl, h_target_lengths, n_tgts*4, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_pr, h_pairs, num_pairs*2*4, cudaMemcpyHostToDevice));

    rc = hmm_viterbi_cuda_launch(
        d_fme, d_ie, d_tr, d_po, d_pl, d_tf, d_to, d_tl, d_pr,
        num_pairs, band_width, max_T, d_out);

    if (rc == 0) {
        CHECK(cudaMemcpy(h_out_scores, d_out, num_pairs*4, cudaMemcpyDeviceToHost));
    }

cleanup:
    if (d_fme) cudaFree(d_fme);
    if (d_ie)  cudaFree(d_ie);
    if (d_tr)  cudaFree(d_tr);
    if (d_po)  cudaFree(d_po);
    if (d_pl)  cudaFree(d_pl);
    if (d_tf)  cudaFree(d_tf);
    if (d_to)  cudaFree(d_to);
    if (d_tl)  cudaFree(d_tl);
    if (d_pr)  cudaFree(d_pr);
    if (d_out) cudaFree(d_out);
    return rc;
#undef CHECK
}


int hmm_viterbi_cuda_device_count(void) {
    int n = 0;
    cudaGetDeviceCount(&n);
    return n;
}

}  /* extern "C" */
