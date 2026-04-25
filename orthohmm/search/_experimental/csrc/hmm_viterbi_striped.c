/*
 * Farrar-striped 3-state profile HMM Viterbi with lazy-D correction.
 *
 * Profile is stored in striped layout: position i = lane v * Q + seg p
 * (with Q = ceil(L/V), V = 8 for AVX2 int32 lanes). At segment position
 * p, the 8-lane SIMD register holds values for profile positions
 * {p, p+Q, p+2Q, ..., p+7Q} — widely spaced.
 *
 * Recurrences (our convention):
 *   M[i,j] = emit_M + max(Mp[i-1,j-1]+tMM, Ip[i-1,j-1]+tIM,
 *                         Dp[i-1,j-1]+tDM, 0)
 *   I[i,j] = emit_I + max(Mp[i,j-1]+tMI, Ip[i,j-1]+tII)
 *   D[i,j] = max(M[i-1,j]+tMD, D[i-1,j]+tDD)
 *
 * M and I use j-1 (previous target) data — finalized, no lazy issue.
 * D uses CURRENT-j's previous profile position, which wraps at seg 0.
 * Farrar's lazy-F (here lazy-D): assume NEG_INF wrap, run all segments,
 * check wrap, re-run if it propagates. Converges in 1-2 iterations
 * in practice.
 *
 * UNBANDED. Banding with striping is complex; see multipair kernel for
 * banded SIMD.
 *
 * Build: gcc -O3 -march=native -mavx2 -fopenmp -shared -fPIC
 *          -o hmm_viterbi_striped.so hmm_viterbi_striped.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <immintrin.h>

#define V 8                        /* AVX2 int32 lane count */
#define NEG_INF (-1000000)

void striped_set_num_threads(int32_t n) { omp_set_num_threads(n); }


/* Shift-right lanes by 1: [v0,v1,...,v7] → [X, v0, v1, ..., v6].
 * X is the value from lane 7 (wrap-around), caller must overwrite
 * lane 0 if a different value is wanted. */
static inline __m256i shift_right_1_lane(__m256i v) {
    /* Use permutevar8x32 with indices [7,0,1,2,3,4,5,6]: result[l] = v[idx[l]].
     * So result[0] = v[7] (wrap), result[1] = v[0], ..., result[7] = v[6]. */
    const __m256i idx = _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6);
    return _mm256_permutevar8x32_epi32(v, idx);
}

/* Scalar-pack ordinary (L, 20) int8 profile into striped int32 layout.
 * Output: emit_striped[p * 21 * V + a * V + v] where
 *   lane v, seg p holds the value for profile position v*Q + p, amino acid a.
 * We pad to 21 AA slots per seg so index a=20 (unknown/X) safely reads
 * as zero rather than spilling into the next segment's table. */
static void build_striped_profile(
    const int8_t* match_emit,  /* (L, 20) int8 row-major */
    int32_t L, int32_t Q,
    int32_t* emit_striped       /* (Q, 21, V) int32 — V=8, slot 20 is zeros */
) {
    for (int32_t p = 0; p < Q; p++) {
        for (int32_t a = 0; a < 20; a++) {
            for (int32_t v = 0; v < V; v++) {
                int32_t i = v * Q + p;
                emit_striped[p * 21 * V + a * V + v] =
                    (i < L) ? (int32_t)match_emit[i * 20 + a] : 0;
            }
        }
        /* Unknown-AA slot: emit 0 for all lanes */
        for (int32_t v = 0; v < V; v++) {
            emit_striped[p * 21 * V + 20 * V + v] = 0;
        }
    }
}


/* Core striped Viterbi for one pair. Returns max score.
 * Workspace: 6 __m256i arrays of size Q each. */
int32_t hmm_viterbi_striped_one(
    const int8_t*   match_emit,       /* (L, 20) */
    const int8_t*   insert_emit,      /* (20,) */
    const int32_t*  transitions,      /* (7,) */
    int32_t         L,
    const uint8_t*  target,
    int32_t         T,
    int32_t         band_width,       /* 0 or <=0 = full DP */
    int32_t*        emit_striped_ws,  /* (Q, 21, V) int32 scratch */
    __m256i*        Mp, __m256i* Mc,
    __m256i*        Ip, __m256i* Ic,
    __m256i*        Dp, __m256i* Dc,
    int32_t         Q                  /* = ceil(L/V) */
) {
    if (L == 0 || T == 0) return 0;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    int32_t banded = (bw > 0 && max_dim > 50);

    /* Build striped profile into scratch */
    build_striped_profile(match_emit, L, Q, emit_striped_ws);

    /* Widen insert_emit to int32[32] (pad to 32 for easy gather index) */
    int32_t i_emit_i32[32];
    for (int32_t a = 0; a < 20; a++) i_emit_i32[a] = (int32_t)insert_emit[a];
    for (int32_t a = 20; a < 32; a++) i_emit_i32[a] = -1;

    const __m256i v_tmm = _mm256_set1_epi32(transitions[0]);
    const __m256i v_tmi = _mm256_set1_epi32(transitions[1]);
    const __m256i v_tmd = _mm256_set1_epi32(transitions[2]);
    const __m256i v_tim = _mm256_set1_epi32(transitions[3]);
    const __m256i v_tii = _mm256_set1_epi32(transitions[4]);
    const __m256i v_tdm = _mm256_set1_epi32(transitions[5]);
    const __m256i v_tdd = _mm256_set1_epi32(transitions[6]);
    const __m256i v_neg_inf = _mm256_set1_epi32(NEG_INF);
    const __m256i v_zero = _mm256_setzero_si256();

    /* Initialize: all states NEG_INF */
    for (int32_t p = 0; p < Q; p++) {
        Mp[p] = Ip[p] = Dp[p] = v_neg_inf;
    }

    __m256i v_max = v_zero;

    for (int32_t j = 1; j <= T; j++) {
        uint8_t ta = target[j - 1];
        int32_t t_aa = (ta < 20) ? (int32_t)ta : 20;

        /* Wraps for M state (from FINALIZED Mp, Ip, Dp — previous j) */
        __m256i M_wrap_prev = shift_right_1_lane(Mp[Q - 1]);
        __m256i I_wrap_prev = shift_right_1_lane(Ip[Q - 1]);
        __m256i D_wrap_prev = shift_right_1_lane(Dp[Q - 1]);
        /* Lane 0 of these wraps represents i=-1 (invalid); set NEG_INF. */
        M_wrap_prev = _mm256_blend_epi32(M_wrap_prev, v_neg_inf, 0x01);
        I_wrap_prev = _mm256_blend_epi32(I_wrap_prev, v_neg_inf, 0x01);
        D_wrap_prev = _mm256_blend_epi32(D_wrap_prev, v_neg_inf, 0x01);

        /* Lazy wrap for D's i-1 within CURRENT j */
        __m256i M_wrap_curr = v_neg_inf;
        __m256i D_wrap_curr = v_neg_inf;

        /* Insert emission for this target residue (broadcast) */
        __m256i i_emit = _mm256_set1_epi32(i_emit_i32[t_aa]);

        /* Per-lane band bounds. Profile positions in band for this j:
         *   (j-bw) * L / T <= i <= (j+bw) * L / T
         * Lane v at seg p holds position v*Q + p. In-band seg range per
         * lane: [max(0, i_lo - v*Q), min(Q-1, i_hi - v*Q)]. */
        int32_t p_union_lo = 0, p_union_hi = Q - 1;
        int32_t plo_per_lane[V] __attribute__((aligned(32)));
        int32_t phi_per_lane[V] __attribute__((aligned(32)));
        __m256i v_active_lane = _mm256_set1_epi32(-1);  /* default: all active */

        if (banded) {
            int64_t i_lo_64 = (int64_t)(j - bw) * L / T;
            int64_t i_hi_64 = (int64_t)(j + bw) * L / T;
            if (i_lo_64 < 0) i_lo_64 = 0;
            if (i_hi_64 >= L) i_hi_64 = L - 1;
            int32_t i_lo = (int32_t)i_lo_64;
            int32_t i_hi = (int32_t)i_hi_64;

            p_union_lo = Q;
            p_union_hi = -1;
            int32_t active_mask[V] __attribute__((aligned(32)));
            for (int32_t v = 0; v < V; v++) {
                int32_t p_lo = i_lo - v * Q;
                int32_t p_hi = i_hi - v * Q;
                if (p_lo < 0) p_lo = 0;
                if (p_hi > Q - 1) p_hi = Q - 1;
                if (p_lo <= p_hi) {
                    plo_per_lane[v] = p_lo;
                    phi_per_lane[v] = p_hi;
                    active_mask[v] = -1;
                    if (p_lo < p_union_lo) p_union_lo = p_lo;
                    if (p_hi > p_union_hi) p_union_hi = p_hi;
                } else {
                    plo_per_lane[v] = Q;       /* impossible condition */
                    phi_per_lane[v] = -1;
                    active_mask[v] = 0;
                }
            }
            if (p_union_lo > p_union_hi) {
                /* No lane in band — skip entire j */
                continue;
            }
            v_active_lane = _mm256_load_si256((const __m256i*)active_mask);

            /* Reset out-of-union segs to NEG_INF so stale data from prior
             * rows doesn't leak into next j's wrap calculations. */
            for (int32_t p = 0; p < p_union_lo; p++) {
                Mc[p] = v_neg_inf; Ic[p] = v_neg_inf; Dc[p] = v_neg_inf;
            }
            for (int32_t p = p_union_hi + 1; p < Q; p++) {
                Mc[p] = v_neg_inf; Ic[p] = v_neg_inf; Dc[p] = v_neg_inf;
            }
        }

        __m256i v_plo_lane = _mm256_load_si256((const __m256i*)plo_per_lane);
        __m256i v_phi_lane = _mm256_load_si256((const __m256i*)phi_per_lane);
        (void)v_plo_lane; (void)v_phi_lane;  /* used below */

        int32_t correction_passes = 0;
        const int32_t max_corrections = 16;

        while (1) {
            for (int32_t p = p_union_lo; p <= p_union_hi; p++) {
                /* Per-lane in-band mask for this p: plo[v] <= p <= phi[v] */
                __m256i v_p = _mm256_set1_epi32(p);
                __m256i in_lo = _mm256_cmpgt_epi32(
                    _mm256_add_epi32(v_p, _mm256_set1_epi32(1)), v_plo_lane);
                __m256i in_hi = _mm256_cmpgt_epi32(
                    _mm256_add_epi32(v_phi_lane, _mm256_set1_epi32(1)), v_p);
                __m256i in_band = _mm256_and_si256(
                    _mm256_and_si256(in_lo, in_hi), v_active_lane);
                /* Load 8 emission values for (seg p, t_aa) in one aligned
                 * load — that's our striping's payoff. Table stride is
                 * 21 AA slots so unknown (t_aa=20) reads the zero row. */
                __m256i m_emit = _mm256_load_si256(
                    (const __m256i*)(emit_striped_ws + p * 21 * V + t_aa * V));

                /* M sources: i-1 in previous j, with wrap at seg 0 */
                __m256i Mp_src = (p > 0) ? Mp[p - 1] : M_wrap_prev;
                __m256i Ip_src = (p > 0) ? Ip[p - 1] : I_wrap_prev;
                __m256i Dp_src = (p > 0) ? Dp[p - 1] : D_wrap_prev;

                __m256i best_m = _mm256_max_epi32(
                    _mm256_max_epi32(
                        _mm256_add_epi32(Mp_src, v_tmm),
                        _mm256_add_epi32(Ip_src, v_tim)),
                    _mm256_max_epi32(
                        _mm256_add_epi32(Dp_src, v_tdm),
                        v_zero));
                __m256i mc_raw = _mm256_add_epi32(best_m, m_emit);
                Mc[p] = _mm256_blendv_epi8(v_neg_inf, mc_raw, in_band);
                v_max = _mm256_max_epi32(v_max, Mc[p]);

                /* I state: same profile pos, previous j */
                __m256i best_i = _mm256_max_epi32(
                    _mm256_add_epi32(Mp[p], v_tmi),
                    _mm256_add_epi32(Ip[p], v_tii));
                __m256i i_valid = _mm256_cmpgt_epi32(best_i,
                    _mm256_set1_epi32(NEG_INF / 2));
                __m256i ic_raw = _mm256_blendv_epi8(v_neg_inf,
                    _mm256_add_epi32(best_i, i_emit), i_valid);
                Ic[p] = _mm256_blendv_epi8(v_neg_inf, ic_raw, in_band);

                /* D state: uses CURRENT j's i-1 values */
                __m256i Mc_src = (p > 0) ? Mc[p - 1] : M_wrap_curr;
                __m256i Dc_src = (p > 0) ? Dc[p - 1] : D_wrap_curr;
                __m256i dc_raw = _mm256_max_epi32(
                    _mm256_add_epi32(Mc_src, v_tmd),
                    _mm256_add_epi32(Dc_src, v_tdd));
                Dc[p] = _mm256_blendv_epi8(v_neg_inf, dc_raw, in_band);
            }

            /* Lazy-D correction check: what wrap values DID we just produce? */
            __m256i new_M_wrap = shift_right_1_lane(Mc[Q - 1]);
            __m256i new_D_wrap = shift_right_1_lane(Dc[Q - 1]);
            new_M_wrap = _mm256_blend_epi32(new_M_wrap, v_neg_inf, 0x01);
            new_D_wrap = _mm256_blend_epi32(new_D_wrap, v_neg_inf, 0x01);

            /* Check if new wraps would change Dc[0]. If the actual
             * wrap contribution is larger than what we assumed
             * (NEG_INF on first pass, or the last iteration's wrap),
             * we need to re-run with the updated wrap. */
            __m256i actual_dc0 = _mm256_max_epi32(
                _mm256_add_epi32(new_M_wrap, v_tmd),
                _mm256_add_epi32(new_D_wrap, v_tdd));
            __m256i have_dc0 = Dc[0];

            /* Did the actual wrap produce a bigger Dc[0]? */
            __m256i bigger = _mm256_cmpgt_epi32(actual_dc0, have_dc0);
            int32_t mask = _mm256_movemask_epi8(bigger);
            if (mask == 0 || correction_passes >= max_corrections) break;

            /* Need correction: update wraps and re-run */
            M_wrap_curr = new_M_wrap;
            D_wrap_curr = new_D_wrap;
            correction_passes++;
        }

        /* Row swap: current becomes previous */
        __m256i* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    /* Horizontal max of v_max */
    int32_t scores[V] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*)scores, v_max);
    int32_t best = 0;
    for (int32_t l = 0; l < V; l++) if (scores[l] > best) best = scores[l];
    return best;
}


/* Batch driver — 1 pair per thread, no SIMD across pairs (each pair is
 * itself SIMD-striped). */
void batch_hmm_viterbi_striped_c(
    const int8_t*   flat_match_emit,
    const int8_t*   insert_emit,
    const int32_t*  transitions,
    const int64_t*  prof_offsets,
    const int32_t*  prof_lengths,
    const uint8_t*  target_flat,
    const int64_t*  target_offsets,
    const int32_t*  target_lengths,
    const int32_t*  pairs,
    int32_t         M,
    int32_t         band_width,
    int32_t*        out_scores
) {

    /* Find max Q for workspace sizing */
    int32_t max_Q = 0;
    int32_t max_L = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t L = prof_lengths[pairs[p * 2]];
        if (L > max_L) max_L = L;
    }
    max_Q = (max_L + V - 1) / V;

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Per-thread workspace: 6 __m256i arrays of size max_Q + 1 striped emit table */
    __m256i** ws = (__m256i**)calloc(n_threads_actual * 6, sizeof(__m256i*));
    int32_t** emit_ws = (int32_t**)calloc(n_threads_actual, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            ws[t * 6 + s] = (__m256i*)aligned_alloc(64, (max_Q + 1) * sizeof(__m256i));
        }
        emit_ws[t] = (int32_t*)aligned_alloc(64, max_Q * 21 * V * sizeof(int32_t));
    }

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        __m256i* Mp = ws[tid * 6 + 0];
        __m256i* Mc = ws[tid * 6 + 1];
        __m256i* Ip = ws[tid * 6 + 2];
        __m256i* Ic = ws[tid * 6 + 3];
        __m256i* Dp = ws[tid * 6 + 4];
        __m256i* Dc = ws[tid * 6 + 5];
        int32_t* e_ws = emit_ws[tid];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];
            int32_t L = prof_lengths[qi];
            int32_t T = target_lengths[ti];
            int32_t Q = (L + V - 1) / V;
            if (Q == 0) { out_scores[p] = 0; continue; }

            out_scores[p] = hmm_viterbi_striped_one(
                flat_match_emit + prof_offsets[qi] * 20,
                insert_emit, transitions, L,
                target_flat + target_offsets[ti], T, band_width,
                e_ws, Mp, Mc, Ip, Ic, Dp, Dc, Q);
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) free(ws[t * 6 + s]);
        free(emit_ws[t]);
    }
    free(ws); free(emit_ws);
}
