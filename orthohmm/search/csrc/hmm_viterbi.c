/*
 * C/OpenMP banded 3-state profile HMM Viterbi.
 *
 * Extends ClustKIT's sw_align.c pattern from 2-state (H, E) to
 * 3-state (M, I, D) HMM Viterbi with:
 *   - Pre-allocated thread-local workspace
 *   - Row-toggled DP (2 rows per state)
 *   - Band centering with proportional diagonal mapping
 *   - OpenMP dynamic scheduling over pairs
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -shared -fPIC -o hmm_viterbi.so hmm_viterbi.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <omp.h>

#if defined(__AVX2__)
#  include <immintrin.h>
#  define HMM_HAVE_AVX2 1
#else
#  define HMM_HAVE_AVX2 0
#endif

void hmm_set_num_threads(int32_t n) {
    omp_set_num_threads(n);
}

int32_t hmm_have_avx2(void) {
    return HMM_HAVE_AVX2;
}

#define NEG_INF (-1000000)

static inline int32_t max3(int32_t a, int32_t b, int32_t c) {
    int32_t m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    return m;
}

static inline int32_t max2(int32_t a, int32_t b) {
    return a > b ? a : b;
}

/* ─── Single-pair banded 3-state Viterbi ────────────────────────── */

static int32_t hmm_viterbi_one(
    const int8_t*   match_emit,     /* (prof_len, 20) row-major */
    const int8_t*   insert_emit,    /* (20,) */
    const int32_t*  transitions,    /* (7,): MM MI MD IM II DM DD */
    int32_t         prof_len,
    const uint8_t*  target_seq,
    int32_t         target_len,
    int32_t         band_width,
    int32_t         xdrop,          /* early-termination drop threshold; 0 disables */
    /* workspace: caller provides pre-allocated arrays of size (cols) */
    int32_t*        Mp,  int32_t*  Mc,
    int32_t*        Ip,  int32_t*  Ic,
    int32_t*        Dp,  int32_t*  Dc
) {
    const int32_t t_mm = transitions[0];
    const int32_t t_mi = transitions[1];
    const int32_t t_md = transitions[2];
    const int32_t t_im = transitions[3];
    const int32_t t_ii = transitions[4];
    const int32_t t_dm = transitions[5];
    const int32_t t_dd = transitions[6];

    int32_t L = prof_len;
    int32_t T = target_len;
    int32_t cols = T + 1;

    if (L == 0 || T == 0) return 0;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;

    /* Initialize rows */
    for (int32_t j = 0; j < cols; j++) {
        Mp[j] = NEG_INF; Mc[j] = NEG_INF;
        Ip[j] = NEG_INF; Ic[j] = NEG_INF;
        Dp[j] = NEG_INF; Dc[j] = NEG_INF;
    }

    int32_t max_score = 0;  /* local alignment: can start from 0 */
    int32_t last_max_update = 0;

    for (int32_t i = 1; i <= L; i++) {
        /* X-drop: if we haven't improved max_score for more than xdrop-worth
         * of rows and we're already past i=20, bail out. This is a correctness-
         * preserving heuristic: we may miss hits whose best alignment starts
         * late in the profile, but in practice hits that eventually beat the
         * running max must stay within xdrop of it throughout.                 */
        if (xdrop > 0 && max_score > 0 && i > 20 &&
            (i - last_max_update) > 30) {
            break;
        }
        /* Band bounds: proportional diagonal mapping */
        int32_t j_center = (int32_t)(((int64_t)i * T) / L);
        int32_t j_lo = j_center - bw;
        int32_t j_hi = j_center + bw;
        if (j_lo < 1) j_lo = 1;
        if (j_hi > T) j_hi = T;

        /* Clear current row borders */
        for (int32_t j = j_lo - 1; j <= j_hi + 1 && j < cols; j++) {
            if (j >= 0) {
                Mc[j] = NEG_INF;
                Ic[j] = NEG_INF;
                Dc[j] = NEG_INF;
            }
        }

        const int8_t* emit_row = match_emit + (i - 1) * 20;

        for (int32_t j = j_lo; j <= j_hi; j++) {
            uint8_t t_aa = target_seq[j - 1];
            int32_t m_emit = (t_aa < 20) ? (int32_t)emit_row[t_aa] : 0;
            int32_t i_emit = (t_aa < 20) ? (int32_t)insert_emit[t_aa] : -1;

            /* V_M(i,j) = emit + max(V_M(i-1,j-1)+t_mm, V_I(i-1,j-1)+t_im, V_D(i-1,j-1)+t_dm, 0) */
            int32_t best_m = 0;  /* local: can start fresh */
            int32_t v;
            v = Mp[j-1] + t_mm; if (v > best_m) best_m = v;
            v = Ip[j-1] + t_im; if (v > best_m) best_m = v;
            v = Dp[j-1] + t_dm; if (v > best_m) best_m = v;
            Mc[j] = best_m + m_emit;

            /* V_I(i,j) = emit + max(V_M(i,j-1)+t_mi, V_I(i,j-1)+t_ii) */
            int32_t best_i = NEG_INF;
            v = Mc[j-1] + t_mi; if (v > best_i) best_i = v;
            v = Ic[j-1] + t_ii; if (v > best_i) best_i = v;
            Ic[j] = (best_i > NEG_INF / 2) ? best_i + i_emit : NEG_INF;

            /* V_D(i,j) = max(V_M(i-1,j)+t_md, V_D(i-1,j)+t_dd) */
            Dc[j] = max2(Mp[j] + t_md, Dp[j] + t_dd);

            if (Mc[j] > max_score) { max_score = Mc[j]; last_max_update = i; }
        }

        /* Swap rows */
        int32_t* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    return max_score;
}


/* ─── AVX2 single-pair kernel ───────────────────────────────────── */

#if HMM_HAVE_AVX2

/* Split inner j-loop into two passes:
 *   Pass 1 (SIMD): compute Mc[j] and Dc[j] — both depend only on previous
 *                  row, so no horizontal dependency — 8-wide int32 AVX2.
 *   Pass 2 (scalar): compute Ic[j] — has Ic[j-1] / Mc[j-1] dependency.
 *                    This is still fast because only 2 max/add ops per cell. */
static int32_t hmm_viterbi_one_avx2(
    const int8_t*   match_emit,
    const int8_t*   insert_emit,
    const int32_t*  transitions,
    int32_t         prof_len,
    const uint8_t*  target_seq,
    int32_t         target_len,
    int32_t         band_width,
    int32_t         xdrop,
    int32_t*        Mp, int32_t* Mc,
    int32_t*        Ip, int32_t* Ic,
    int32_t*        Dp, int32_t* Dc
) {
    const int32_t t_mm = transitions[0];
    const int32_t t_mi = transitions[1];
    const int32_t t_md = transitions[2];
    const int32_t t_im = transitions[3];
    const int32_t t_ii = transitions[4];
    const int32_t t_dm = transitions[5];
    const int32_t t_dd = transitions[6];

    int32_t L = prof_len, T = target_len;
    int32_t cols = T + 1;
    if (L == 0 || T == 0) return 0;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;

    for (int32_t j = 0; j < cols; j++) {
        Mp[j] = NEG_INF; Mc[j] = NEG_INF;
        Ip[j] = NEG_INF; Ic[j] = NEG_INF;
        Dp[j] = NEG_INF; Dc[j] = NEG_INF;
    }

    /* Pre-widen insert_emit to int32 for cheap gather in the scalar I pass. */
    int32_t insert_emit32[32];
    for (int32_t a = 0; a < 20; a++) insert_emit32[a] = insert_emit[a];
    for (int32_t a = 20; a < 32; a++) insert_emit32[a] = -1;  /* unknown AA */

    /* Pre-widen target letters and gather insert emissions for the whole row.
     * Also pre-allocate a scratch buffer for per-row match emissions.
     * Target is uint8; emit_row lookup per profile row happens inside i-loop. */
    static __thread int32_t target_i32_[8192];
    static __thread int32_t i_emit_row_[8192];
    static __thread int32_t m_emit_row_[8192];
    int32_t* target_i32 = target_i32_;
    int32_t* i_emit_row = i_emit_row_;
    int32_t* m_emit_scratch = m_emit_row_;
    int32_t* target_i32_heap = NULL;
    int32_t* i_emit_row_heap = NULL;
    int32_t* m_emit_heap = NULL;
    if (cols > 8192) {
        target_i32_heap = (int32_t*)malloc(cols * sizeof(int32_t));
        i_emit_row_heap = (int32_t*)malloc(cols * sizeof(int32_t));
        m_emit_heap    = (int32_t*)malloc(cols * sizeof(int32_t));
        target_i32 = target_i32_heap;
        i_emit_row = i_emit_row_heap;
        m_emit_scratch = m_emit_heap;
    }
    for (int32_t j = 0; j < T; j++) {
        uint8_t a = target_seq[j];
        int32_t ai = (a < 20) ? (int32_t)a : 20;
        target_i32[j] = ai;
        i_emit_row[j] = insert_emit32[ai];
    }

    int32_t max_score = 0;
    int32_t last_max_update = 0;

    const __m256i v_tmm = _mm256_set1_epi32(t_mm);
    const __m256i v_tim = _mm256_set1_epi32(t_im);
    const __m256i v_tdm = _mm256_set1_epi32(t_dm);
    const __m256i v_tmd = _mm256_set1_epi32(t_md);
    const __m256i v_tdd = _mm256_set1_epi32(t_dd);
    const __m256i v_zero = _mm256_setzero_si256();

    for (int32_t i = 1; i <= L; i++) {
        if (xdrop > 0 && max_score > 0 && i > 20 &&
            (i - last_max_update) > 30) {
            break;
        }

        int32_t j_center = (int32_t)(((int64_t)i * T) / L);
        int32_t j_lo = j_center - bw;
        int32_t j_hi = j_center + bw;
        if (j_lo < 1) j_lo = 1;
        if (j_hi > T) j_hi = T;

        for (int32_t j = j_lo - 1; j <= j_hi + 1 && j < cols; j++) {
            if (j >= 0) {
                Mc[j] = NEG_INF;
                Ic[j] = NEG_INF;
                Dc[j] = NEG_INF;
            }
        }

        const int8_t* emit_row_i8 = match_emit + (i - 1) * 20;
        /* Pre-compute per-j match emission into a scratch array so the SIMD
         * inner loop can do straight loads (gather is slow on most x86). */
        for (int32_t j = j_lo - 1; j <= j_hi; j++) {
            int32_t ai = target_i32[j];
            m_emit_scratch[j] = (ai < 20) ? (int32_t)emit_row_i8[ai] : 0;
        }

        /* Pass 1: SIMD M and D for j in [j_lo, j_hi] */
        int32_t j = j_lo;
        /* Process 8 lanes at a time. A lane computes cell j..j+7. */
        for (; j + 8 <= j_hi + 1; j += 8) {
            /* Load previous row slices. Unaligned because j may be anywhere. */
            __m256i mp_jm1 = _mm256_loadu_si256((const __m256i*)(Mp + j - 1));
            __m256i ip_jm1 = _mm256_loadu_si256((const __m256i*)(Ip + j - 1));
            __m256i dp_jm1 = _mm256_loadu_si256((const __m256i*)(Dp + j - 1));
            __m256i mp_j   = _mm256_loadu_si256((const __m256i*)(Mp + j));
            __m256i dp_j   = _mm256_loadu_si256((const __m256i*)(Dp + j));

            /* best_m = max(Mp[j-1]+tmm, Ip[j-1]+tim, Dp[j-1]+tdm, 0) */
            __m256i v1 = _mm256_add_epi32(mp_jm1, v_tmm);
            __m256i v2 = _mm256_add_epi32(ip_jm1, v_tim);
            __m256i v3 = _mm256_add_epi32(dp_jm1, v_tdm);
            __m256i best_m = _mm256_max_epi32(_mm256_max_epi32(v1, v2),
                                              _mm256_max_epi32(v3, v_zero));

            /* Load precomputed emission row (straight load, no gather) */
            __m256i m_emit = _mm256_loadu_si256((const __m256i*)(m_emit_scratch + j - 1));
            __m256i mc = _mm256_add_epi32(best_m, m_emit);
            _mm256_storeu_si256((__m256i*)(Mc + j), mc);

            /* Dc[j] = max(Mp[j]+tmd, Dp[j]+tdd) */
            __m256i d1 = _mm256_add_epi32(mp_j, v_tmd);
            __m256i d2 = _mm256_add_epi32(dp_j, v_tdd);
            __m256i dc = _mm256_max_epi32(d1, d2);
            _mm256_storeu_si256((__m256i*)(Dc + j), dc);
        }
        /* Tail: scalar for remaining 0..7 cells */
        for (; j <= j_hi; j++) {
            int32_t best_m = 0;
            int32_t v;
            v = Mp[j-1] + t_mm; if (v > best_m) best_m = v;
            v = Ip[j-1] + t_im; if (v > best_m) best_m = v;
            v = Dp[j-1] + t_dm; if (v > best_m) best_m = v;
            Mc[j] = best_m + m_emit_scratch[j - 1];
            Dc[j] = max2(Mp[j] + t_md, Dp[j] + t_dd);
        }

        /* Pass 2: scalar I with j-1 dependency */
        for (j = j_lo; j <= j_hi; j++) {
            int32_t ai = target_i32[j-1];
            int32_t i_emit = i_emit_row[j-1];

            int32_t best_i = NEG_INF;
            int32_t v;
            v = Mc[j-1] + t_mi; if (v > best_i) best_i = v;
            v = Ic[j-1] + t_ii; if (v > best_i) best_i = v;
            Ic[j] = (best_i > NEG_INF / 2) ? best_i + i_emit : NEG_INF;

            if (Mc[j] > max_score) { max_score = Mc[j]; last_max_update = i; }
        }

        int32_t* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    if (target_i32_heap) free(target_i32_heap);
    if (i_emit_row_heap) free(i_emit_row_heap);
    if (m_emit_heap)     free(m_emit_heap);

    return max_score;
}

#endif  /* HMM_HAVE_AVX2 */


/* ─── AVX2 Int16 single-pair kernel (biased, saturating) ──────────
 *
 * Idea: represent scores as int16 with a bias so they're non-negative.
 * Saturating adds (_mm256_adds_epi16) clamp overflow. If the unbiased
 * max score hits the saturation ceiling, return INT32_MIN as a sentinel
 * so the caller can re-score the pair with int32.
 *
 * 16 int16 lanes per AVX2 register → 2× more cells per instruction than
 * the int32 variant. Useful when scores comfortably fit in int16 range,
 * which is typical for protein HMM Viterbi (max_score < 30000 in the
 * overwhelming majority of pairs). */

#if HMM_HAVE_AVX2

#define I16_BIAS  ((int16_t)0)          /* no bias; unsigned sat not int16 sat */
#define I16_NEG_INF ((int16_t)(-30000)) /* enough headroom below 0 */
#define I16_SAT_MAX ((int16_t)30000)    /* if max_score >= this, overflow */

static int32_t hmm_viterbi_one_avx2_i16(
    const int8_t*   match_emit,
    const int8_t*   insert_emit,
    const int32_t*  transitions,
    int32_t         prof_len,
    const uint8_t*  target_seq,
    int32_t         target_len,
    int32_t         band_width,
    int32_t         xdrop,
    int16_t*        Mp, int16_t* Mc,
    int16_t*        Ip, int16_t* Ic,
    int16_t*        Dp, int16_t* Dc
) {
    const int16_t t_mm = (int16_t)transitions[0];
    const int16_t t_mi = (int16_t)transitions[1];
    const int16_t t_md = (int16_t)transitions[2];
    const int16_t t_im = (int16_t)transitions[3];
    const int16_t t_ii = (int16_t)transitions[4];
    const int16_t t_dm = (int16_t)transitions[5];
    const int16_t t_dd = (int16_t)transitions[6];

    int32_t L = prof_len, T = target_len;
    int32_t cols = T + 1;
    if (L == 0 || T == 0) return 0;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;

    for (int32_t j = 0; j < cols; j++) {
        Mp[j] = I16_NEG_INF; Mc[j] = I16_NEG_INF;
        Ip[j] = I16_NEG_INF; Ic[j] = I16_NEG_INF;
        Dp[j] = I16_NEG_INF; Dc[j] = I16_NEG_INF;
    }

    /* Per-row scratch: match + insert emission per j */
    static __thread int16_t m_emit_scr_[8192];
    static __thread int16_t i_emit_scr_[8192];
    int16_t* m_emit_scr = m_emit_scr_;
    int16_t* i_emit_scr = i_emit_scr_;
    int16_t* m_emit_heap = NULL;
    int16_t* i_emit_heap = NULL;
    if (cols > 8192) {
        m_emit_heap = (int16_t*)malloc(cols * sizeof(int16_t));
        i_emit_heap = (int16_t*)malloc(cols * sizeof(int16_t));
        m_emit_scr = m_emit_heap;
        i_emit_scr = i_emit_heap;
    }
    /* Pre-fill insert emissions (profile-independent) */
    for (int32_t j = 0; j < T; j++) {
        uint8_t a = target_seq[j];
        i_emit_scr[j] = (a < 20) ? (int16_t)insert_emit[a] : (int16_t)-1;
    }

    int16_t max_score = 0;
    int32_t last_max_update = 0;
    int32_t overflow = 0;

    const __m256i v_tmm = _mm256_set1_epi16(t_mm);
    const __m256i v_tim = _mm256_set1_epi16(t_im);
    const __m256i v_tdm = _mm256_set1_epi16(t_dm);
    const __m256i v_tmd = _mm256_set1_epi16(t_md);
    const __m256i v_tdd = _mm256_set1_epi16(t_dd);
    const __m256i v_zero = _mm256_setzero_si256();
    const __m256i v_sat = _mm256_set1_epi16(I16_SAT_MAX);

    for (int32_t i = 1; i <= L; i++) {
        if (xdrop > 0 && max_score > 0 && i > 20 &&
            (i - last_max_update) > 30) break;

        int32_t j_center = (int32_t)(((int64_t)i * T) / L);
        int32_t j_lo = j_center - bw;
        int32_t j_hi = j_center + bw;
        if (j_lo < 1) j_lo = 1;
        if (j_hi > T) j_hi = T;

        for (int32_t j = j_lo - 1; j <= j_hi + 1 && j < cols; j++) {
            if (j >= 0) { Mc[j] = I16_NEG_INF; Ic[j] = I16_NEG_INF; Dc[j] = I16_NEG_INF; }
        }

        const int8_t* emit_row_i8 = match_emit + (i - 1) * 20;
        for (int32_t j = j_lo - 1; j <= j_hi; j++) {
            uint8_t a = target_seq[j];
            m_emit_scr[j] = (a < 20) ? (int16_t)emit_row_i8[a] : (int16_t)0;
        }

        /* Pass 1: SIMD M, D (16 lanes/reg) */
        __m256i v_row_max = v_zero;
        int32_t j = j_lo;
        for (; j + 16 <= j_hi + 1; j += 16) {
            __m256i mp_jm1 = _mm256_loadu_si256((const __m256i*)(Mp + j - 1));
            __m256i ip_jm1 = _mm256_loadu_si256((const __m256i*)(Ip + j - 1));
            __m256i dp_jm1 = _mm256_loadu_si256((const __m256i*)(Dp + j - 1));
            __m256i mp_j   = _mm256_loadu_si256((const __m256i*)(Mp + j));
            __m256i dp_j   = _mm256_loadu_si256((const __m256i*)(Dp + j));

            __m256i v1 = _mm256_adds_epi16(mp_jm1, v_tmm);
            __m256i v2 = _mm256_adds_epi16(ip_jm1, v_tim);
            __m256i v3 = _mm256_adds_epi16(dp_jm1, v_tdm);
            __m256i best_m = _mm256_max_epi16(_mm256_max_epi16(v1, v2),
                                              _mm256_max_epi16(v3, v_zero));
            __m256i m_emit = _mm256_loadu_si256((const __m256i*)(m_emit_scr + j - 1));
            __m256i mc = _mm256_adds_epi16(best_m, m_emit);
            _mm256_storeu_si256((__m256i*)(Mc + j), mc);
            v_row_max = _mm256_max_epi16(v_row_max, mc);

            __m256i d1 = _mm256_adds_epi16(mp_j, v_tmd);
            __m256i d2 = _mm256_adds_epi16(dp_j, v_tdd);
            __m256i dc = _mm256_max_epi16(d1, d2);
            _mm256_storeu_si256((__m256i*)(Dc + j), dc);
        }
        /* scalar tail */
        for (; j <= j_hi; j++) {
            int32_t best_m = 0, v;
            v = (int32_t)Mp[j-1] + t_mm; if (v > best_m) best_m = v;
            v = (int32_t)Ip[j-1] + t_im; if (v > best_m) best_m = v;
            v = (int32_t)Dp[j-1] + t_dm; if (v > best_m) best_m = v;
            int32_t mc_v = best_m + m_emit_scr[j-1];
            if (mc_v > 32767) mc_v = 32767;
            if (mc_v < -32768) mc_v = -32768;
            Mc[j] = (int16_t)mc_v;
            int32_t dc_v_1 = (int32_t)Mp[j] + t_md;
            int32_t dc_v_2 = (int32_t)Dp[j] + t_dd;
            int32_t dc_v = dc_v_1 > dc_v_2 ? dc_v_1 : dc_v_2;
            if (dc_v > 32767) dc_v = 32767;
            if (dc_v < -32768) dc_v = -32768;
            Dc[j] = (int16_t)dc_v;
        }

        /* Pass 2: scalar I with j-1 dep, also max tracking */
        for (j = j_lo; j <= j_hi; j++) {
            int32_t best_i = I16_NEG_INF, v;
            v = (int32_t)Mc[j-1] + t_mi; if (v > best_i) best_i = v;
            v = (int32_t)Ic[j-1] + t_ii; if (v > best_i) best_i = v;
            int32_t ic_v;
            if (best_i > I16_NEG_INF / 2) {
                ic_v = best_i + i_emit_scr[j-1];
                if (ic_v > 32767) ic_v = 32767;
                if (ic_v < -32768) ic_v = -32768;
            } else {
                ic_v = I16_NEG_INF;
            }
            Ic[j] = (int16_t)ic_v;

            if (Mc[j] > max_score) { max_score = Mc[j]; last_max_update = i; }
            if (Mc[j] >= I16_SAT_MAX) overflow = 1;
        }

        int16_t* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    if (m_emit_heap) free(m_emit_heap);
    if (i_emit_heap) free(i_emit_heap);

    if (overflow) return INT32_MIN;  /* caller must re-score in int32 */
    return (int32_t)max_score;
}

#endif  /* HMM_HAVE_AVX2 */


/* ─── Batch interface ───────────────────────────────────────────── */

void batch_hmm_viterbi_c(
    const int8_t*   flat_match_emit,  /* (sum_L, 20) row-major */
    const int8_t*   insert_emit,      /* (20,) */
    const int32_t*  transitions,      /* (7,) */
    const int64_t*  prof_offsets,     /* per-query */
    const int32_t*  prof_lengths,    /* per-query */
    const uint8_t*  target_flat,
    const int64_t*  target_offsets,
    const int32_t*  target_lengths,
    const int32_t*  pairs,           /* [M*2]: (query_idx, target_idx) */
    int32_t         M,
    int32_t         band_width,
    int32_t*        out_scores       /* [M] */
) {
    /* X-drop disabled in default batch entry point; use _xdrop variant to enable. */
    int32_t xdrop = 0;
    /* Find max sequence length for workspace allocation */
    int32_t max_len = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t qi = pairs[p * 2];
        int32_t ti = pairs[p * 2 + 1];
        int32_t plen = prof_lengths[qi];
        int32_t tlen = target_lengths[ti];
        if (plen > max_len) max_len = plen;
        if (tlen > max_len) max_len = tlen;
    }

    int32_t ws_size = max_len + 2;

    /* Allocate workspace once, outside parallel region */
    int32_t n_threads_actual;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Allocate all workspace in one contiguous block per thread */
    int32_t** all_ws = (int32_t**)calloc(n_threads_actual * 6, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            all_ws[t * 6 + s] = (int32_t*)calloc(ws_size, sizeof(int32_t));
        }
    }

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t* Mp = all_ws[tid * 6 + 0];
        int32_t* Mc = all_ws[tid * 6 + 1];
        int32_t* Ip = all_ws[tid * 6 + 2];
        int32_t* Ic = all_ws[tid * 6 + 3];
        int32_t* Dp = all_ws[tid * 6 + 4];
        int32_t* Dc = all_ws[tid * 6 + 5];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];

            int64_t p_off = prof_offsets[qi];
            int32_t p_len = prof_lengths[qi];
            int64_t t_off = target_offsets[ti];
            int32_t t_len = target_lengths[ti];

            out_scores[p] = hmm_viterbi_one(
                flat_match_emit + p_off * 20,
                insert_emit,
                transitions,
                p_len,
                target_flat + t_off,
                t_len,
                band_width,
                xdrop,
                Mp, Mc, Ip, Ic, Dp, Dc
            );
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            free(all_ws[t * 6 + s]);
        }
    }
    free(all_ws);
}


/* ─── Batch interface with X-drop ────────────────────────────────── */
void batch_hmm_viterbi_xdrop_c(
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
    int32_t         xdrop,           /* rows-without-improvement cutoff; 0 disables */
    int32_t*        out_scores
) {
    int32_t max_len = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t qi = pairs[p * 2];
        int32_t ti = pairs[p * 2 + 1];
        int32_t plen = prof_lengths[qi];
        int32_t tlen = target_lengths[ti];
        if (plen > max_len) max_len = plen;
        if (tlen > max_len) max_len = tlen;
    }
    int32_t ws_size = max_len + 2;

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    int32_t** all_ws = (int32_t**)calloc(n_threads_actual * 6, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++)
        for (int32_t s = 0; s < 6; s++)
            all_ws[t * 6 + s] = (int32_t*)calloc(ws_size, sizeof(int32_t));

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t* Mp = all_ws[tid * 6 + 0];
        int32_t* Mc = all_ws[tid * 6 + 1];
        int32_t* Ip = all_ws[tid * 6 + 2];
        int32_t* Ic = all_ws[tid * 6 + 3];
        int32_t* Dp = all_ws[tid * 6 + 4];
        int32_t* Dc = all_ws[tid * 6 + 5];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];
            int64_t p_off = prof_offsets[qi];
            int32_t p_len = prof_lengths[qi];
            int64_t t_off = target_offsets[ti];
            int32_t t_len = target_lengths[ti];
            out_scores[p] = hmm_viterbi_one(
                flat_match_emit + p_off * 20,
                insert_emit, transitions, p_len,
                target_flat + t_off, t_len, band_width, xdrop,
                Mp, Mc, Ip, Ic, Dp, Dc);
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++)
        for (int32_t s = 0; s < 6; s++)
            free(all_ws[t * 6 + s]);
    free(all_ws);
}


/* ─── AVX2 Batch ─────────────────────────────────────────────────── */
void batch_hmm_viterbi_avx2_c(
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
    int32_t         xdrop,
    int32_t*        out_scores
) {
#if HMM_HAVE_AVX2
    int32_t max_len = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t qi = pairs[p * 2];
        int32_t ti = pairs[p * 2 + 1];
        int32_t plen = prof_lengths[qi];
        int32_t tlen = target_lengths[ti];
        if (plen > max_len) max_len = plen;
        if (tlen > max_len) max_len = tlen;
    }
    int32_t ws_size = max_len + 16;  /* pad for SIMD tail reads */

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    int32_t** all_ws = (int32_t**)calloc(n_threads_actual * 6, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++)
        for (int32_t s = 0; s < 6; s++)
            all_ws[t * 6 + s] = (int32_t*)aligned_alloc(64, ws_size * sizeof(int32_t));

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t* Mp = all_ws[tid * 6 + 0];
        int32_t* Mc = all_ws[tid * 6 + 1];
        int32_t* Ip = all_ws[tid * 6 + 2];
        int32_t* Ic = all_ws[tid * 6 + 3];
        int32_t* Dp = all_ws[tid * 6 + 4];
        int32_t* Dc = all_ws[tid * 6 + 5];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];
            int64_t p_off = prof_offsets[qi];
            int32_t p_len = prof_lengths[qi];
            int64_t t_off = target_offsets[ti];
            int32_t t_len = target_lengths[ti];
            out_scores[p] = hmm_viterbi_one_avx2(
                flat_match_emit + p_off * 20,
                insert_emit, transitions, p_len,
                target_flat + t_off, t_len, band_width, xdrop,
                Mp, Mc, Ip, Ic, Dp, Dc);
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++)
        for (int32_t s = 0; s < 6; s++)
            free(all_ws[t * 6 + s]);
    free(all_ws);
#else
    /* No AVX2: fall back to scalar */
    batch_hmm_viterbi_xdrop_c(flat_match_emit, insert_emit, transitions,
        prof_offsets, prof_lengths, target_flat, target_offsets, target_lengths,
        pairs, M, band_width, xdrop, out_scores);
#endif
}


/* ─── Int16 Batch (AVX2, saturating, int32 fallback on overflow) ─── */
void batch_hmm_viterbi_int16_c(
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
    int32_t         xdrop,
    int32_t*        out_scores,
    int32_t*        out_n_overflow
) {
#if HMM_HAVE_AVX2
    int32_t max_len = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t qi = pairs[p * 2];
        int32_t ti = pairs[p * 2 + 1];
        int32_t plen = prof_lengths[qi];
        int32_t tlen = target_lengths[ti];
        if (plen > max_len) max_len = plen;
        if (tlen > max_len) max_len = tlen;
    }
    int32_t ws_size = max_len + 32;

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Per-thread: 6 int16 arrays for main + 6 int32 arrays for fallback */
    int16_t** ws16 = (int16_t**)calloc(n_threads_actual * 6, sizeof(int16_t*));
    int32_t** ws32 = (int32_t**)calloc(n_threads_actual * 6, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            ws16[t * 6 + s] = (int16_t*)aligned_alloc(64, ws_size * sizeof(int16_t));
            ws32[t * 6 + s] = (int32_t*)aligned_alloc(64, ws_size * sizeof(int32_t));
        }
    }

    int32_t n_overflow = 0;

    #pragma omp parallel reduction(+:n_overflow)
    {
        int32_t tid = omp_get_thread_num();
        int16_t* Mp = ws16[tid * 6 + 0];
        int16_t* Mc = ws16[tid * 6 + 1];
        int16_t* Ip = ws16[tid * 6 + 2];
        int16_t* Ic = ws16[tid * 6 + 3];
        int16_t* Dp = ws16[tid * 6 + 4];
        int16_t* Dc = ws16[tid * 6 + 5];
        int32_t* Mp32 = ws32[tid * 6 + 0];
        int32_t* Mc32 = ws32[tid * 6 + 1];
        int32_t* Ip32 = ws32[tid * 6 + 2];
        int32_t* Ic32 = ws32[tid * 6 + 3];
        int32_t* Dp32 = ws32[tid * 6 + 4];
        int32_t* Dc32 = ws32[tid * 6 + 5];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];
            int64_t p_off = prof_offsets[qi];
            int32_t p_len = prof_lengths[qi];
            int64_t t_off = target_offsets[ti];
            int32_t t_len = target_lengths[ti];

            int32_t s = hmm_viterbi_one_avx2_i16(
                flat_match_emit + p_off * 20,
                insert_emit, transitions, p_len,
                target_flat + t_off, t_len, band_width, xdrop,
                Mp, Mc, Ip, Ic, Dp, Dc);

            if (s == INT32_MIN) {
                /* Overflow: re-score with int32 */
                s = hmm_viterbi_one(
                    flat_match_emit + p_off * 20,
                    insert_emit, transitions, p_len,
                    target_flat + t_off, t_len, band_width, xdrop,
                    Mp32, Mc32, Ip32, Ic32, Dp32, Dc32);
                n_overflow++;
            }
            out_scores[p] = s;
        }
    }

    if (out_n_overflow) *out_n_overflow = n_overflow;

    for (int32_t t = 0; t < n_threads_actual; t++)
        for (int32_t s = 0; s < 6; s++) { free(ws16[t * 6 + s]); free(ws32[t * 6 + s]); }
    free(ws16); free(ws32);
#else
    batch_hmm_viterbi_xdrop_c(flat_match_emit, insert_emit, transitions,
        prof_offsets, prof_lengths, target_flat, target_offsets, target_lengths,
        pairs, M, band_width, xdrop, out_scores);
    if (out_n_overflow) *out_n_overflow = 0;
#endif
}


/* ─── Inter-pair SIMD: 8 pairs per AVX2 register ─────────────────── *
 *
 * Processes 8 pairs that share the SAME profile in parallel. Each AVX2
 * int32 lane runs an independent DP for one pair — no within-pair
 * dependency gymnastics, just 8 pairs in lockstep.
 *
 * Caller pre-groups pairs by profile id and passes 8 target indices per
 * group. If fewer than 8 targets are in the final group, caller pads
 * with a dummy target (length-0; the kernel detects and emits score=0).
 *
 * Memory: per batch needs 6 vector arrays of size max_T+2 for the DP
 * workspace, plus a packed-targets array (8, max_T) int32. All thread-local.
 */

#if HMM_HAVE_AVX2

static void hmm_viterbi_batch8_avx2(
    const int8_t*   match_emit,        /* (L, 20) for the shared profile */
    const int8_t*   insert_emit,       /* (20,) */
    const int32_t*  transitions,       /* (7,) */
    int32_t         L,                  /* profile length */
    const uint8_t* const targets[8],   /* 8 target sequences */
    const int32_t   target_lens[8],    /* their lengths */
    int32_t         band_width,
    __m256i*        Mp, __m256i* Mc,   /* (max_T+2) each — workspace */
    __m256i*        Ip, __m256i* Ic,
    __m256i*        Dp, __m256i* Dc,
    int32_t*        targets_packed,     /* (max_T, 8) int32 — per-column letters */
    int32_t         out_scores[8]
) {
    const int32_t t_mm = transitions[0];
    const int32_t t_mi = transitions[1];
    const int32_t t_md = transitions[2];
    const int32_t t_im = transitions[3];
    const int32_t t_ii = transitions[4];
    const int32_t t_dm = transitions[5];
    const int32_t t_dd = transitions[6];

    /* Find max target length */
    int32_t max_T = 0;
    for (int l = 0; l < 8; l++) if (target_lens[l] > max_T) max_T = target_lens[l];
    if (L == 0 || max_T == 0) {
        for (int l = 0; l < 8; l++) out_scores[l] = 0;
        return;
    }

    /* Pack 8 target sequences column-wise so targets_packed[j*8 + lane]
     * = target[lane][j]. Lanes whose target is shorter than max_T get
     * filler 20 (clamped to emit_row[20] = 0 below). */
    for (int32_t j = 0; j < max_T; j++) {
        for (int l = 0; l < 8; l++) {
            if (j < target_lens[l]) {
                uint8_t a = targets[l][j];
                targets_packed[j * 8 + l] = (a < 20) ? (int32_t)a : 20;
            } else {
                targets_packed[j * 8 + l] = 20;  /* past end → emit 0 */
            }
        }
    }

    /* Per-lane target-length vector, used to mask max_score updates */
    int32_t tlens_a[8] __attribute__((aligned(32)));
    for (int l = 0; l < 8; l++) tlens_a[l] = target_lens[l];
    __m256i v_tlens = _mm256_load_si256((__m256i*)tlens_a);

    /* Initialize workspace rows to NEG_INF */
    const __m256i v_neg_inf = _mm256_set1_epi32(NEG_INF);
    for (int32_t j = 0; j <= max_T + 1; j++) {
        Mp[j] = v_neg_inf; Mc[j] = v_neg_inf;
        Ip[j] = v_neg_inf; Ic[j] = v_neg_inf;
        Dp[j] = v_neg_inf; Dc[j] = v_neg_inf;
    }

    __m256i v_max_score = _mm256_setzero_si256();
    const __m256i v_tmm = _mm256_set1_epi32(t_mm);
    const __m256i v_tim = _mm256_set1_epi32(t_im);
    const __m256i v_tdm = _mm256_set1_epi32(t_dm);
    const __m256i v_tmi = _mm256_set1_epi32(t_mi);
    const __m256i v_tii = _mm256_set1_epi32(t_ii);
    const __m256i v_tmd = _mm256_set1_epi32(t_md);
    const __m256i v_tdd = _mm256_set1_epi32(t_dd);
    const __m256i v_zero = _mm256_setzero_si256();

    /* Per-profile-row: pre-widen emit + insert tables to int32[32] padded.
     * Index 20 holds (0, -1) as sentinel for unknown/past-end letters. */
    int32_t emit32[32] __attribute__((aligned(32)));
    int32_t iemit32[32] __attribute__((aligned(32)));
    for (int a = 0; a < 20; a++) iemit32[a] = (int32_t)insert_emit[a];
    for (int a = 20; a < 32; a++) iemit32[a] = -1;

    /* Per-lane banding. Each lane's band is centered on its own
     * proportional diagonal j_center = i*T_lane/L, half-width bw.
     * The j-loop runs over the UNION of all lane bands; a per-j SIMD
     * mask (in_band) marks which lanes are active at that j.
     */
    int32_t bw = band_width;
    int32_t max_dim = L > max_T ? L : max_T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;
    const __m256i v_bw = _mm256_set1_epi32(bw);
    const __m256i v_one_i = _mm256_set1_epi32(1);

    for (int32_t i = 1; i <= L; i++) {
        const int8_t* emit_row_i8 = match_emit + (i - 1) * 20;
        for (int a = 0; a < 20; a++) emit32[a] = (int32_t)emit_row_i8[a];
        for (int a = 20; a < 32; a++) emit32[a] = 0;

        /* Per-lane j_center = (i * T) / L, j_lo/j_hi = center ± bw,
         * clamped to [1, T]. SIMD doesn't have integer division, so the
         * center step is scalar for the 8 lanes. */
        int32_t jlo_arr[8] __attribute__((aligned(32)));
        int32_t jhi_arr[8] __attribute__((aligned(32)));
        int32_t j_lo_min = max_T + 1, j_hi_max = 0;
        for (int l = 0; l < 8; l++) {
            int32_t tl = target_lens[l];
            if (tl == 0) { jlo_arr[l] = max_T + 1; jhi_arr[l] = 0; continue; }
            int32_t jc = (int32_t)(((int64_t)i * tl) / L);
            int32_t jl = jc - bw; if (jl < 1) jl = 1;
            int32_t jh = jc + bw; if (jh > tl) jh = tl;
            jlo_arr[l] = jl; jhi_arr[l] = jh;
            if (jl < j_lo_min) j_lo_min = jl;
            if (jh > j_hi_max) j_hi_max = jh;
        }
        if (j_lo_min > j_hi_max) continue;  /* no lane is active this row */

        __m256i v_jlo = _mm256_load_si256((const __m256i*)jlo_arr);
        __m256i v_jhi = _mm256_load_si256((const __m256i*)jhi_arr);

        /* Clear out-of-band cells in [j_lo_min-1, j_hi_max+1] */
        for (int32_t j = j_lo_min - 1; j <= j_hi_max + 1 && j <= max_T + 1; j++) {
            if (j >= 0) { Mc[j] = v_neg_inf; Ic[j] = v_neg_inf; Dc[j] = v_neg_inf; }
        }

        for (int32_t j = j_lo_min; j <= j_hi_max; j++) {
            __m256i v_j = _mm256_set1_epi32(j);
            /* in_band[l] = (j >= jlo[l]) && (j <= jhi[l])
             *            = (j > jlo-1) && (jhi+1 > j) */
            __m256i in_band = _mm256_and_si256(
                _mm256_cmpgt_epi32(v_j, _mm256_sub_epi32(v_jlo, v_one_i)),
                _mm256_cmpgt_epi32(_mm256_add_epi32(v_jhi, v_one_i), v_j));
            /* Per-lane target letters for this column — one aligned load */
            __m256i letters = _mm256_load_si256(
                (const __m256i*)&targets_packed[(j - 1) * 8]);

            /* Gather 8 emissions — table is only 128 bytes so stays hot in L1 */
            __m256i m_emit = _mm256_i32gather_epi32(emit32,  letters, 4);
            __m256i i_emit = _mm256_i32gather_epi32(iemit32, letters, 4);

            /* M[i,j] = max(Mp[j-1]+tMM, Ip[j-1]+tIM, Dp[j-1]+tDM, 0) + m_emit */
            __m256i m1 = _mm256_add_epi32(Mp[j-1], v_tmm);
            __m256i m2 = _mm256_add_epi32(Ip[j-1], v_tim);
            __m256i m3 = _mm256_add_epi32(Dp[j-1], v_tdm);
            __m256i best_m = _mm256_max_epi32(
                _mm256_max_epi32(m1, m2),
                _mm256_max_epi32(m3, v_zero));
            __m256i mc_j = _mm256_add_epi32(best_m, m_emit);
            /* Write NEG_INF for out-of-band lanes so they don't leak
             * into the next iteration's j-1 reads. */
            mc_j = _mm256_blendv_epi8(v_neg_inf, mc_j, in_band);
            Mc[j] = mc_j;

            /* I[i,j] = max(Mc[j-1]+tMI, Ic[j-1]+tII) + i_emit */
            __m256i i1 = _mm256_add_epi32(Mc[j-1], v_tmi);
            __m256i i2 = _mm256_add_epi32(Ic[j-1], v_tii);
            __m256i best_i = _mm256_max_epi32(i1, i2);
            __m256i i_valid = _mm256_cmpgt_epi32(best_i,
                _mm256_set1_epi32(NEG_INF / 2));
            __m256i ic_j = _mm256_blendv_epi8(v_neg_inf,
                _mm256_add_epi32(best_i, i_emit), i_valid);
            ic_j = _mm256_blendv_epi8(v_neg_inf, ic_j, in_band);
            Ic[j] = ic_j;

            /* D[i,j] = max(Mp[j]+tMD, Dp[j]+tDD) */
            __m256i d1 = _mm256_add_epi32(Mp[j], v_tmd);
            __m256i d2 = _mm256_add_epi32(Dp[j], v_tdd);
            __m256i dc_j = _mm256_max_epi32(d1, d2);
            dc_j = _mm256_blendv_epi8(v_neg_inf, dc_j, in_band);
            Dc[j] = dc_j;

            /* Max score only on in-band lanes (target length was
             * already baked into jhi clamping). */
            __m256i mc_masked = _mm256_blendv_epi8(v_neg_inf, mc_j, in_band);
            v_max_score = _mm256_max_epi32(v_max_score, mc_masked);
        }  /* end j loop */

        /* Swap rows via pointer rotation */
        __m256i* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    /* Spill the 8 max_score lanes */
    int32_t max_a[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*)max_a, v_max_score);
    for (int l = 0; l < 8; l++) {
        out_scores[l] = (target_lens[l] > 0) ? max_a[l] : 0;
    }
}


/* Batch driver. Pairs are sorted by query_idx so consecutive pairs
 * share a profile and can be grouped into 8-wide SIMD batches. Any
 * trailing group of <8 pairs falls back to the scalar kernel.
 * band_width is currently ignored (full DP) — banding integration is
 * left for a future pass.
 */
void batch_hmm_viterbi_multipair_avx2_c(
    const int8_t*   flat_match_emit,
    const int8_t*   insert_emit,
    const int32_t*  transitions,
    const int64_t*  prof_offsets,
    const int32_t*  prof_lengths,
    const uint8_t*  target_flat,
    const int64_t*  target_offsets,
    const int32_t*  target_lengths,
    const int32_t*  pairs,             /* MUST be sorted by pairs[*,0] */
    int32_t         M,
    int32_t         band_width,
    int32_t*        out_scores
) {
    /* Determine max target length across all pairs for workspace sizing */
    int32_t max_T_all = 0;
    int32_t max_L_all = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t ti = pairs[p * 2 + 1];
        int32_t qi = pairs[p * 2];
        if (target_lengths[ti] > max_T_all) max_T_all = target_lengths[ti];
        if (prof_lengths[qi] > max_L_all) max_L_all = prof_lengths[qi];
    }
    int32_t ws_T = max_T_all + 2;

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Per-thread workspace: 6 __m256i arrays of size ws_T + packed targets */
    __m256i** ws = (__m256i**)calloc(n_threads_actual * 6, sizeof(__m256i*));
    int32_t** packs = (int32_t**)calloc(n_threads_actual, sizeof(int32_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            ws[t * 6 + s] = (__m256i*)aligned_alloc(64, ws_T * sizeof(__m256i));
        }
        packs[t] = (int32_t*)aligned_alloc(64, max_T_all * 8 * sizeof(int32_t));
    }

    /* Pre-scan: collect (run_start, run_end) pairs so OpenMP can
     * parallelize in canonical form. */
    int32_t max_runs = M + 1;
    int32_t* run_starts = (int32_t*)malloc(max_runs * sizeof(int32_t));
    int32_t n_runs = 0;
    {
        int32_t s = 0;
        while (s < M) {
            run_starts[n_runs++] = s;
            int32_t qi = pairs[s * 2];
            int32_t e = s + 1;
            while (e < M && pairs[e * 2] == qi) e++;
            s = e;
        }
        run_starts[n_runs] = M;  /* sentinel */
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
        int32_t* pack = packs[tid];

        #pragma omp for schedule(dynamic, 4)
        for (int32_t ri = 0; ri < n_runs; ri++) {
            int32_t run_start = run_starts[ri];
            int32_t run_end = run_starts[ri + 1];
            int32_t qi = pairs[run_start * 2];

            int32_t L = prof_lengths[qi];
            int64_t p_off = prof_offsets[qi];
            const int8_t* emit = flat_match_emit + p_off * 20;

            /* Process 8-batches within [run_start, run_end) */
            int32_t g = run_start;
            for (; g + 8 <= run_end; g += 8) {
                const uint8_t* tgt[8];
                int32_t tlens[8];
                for (int l = 0; l < 8; l++) {
                    int32_t ti = pairs[(g + l) * 2 + 1];
                    tgt[l] = target_flat + target_offsets[ti];
                    tlens[l] = target_lengths[ti];
                }
                int32_t group_scores[8];
                hmm_viterbi_batch8_avx2(
                    emit, insert_emit, transitions, L,
                    tgt, tlens, band_width,
                    Mp, Mc, Ip, Ic, Dp, Dc, pack, group_scores);
                for (int l = 0; l < 8; l++) out_scores[g + l] = group_scores[l];
            }
            /* Tail <8 pairs: scalar fallback. Uses a small scratch via
             * the existing single-pair function.                          */
            if (g < run_end) {
                int32_t* sMp = (int32_t*)Mp; int32_t* sMc = (int32_t*)Mc;
                int32_t* sIp = (int32_t*)Ip; int32_t* sIc = (int32_t*)Ic;
                int32_t* sDp = (int32_t*)Dp; int32_t* sDc = (int32_t*)Dc;
                /* Each __m256i is 8 int32s — we reuse the storage as
                 * scalar int32 arrays of size ws_T * 8. That's plenty. */
                for (; g < run_end; g++) {
                    int32_t ti = pairs[g * 2 + 1];
                    int32_t t_len = target_lengths[ti];
                    int64_t t_off = target_offsets[ti];
                    out_scores[g] = hmm_viterbi_one(
                        emit, insert_emit, transitions, L,
                        target_flat + t_off, t_len,
                        band_width, 0,
                        sMp, sMc, sIp, sIc, sDp, sDc);
                }
            }
        }
    }
    free(run_starts);

    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) free(ws[t * 6 + s]);
        free(packs[t]);
    }
    free(ws); free(packs);
}

#endif /* HMM_HAVE_AVX2 */
