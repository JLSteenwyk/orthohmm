/*
 * C/OpenMP banded 3-state profile HMM Forward algorithm.
 *
 * Mirrors hmm_viterbi.c structure but uses log-sum-exp instead of max
 * to sum over all alignment paths (Forward) rather than taking the
 * best single path (Viterbi).
 *
 * Forward scoring is more sensitive for distant homologs where
 * multiple plausible alignments exist — it captures alignment
 * uncertainty that Viterbi discards.
 *
 * Scores are in log-space (natural log, integer-scaled *100 for int32).
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -shared -fPIC -o hmm_forward.so hmm_forward.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

void forward_set_num_threads(int32_t n) {
    omp_set_num_threads(n);
}

#define NEG_INF_F  (-1.0e9f)
#define LOG2E       1.4426950408889634f
#define LN2         0.6931471805599453f

/* Numerically stable log-sum-exp for two values (natural log) */
static inline float logsumexp2(float a, float b) {
    if (a > b) {
        return a + logf(1.0f + expf(b - a));
    } else {
        return b + logf(1.0f + expf(a - b));
    }
}

static inline float logsumexp3(float a, float b, float c) {
    return logsumexp2(a, logsumexp2(b, c));
}


/* ─── Single-pair banded 3-state Forward ─────────────────────────── */

static float hmm_forward_one(
    const int8_t*   match_emit,     /* (prof_len, 20) row-major */
    const int8_t*   insert_emit,    /* (20,) */
    const int32_t*  transitions,    /* (7,): MM MI MD IM II DM DD */
    int32_t         prof_len,
    const uint8_t*  target_seq,
    int32_t         target_len,
    int32_t         band_width,
    float*          Mp, float* Mc,
    float*          Ip, float* Ic,
    float*          Dp, float* Dc
) {
    /* Convert integer transitions (bits*2) to natural log */
    /* Our scores are in half-bits; convert to nats: nats = hbits * ln(2)/2 */
    const float t_mm = transitions[0] * LN2 * 0.5f;
    const float t_mi = transitions[1] * LN2 * 0.5f;
    const float t_md = transitions[2] * LN2 * 0.5f;
    const float t_im = transitions[3] * LN2 * 0.5f;
    const float t_ii = transitions[4] * LN2 * 0.5f;
    const float t_dm = transitions[5] * LN2 * 0.5f;
    const float t_dd = transitions[6] * LN2 * 0.5f;

    int32_t L = prof_len;
    int32_t T = target_len;
    int32_t cols = T + 1;

    if (L == 0 || T == 0) return 0.0f;

    int32_t bw = band_width;
    int32_t max_dim = L > T ? L : T;
    if (bw <= 0 || max_dim <= 50) bw = max_dim;

    for (int32_t j = 0; j < cols; j++) {
        Mp[j] = NEG_INF_F; Mc[j] = NEG_INF_F;
        Ip[j] = NEG_INF_F; Ic[j] = NEG_INF_F;
        Dp[j] = NEG_INF_F; Dc[j] = NEG_INF_F;
    }

    /* Use local-alignment Forward: sum over all alignment "start" positions.
     * At each (i,j), M can start fresh with emission score.
     * Final score: logsumexp over all M[i,j] values (marginalize start and end).
     */
    float max_score = 0.0f;
    float total_sum = NEG_INF_F;  /* log-sum of all M scores = log(sum of prob over alignments) */

    for (int32_t i = 1; i <= L; i++) {
        int32_t j_center = (int32_t)(((int64_t)i * T) / L);
        int32_t j_lo = j_center - bw;
        int32_t j_hi = j_center + bw;
        if (j_lo < 1) j_lo = 1;
        if (j_hi > T) j_hi = T;

        for (int32_t j = j_lo - 1; j <= j_hi + 1 && j < cols; j++) {
            if (j >= 0) {
                Mc[j] = NEG_INF_F;
                Ic[j] = NEG_INF_F;
                Dc[j] = NEG_INF_F;
            }
        }

        const int8_t* emit_row = match_emit + (i - 1) * 20;

        for (int32_t j = j_lo; j <= j_hi; j++) {
            uint8_t t_aa = target_seq[j - 1];
            float m_emit = (t_aa < 20) ? ((float)emit_row[t_aa] * LN2 * 0.5f) : 0.0f;
            float i_emit = (t_aa < 20) ? ((float)insert_emit[t_aa] * LN2 * 0.5f) : -LN2;

            /* M[i,j] = m_emit + logsumexp(M[i-1,j-1]+t_mm, I[i-1,j-1]+t_im, D[i-1,j-1]+t_dm, 0) */
            float m_prev = logsumexp3(Mp[j-1] + t_mm, Ip[j-1] + t_im, Dp[j-1] + t_dm);
            /* Include "local start" term: probability of starting alignment here */
            m_prev = logsumexp2(m_prev, 0.0f);
            Mc[j] = m_prev + m_emit;

            /* I[i,j] = i_emit + logsumexp(M[i,j-1]+t_mi, I[i,j-1]+t_ii) */
            float i_prev = logsumexp2(Mc[j-1] + t_mi, Ic[j-1] + t_ii);
            Ic[j] = (i_prev > NEG_INF_F / 2) ? (i_prev + i_emit) : NEG_INF_F;

            /* D[i,j] = logsumexp(M[i-1,j]+t_md, D[i-1,j]+t_dd) */
            Dc[j] = logsumexp2(Mp[j] + t_md, Dp[j] + t_dd);

            /* Accumulate total (marginalize over end positions) */
            if (Mc[j] > NEG_INF_F / 2) {
                total_sum = logsumexp2(total_sum, Mc[j]);
                if (Mc[j] > max_score) max_score = Mc[j];
            }
        }

        float* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    /* Return Forward score in half-bits (2/ln2) for comparability with Viterbi scores */
    return total_sum * 2.0f / LN2;
}


/* ─── Batch interface ───────────────────────────────────────────── */

void batch_hmm_forward_c(
    const int8_t*   flat_match_emit,  /* (sum_L, 20) row-major */
    const int8_t*   insert_emit,      /* (20,) */
    const int32_t*  transitions,      /* (7,) */
    const int64_t*  prof_offsets,
    const int32_t*  prof_lengths,
    const uint8_t*  target_flat,
    const int64_t*  target_offsets,
    const int32_t*  target_lengths,
    const int32_t*  pairs,           /* [M*2]: (query_idx, target_idx) */
    int32_t         M,
    int32_t         band_width,
    float*          out_scores       /* [M] */
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

    float** all_ws = (float**)calloc(n_threads_actual * 6, sizeof(float*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        for (int32_t s = 0; s < 6; s++) {
            all_ws[t * 6 + s] = (float*)calloc(ws_size, sizeof(float));
        }
    }

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        float* Mp = all_ws[tid * 6 + 0];
        float* Mc = all_ws[tid * 6 + 1];
        float* Ip = all_ws[tid * 6 + 2];
        float* Ic = all_ws[tid * 6 + 3];
        float* Dp = all_ws[tid * 6 + 4];
        float* Dc = all_ws[tid * 6 + 5];

        #pragma omp for schedule(dynamic, 64)
        for (int32_t p = 0; p < M; p++) {
            int32_t qi = pairs[p * 2];
            int32_t ti = pairs[p * 2 + 1];

            int64_t p_off = prof_offsets[qi];
            int32_t p_len = prof_lengths[qi];
            int64_t t_off = target_offsets[ti];
            int32_t t_len = target_lengths[ti];

            out_scores[p] = hmm_forward_one(
                flat_match_emit + p_off * 20,
                insert_emit,
                transitions,
                p_len,
                target_flat + t_off,
                t_len,
                band_width,
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
