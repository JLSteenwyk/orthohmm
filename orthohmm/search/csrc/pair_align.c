/*
 * Banded Needleman-Wunsch pairwise aligner with traceback.
 *
 * Used as an in-tree MAFFT replacement for center-star MSA construction
 * during Pass 3 profile building. Only needs to be fast enough for
 * thousands of short-ish pairs; quality matters more than raw speed since
 * the results feed PSSM construction.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -shared -fPIC -o pair_align.so pair_align.c
 *
 * Interface:
 *   batch_pair_align_c(
 *       seqs_flat, seq_offsets, seq_lengths,      # all sequences
 *       sub_matrix,                                # 20x20 int8
 *       gap_open, gap_extend,                     # int32
 *       pairs, M,                                  # (M,2) query_idx, target_idx
 *       band_width,                                # int32 (0 = full DP)
 *       out_align_a, out_align_b,                 # (M * max_align_len) uint8
 *                                                  #  values 0-19 AA, 20 gap
 *       out_align_lens,                           # (M,) int32 alignment length
 *       max_align_len                             # capacity
 *   )
 *
 * Returns alignment pair for each input pair: out_align_a[m, :L] and
 * out_align_b[m, :L] are parallel strings of length L = out_align_lens[m].
 * Gap is 20.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void pair_align_set_num_threads(int32_t n) {
    omp_set_num_threads(n);
}

#define NEG_INF (-1000000)

/* Traceback directions */
#define TB_STOP 0
#define TB_DIAG 1   /* match/mismatch */
#define TB_UP   2   /* gap in seq b (consume seq a) */
#define TB_LEFT 3   /* gap in seq a (consume seq b) */

static int32_t pair_align_one(
    const uint8_t*  seq_a, int32_t la,
    const uint8_t*  seq_b, int32_t lb,
    const int8_t*   sub_matrix,   /* 20x20 row-major */
    int32_t         gap_open,
    int32_t         gap_extend,
    int32_t         band_width,   /* 0 or <=0 means full DP */
    int32_t*        score_ws,     /* (la+1)*(lb+1) int32 */
    uint8_t*        tb_ws,        /* (la+1)*(lb+1) uint8 */
    uint8_t*        out_a,        /* caller buffer, at least la+lb long */
    uint8_t*        out_b,        /* same */
    int32_t         max_align_len,
    int32_t*        out_len
) {
    int32_t cols = lb + 1;

    /* Initialize row 0 and column 0. We use NW-style scoring (global);
     * for orthologs this is appropriate because we expect them to be
     * alignable end-to-end. */
    int32_t* H = score_ws;
    uint8_t* TB = tb_ws;

    H[0] = 0;
    TB[0] = TB_STOP;
    for (int32_t j = 1; j <= lb; j++) {
        H[j] = gap_open + (j - 1) * gap_extend;
        TB[j] = TB_LEFT;
    }
    for (int32_t i = 1; i <= la; i++) {
        int32_t row = i * cols;
        H[row] = gap_open + (i - 1) * gap_extend;
        TB[row] = TB_UP;
    }

    /* Banded fill */
    for (int32_t i = 1; i <= la; i++) {
        int32_t row = i * cols;
        int32_t prev_row = (i - 1) * cols;

        int32_t j_lo = 1, j_hi = lb;
        if (band_width > 0) {
            /* Proportional diagonal band */
            int32_t j_center = (int32_t)(((int64_t)i * lb) / la);
            j_lo = j_center - band_width;
            j_hi = j_center + band_width;
            if (j_lo < 1) j_lo = 1;
            if (j_hi > lb) j_hi = lb;
            /* Zero out out-of-band cells in this row so traceback doesn't
             * wander off. NEG_INF prevents them from being chosen. */
            for (int32_t j = 1; j < j_lo; j++) { H[row + j] = NEG_INF; TB[row + j] = TB_STOP; }
            for (int32_t j = j_hi + 1; j <= lb; j++) { H[row + j] = NEG_INF; TB[row + j] = TB_STOP; }
        }

        uint8_t a_aa = seq_a[i - 1];
        const int8_t* sub_row = (a_aa < 20) ? (sub_matrix + a_aa * 20) : NULL;

        for (int32_t j = j_lo; j <= j_hi; j++) {
            uint8_t b_aa = seq_b[j - 1];
            int32_t match_score;
            if (sub_row && b_aa < 20) {
                match_score = (int32_t)sub_row[b_aa];
            } else {
                match_score = 0;  /* unknown AA: neutral */
            }

            int32_t diag = H[prev_row + j - 1] + match_score;

            /* Affine gap. Simplified: any gap continuation uses gap_extend,
             * a gap opening uses gap_open. We track gap state implicitly via
             * TB (less exact than Gotoh but sufficient for orthologs).       */
            int32_t up_prev = H[prev_row + j];
            int32_t up_cost = (TB[prev_row + j] == TB_UP) ? gap_extend : gap_open;
            int32_t up = up_prev + up_cost;

            int32_t left_prev = H[row + j - 1];
            int32_t left_cost = (TB[row + j - 1] == TB_LEFT) ? gap_extend : gap_open;
            int32_t left = left_prev + left_cost;

            int32_t best = diag;
            uint8_t dir = TB_DIAG;
            if (up > best) { best = up; dir = TB_UP; }
            if (left > best) { best = left; dir = TB_LEFT; }

            H[row + j] = best;
            TB[row + j] = dir;
        }
    }

    /* Traceback from (la, lb) */
    int32_t i = la, j = lb;
    int32_t pos = 0;
    /* Write from the end: we'll reverse at the end. Use the caller's
     * buffer in reverse order. */
    while (i > 0 || j > 0) {
        if (pos >= max_align_len) {
            *out_len = 0;
            return -1;  /* overflow */
        }
        uint8_t dir;
        if (i == 0) dir = TB_LEFT;
        else if (j == 0) dir = TB_UP;
        else dir = TB[i * cols + j];

        if (dir == TB_DIAG) {
            out_a[pos] = seq_a[i - 1];
            out_b[pos] = seq_b[j - 1];
            i--; j--;
        } else if (dir == TB_UP) {
            out_a[pos] = seq_a[i - 1];
            out_b[pos] = 20;  /* gap */
            i--;
        } else if (dir == TB_LEFT) {
            out_a[pos] = 20;
            out_b[pos] = seq_b[j - 1];
            j--;
        } else {
            break;  /* STOP — shouldn't happen in NW */
        }
        pos++;
    }

    /* Reverse in place */
    for (int32_t k = 0; k < pos / 2; k++) {
        uint8_t ta = out_a[k]; out_a[k] = out_a[pos - 1 - k]; out_a[pos - 1 - k] = ta;
        uint8_t tb = out_b[k]; out_b[k] = out_b[pos - 1 - k]; out_b[pos - 1 - k] = tb;
    }

    *out_len = pos;
    return H[la * cols + lb];
}


void batch_pair_align_c(
    const uint8_t*  seqs_flat,        /* all sequences concatenated */
    const int64_t*  seq_offsets,
    const int32_t*  seq_lengths,
    const int8_t*   sub_matrix,        /* (20,20) row-major */
    int32_t         gap_open,
    int32_t         gap_extend,
    const int32_t*  pairs,              /* (M, 2): (a_idx, b_idx) */
    int32_t         M,
    int32_t         band_width,
    uint8_t*        out_align_a,        /* (M, max_align_len) */
    uint8_t*        out_align_b,        /* (M, max_align_len) */
    int32_t*        out_align_lens,     /* (M,) */
    int32_t         max_align_len,
    int32_t*        out_scores          /* (M,) int32 */
) {
    /* Find max (la+1)*(lb+1) across pairs for workspace sizing */
    int64_t max_cells = 0;
    for (int32_t p = 0; p < M; p++) {
        int32_t ai = pairs[p * 2];
        int32_t bi = pairs[p * 2 + 1];
        int64_t la = seq_lengths[ai];
        int64_t lb = seq_lengths[bi];
        int64_t cells = (la + 1) * (lb + 1);
        if (cells > max_cells) max_cells = cells;
    }

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Per-thread workspace */
    int32_t** score_wss = (int32_t**)calloc(n_threads_actual, sizeof(int32_t*));
    uint8_t** tb_wss    = (uint8_t**)calloc(n_threads_actual, sizeof(uint8_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        score_wss[t] = (int32_t*)malloc(max_cells * sizeof(int32_t));
        tb_wss[t]    = (uint8_t*)malloc(max_cells * sizeof(uint8_t));
    }

    #pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t* score_ws = score_wss[tid];
        uint8_t* tb_ws    = tb_wss[tid];

        #pragma omp for schedule(dynamic, 32)
        for (int32_t p = 0; p < M; p++) {
            int32_t ai = pairs[p * 2];
            int32_t bi = pairs[p * 2 + 1];
            int64_t a_off = seq_offsets[ai];
            int64_t b_off = seq_offsets[bi];
            int32_t la = seq_lengths[ai];
            int32_t lb = seq_lengths[bi];

            uint8_t* out_a = out_align_a + (int64_t)p * max_align_len;
            uint8_t* out_b = out_align_b + (int64_t)p * max_align_len;
            int32_t align_len = 0;

            int32_t s = pair_align_one(
                seqs_flat + a_off, la,
                seqs_flat + b_off, lb,
                sub_matrix, gap_open, gap_extend,
                band_width, score_ws, tb_ws,
                out_a, out_b, max_align_len, &align_len
            );

            out_align_lens[p] = align_len;
            if (out_scores) out_scores[p] = s;
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++) {
        free(score_wss[t]);
        free(tb_wss[t]);
    }
    free(score_wss);
    free(tb_wss);
}
