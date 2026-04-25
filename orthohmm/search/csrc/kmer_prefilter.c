/*
 * C/OpenMP k-mer prefilter with ungapped extension re-ranking.
 *
 * Pipeline per query:
 *   Phase A: k-mer counting → identify candidate targets
 *   Phase A.5: Ungapped extension scoring → re-rank candidates by alignment potential
 *   Top-K selection by extension score → high-quality candidate set
 *
 * The key insight: k-mer count is a poor proxy for alignment quality.
 * Ungapped extension (Kadane's algorithm with substitution matrix) on
 * 3 diagonals gives a much better estimate for ~3x the cost of Phase A
 * alone, but produces dramatically better candidates for downstream
 * HMM Viterbi scoring.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -shared -fPIC -o kmer_prefilter.so kmer_prefilter.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void prefilter_set_num_threads(int32_t n) {
    omp_set_num_threads(n);
}

static int cmp_score_desc(const void* a, const void* b) {
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    return (vb > va) - (vb < va);
}

/* ─── Phase A: k-mer counting with prefetch ─────────────────────── */

static void phase_a_score(
    const uint8_t*  q_seq,
    int32_t         q_len,
    int32_t         k,
    int32_t         alpha_size,
    const int64_t*  kmer_offsets,
    const int64_t*  kmer_entries,
    const int32_t*  kmer_freqs,
    int32_t         freq_thresh,
    int32_t         num_db,
    int16_t*        target_counts
) {
    int32_t num_kmers = q_len - k + 1;
    if (num_kmers <= 0) return;

    for (int32_t qpos = 0; qpos < num_kmers; qpos++) {
        int64_t kmer_val = 0;
        int valid = 1;
        for (int32_t j = 0; j < k; j++) {
            uint8_t r = q_seq[qpos + j];
            if (r >= (uint8_t)alpha_size) { valid = 0; break; }
            kmer_val = kmer_val * alpha_size + r;
        }
        if (!valid) continue;
        if (kmer_freqs[kmer_val] > freq_thresh) continue;

        int64_t s = kmer_offsets[kmer_val];
        int64_t e = kmer_offsets[kmer_val + 1];

        /* Prefetch next posting list */
        if (qpos + 1 < num_kmers) {
            int64_t nk = 0; int nv = 1;
            for (int32_t j = 0; j < k; j++) {
                uint8_t r = q_seq[qpos + 1 + j];
                if (r >= (uint8_t)alpha_size) { nv = 0; break; }
                nk = nk * alpha_size + r;
            }
            if (nv && kmer_freqs[nk] <= freq_thresh)
                __builtin_prefetch(&kmer_entries[kmer_offsets[nk]], 0, 1);
        }

        for (int64_t h = s; h < e; h++) {
            int32_t tid = (int32_t)(kmer_entries[h] >> 32);
            if (target_counts[tid] < 32767)
                target_counts[tid]++;
        }
    }
}


/* ─── Ungapped extension (Kadane's algorithm) ───────────────────── */

static int32_t ungapped_extend_diag(
    const uint8_t*  q_seq, int32_t q_len,
    const uint8_t*  t_seq, int32_t t_len,
    int32_t         diag_offset,   /* t_start - q_start */
    const int8_t*   sub_matrix     /* 20x20 row-major */
) {
    int32_t a_start, b_start;
    if (diag_offset >= 0) {
        a_start = 0;
        b_start = diag_offset;
    } else {
        a_start = -diag_offset;
        b_start = 0;
    }

    int32_t a_remain = q_len - a_start;
    int32_t b_remain = t_len - b_start;
    int32_t diag_len = a_remain < b_remain ? a_remain : b_remain;
    if (diag_len <= 0) return 0;

    int32_t max_score = 0;
    int32_t curr_score = 0;
    for (int32_t i = 0; i < diag_len; i++) {
        uint8_t a = q_seq[a_start + i];
        uint8_t b = t_seq[b_start + i];
        int32_t s = (a < 20 && b < 20) ? (int32_t)sub_matrix[a * 20 + b] : -4;
        curr_score += s;
        if (curr_score > max_score) max_score = curr_score;
        if (curr_score < 0) curr_score = 0;
    }
    return max_score;
}

/* Score query vs target on 3 representative diagonals, return best */
static int32_t ungapped_extend_3diag(
    const uint8_t*  q_seq, int32_t q_len,
    const uint8_t*  t_seq, int32_t t_len,
    const int8_t*   sub_matrix
) {
    /* Three diagonals: start-aligned, center-aligned, end-aligned */
    int32_t diag_center = (t_len - q_len) / 2;
    int32_t diag_start = 0;
    int32_t diag_end = t_len - q_len;

    int32_t s1 = ungapped_extend_diag(q_seq, q_len, t_seq, t_len, diag_start, sub_matrix);
    int32_t s2 = ungapped_extend_diag(q_seq, q_len, t_seq, t_len, diag_center, sub_matrix);
    int32_t s3 = ungapped_extend_diag(q_seq, q_len, t_seq, t_len, diag_end, sub_matrix);

    int32_t best = s1;
    if (s2 > best) best = s2;
    if (s3 > best) best = s3;
    return best;
}


/* ─── Full pipeline: Phase A → extension re-rank → top-K ────────── */

static int32_t score_query_full(
    /* Remapped query (for k-mer lookup) */
    const uint8_t*  q_seq_remapped,
    int32_t         q_len,
    int32_t         k,
    int32_t         alpha_size,
    const int64_t*  kmer_offsets,
    const int64_t*  kmer_entries,
    const int32_t*  kmer_freqs,
    int32_t         freq_thresh,
    int32_t         num_db,
    int32_t         min_total_hits,
    int32_t         phase_a_topk,
    /* Original sequences (for extension scoring) */
    const uint8_t*  q_seq_original,
    const uint8_t*  t_flat_original,
    const int64_t*  t_offsets,
    const int32_t*  t_lengths,
    const int8_t*   sub_matrix,
    /* Output */
    int32_t         max_cands,
    int32_t*        out_ids
) {
    int32_t num_kmers = q_len - k + 1;
    if (num_kmers <= 0) return 0;

    /* Phase A: k-mer counting */
    int16_t* counts = (int16_t*)calloc(num_db, sizeof(int16_t));
    if (!counts) return 0;

    phase_a_score(q_seq_remapped, q_len, k, alpha_size,
                  kmer_offsets, kmer_entries, kmer_freqs,
                  freq_thresh, num_db, counts);

    /* Collect Phase A survivors */
    int32_t num_passing = 0;
    for (int32_t i = 0; i < num_db; i++)
        if (counts[i] >= min_total_hits) num_passing++;

    if (num_passing == 0) { free(counts); return 0; }

    /* Select top phase_a_topk by k-mer count */
    int32_t* pass_ids = (int32_t*)malloc(num_passing * sizeof(int32_t));
    int32_t p = 0;
    for (int32_t i = 0; i < num_db; i++)
        if (counts[i] >= min_total_hits)
            pass_ids[p++] = i;

    int32_t topk = phase_a_topk < num_passing ? phase_a_topk : num_passing;
    if (topk < num_passing) {
        int64_t* pk = (int64_t*)malloc(num_passing * sizeof(int64_t));
        for (int32_t i = 0; i < num_passing; i++)
            pk[i] = ((int64_t)(uint16_t)counts[pass_ids[i]] << 32)
                   | (uint32_t)(0x7FFFFFFF - pass_ids[i]);
        qsort(pk, num_passing, sizeof(int64_t), cmp_score_desc);
        for (int32_t i = 0; i < topk; i++)
            pass_ids[i] = (int32_t)(0x7FFFFFFF - (pk[i] & 0xFFFFFFFF));
        free(pk);
    }
    free(counts);

    /* Phase A.5: Ungapped extension re-ranking */
    int64_t* ext_packed = (int64_t*)malloc(topk * sizeof(int64_t));
    for (int32_t i = 0; i < topk; i++) {
        int32_t tid = pass_ids[i];
        const uint8_t* t_seq = t_flat_original + t_offsets[tid];
        int32_t t_len = t_lengths[tid];

        int32_t ext_score = ungapped_extend_3diag(
            q_seq_original, q_len, t_seq, t_len, sub_matrix
        );

        /* Pack as (score << 32 | inverted_id) for descending sort */
        ext_packed[i] = ((int64_t)(uint32_t)ext_score << 32)
                       | (uint32_t)(0x7FFFFFFF - tid);
    }
    free(pass_ids);

    /* Sort by extension score descending */
    qsort(ext_packed, topk, sizeof(int64_t), cmp_score_desc);

    /* Return top max_cands */
    int32_t nc = topk < max_cands ? topk : max_cands;
    for (int32_t i = 0; i < nc; i++)
        out_ids[i] = (int32_t)(0x7FFFFFFF - (ext_packed[i] & 0xFFFFFFFF));

    free(ext_packed);
    return nc;
}


/* ─── Batch entry point (OpenMP parallel over queries) ───────────── */

void batch_prefilter_c(
    /* Remapped queries (for k-mer lookup) */
    const uint8_t*  q_flat,
    const int64_t*  q_offsets,
    const int32_t*  q_lengths,
    int32_t         nq,
    int32_t         k,
    int32_t         alpha_size,
    const int64_t*  kmer_offsets,
    const int64_t*  kmer_entries,
    const int32_t*  kmer_freqs,
    int32_t         freq_thresh,
    int32_t         num_db,
    int32_t         min_total_hits,
    int32_t         min_diag_hits,     /* unused, kept for API compat */
    int32_t         diag_bin_width,    /* unused, kept for API compat */
    int32_t         max_cands,
    int32_t         phase_a_topk,
    int32_t*        out_targets,
    int32_t*        out_counts,
    /* NEW: original sequences + sub matrix for extension scoring */
    const uint8_t*  q_flat_original,   /* NULL = skip extension, use kmer count */
    const uint8_t*  t_flat_original,
    const int64_t*  t_offsets_original,
    const int32_t*  t_lengths_original,
    const int8_t*   sub_matrix
) {
    int use_extension = (q_flat_original != NULL && sub_matrix != NULL);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int32_t qi = 0; qi < nq; qi++) {
        int32_t q_len = q_lengths[qi];
        int32_t* row = out_targets + (int64_t)qi * max_cands;

        if (use_extension) {
            const uint8_t* q_remap = q_flat + q_offsets[qi];
            const uint8_t* q_orig  = q_flat_original + q_offsets[qi];

            out_counts[qi] = score_query_full(
                q_remap, q_len, k, alpha_size,
                kmer_offsets, kmer_entries, kmer_freqs,
                freq_thresh, num_db, min_total_hits, phase_a_topk,
                q_orig, t_flat_original, t_offsets_original,
                t_lengths_original, sub_matrix,
                max_cands, row
            );
        } else {
            /* Fallback: Phase A only, rank by kmer count */
            int32_t num_kmers = q_len - k + 1;
            if (num_kmers <= 0) { out_counts[qi] = 0; continue; }

            const uint8_t* q_seq = q_flat + q_offsets[qi];
            int16_t* counts = (int16_t*)calloc(num_db, sizeof(int16_t));
            if (!counts) { out_counts[qi] = 0; continue; }

            phase_a_score(q_seq, q_len, k, alpha_size,
                          kmer_offsets, kmer_entries, kmer_freqs,
                          freq_thresh, num_db, counts);

            int32_t nc = 0;
            for (int32_t i = 0; i < num_db && nc < max_cands; i++)
                if (counts[i] >= min_total_hits)
                    row[nc++] = i;

            out_counts[qi] = nc;
            free(counts);
        }
    }
}
