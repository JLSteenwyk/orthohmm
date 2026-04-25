/*
 * MMseqs2-style batched sorted-merge prefilter (Phase A only).
 *
 * Instead of looping queries (outer) and for each query looking up its
 * k-mers in the target CSR, we:
 *   1. Extract (kmer, qid) for ALL queries in a chunk → packed int64s
 *   2. Sort by kmer
 *   3. Walk the sorted list, matching against target's CSR posting lists
 *      (already sorted by kmer). For each kmer present in both, emit
 *      count increments for the Q × T cross-product.
 *
 * Query chunking keeps the counts matrix (chunk_size × num_db) fitting
 * in L3 cache. Parallelism across chunks via OpenMP.
 *
 * This is an EXPERIMENT: the current CSR prefilter's per-query counts
 * array fits in L2 and already has good cache behavior. The merge
 * approach has better kmer-side sequential access but a bigger counts
 * footprint — which wins depends on dataset size. We measure.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -shared -fPIC -o kmer_merge.so kmer_merge.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void kmer_merge_set_num_threads(int32_t n) { omp_set_num_threads(n); }

static int cmp_i64(const void* a, const void* b) {
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    return (va > vb) - (va < vb);
}


/* Extract (kmer, qid_in_chunk) packed int64s for a chunk of queries,
 * filter out high-frequency k-mers. qid_in_chunk < chunk_size.
 * Encoding: (kmer << 22) | qid_in_chunk   — supports up to ~4M qid per chunk.
 * Kmer can be up to alpha^k. For alpha=10, k=6: 1M → fits in 22-bit space
 * swap. We use a generous packing: (kmer << 20) | qid so up to 1M qids. */
static int64_t extract_chunk_kmers(
    const uint8_t* q_flat,
    const int64_t* q_offsets,
    const int32_t* q_lengths,
    int32_t chunk_start,
    int32_t chunk_size,
    int32_t k,
    int32_t alpha_size,
    const int32_t* kmer_freqs,
    int32_t freq_thresh,
    int64_t* out /* caller buffer */
) {
    int64_t idx = 0;
    for (int32_t qc = 0; qc < chunk_size; qc++) {
        int32_t qi = chunk_start + qc;
        int32_t q_len = q_lengths[qi];
        int32_t num_kmers = q_len - k + 1;
        if (num_kmers <= 0) continue;
        const uint8_t* q = q_flat + q_offsets[qi];
        for (int32_t qp = 0; qp < num_kmers; qp++) {
            int64_t kmer = 0;
            int valid = 1;
            for (int32_t j = 0; j < k; j++) {
                uint8_t r = q[qp + j];
                if (r >= (uint8_t)alpha_size) { valid = 0; break; }
                kmer = kmer * alpha_size + r;
            }
            if (!valid) continue;
            if (kmer_freqs[kmer] > freq_thresh) continue;
            /* Pack (kmer << 20) | qc; qc < chunk_size <= 2^20 */
            out[idx++] = (kmer << 20) | (int64_t)qc;
        }
    }
    return idx;
}


/* Batched sorted-merge Phase A. Produces per-query sorted-by-count
 * candidate list of size `max_cands`.
 *
 * Interface mirrors score_query_full's output:
 *   out_targets [nq * max_cands] int32
 *   out_counts  [nq] int32 — number of candidates per query
 */
void batch_phase_a_merge_c(
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
    int32_t         max_cands,
    int32_t*        out_targets,
    int32_t*        out_counts
) {
    /* Chunk size tuned so (chunk × num_db × 2B) fits in L3.
     * For num_db=25000, chunk=256 → 12.5 MB (fits in most L3s). */
    int32_t chunk_size = 256;
    if (chunk_size > nq) chunk_size = nq;
    int32_t n_chunks = (nq + chunk_size - 1) / chunk_size;

    int32_t n_threads_actual = 1;
    #pragma omp parallel
    {
        #pragma omp single
        n_threads_actual = omp_get_num_threads();
    }

    /* Per-thread scratch: kmer buffer + counts matrix */
    int64_t max_kmers_per_chunk = 0;
    for (int32_t c = 0; c < n_chunks; c++) {
        int32_t start = c * chunk_size;
        int32_t sz = (start + chunk_size <= nq) ? chunk_size : (nq - start);
        int64_t sum = 0;
        for (int32_t qc = 0; qc < sz; qc++) {
            int32_t ql = q_lengths[start + qc];
            int64_t nk = ql - k + 1;
            if (nk > 0) sum += nk;
        }
        if (sum > max_kmers_per_chunk) max_kmers_per_chunk = sum;
    }

    int64_t** kmer_buffers = (int64_t**)calloc(n_threads_actual, sizeof(int64_t*));
    int16_t** counts_bufs  = (int16_t**)calloc(n_threads_actual, sizeof(int16_t*));
    for (int32_t t = 0; t < n_threads_actual; t++) {
        kmer_buffers[t] = (int64_t*)malloc(max_kmers_per_chunk * sizeof(int64_t));
        counts_bufs[t]  = (int16_t*)calloc((size_t)chunk_size * num_db, sizeof(int16_t));
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (int32_t c = 0; c < n_chunks; c++) {
        int32_t tid = omp_get_thread_num();
        int32_t chunk_start = c * chunk_size;
        int32_t chunk_sz = (chunk_start + chunk_size <= nq) ? chunk_size : (nq - chunk_start);

        int64_t* kbuf = kmer_buffers[tid];
        int16_t* counts = counts_bufs[tid];

        /* Zero out counts for this chunk */
        memset(counts, 0, (size_t)chunk_sz * num_db * sizeof(int16_t));

        /* Step 1: extract */
        int64_t n_entries = extract_chunk_kmers(
            q_flat, q_offsets, q_lengths, chunk_start, chunk_sz,
            k, alpha_size, kmer_freqs, freq_thresh, kbuf);

        /* Step 2: sort by kmer (qsort on packed int64 — kmer is high bits) */
        qsort(kbuf, n_entries, sizeof(int64_t), cmp_i64);

        /* Step 3: merge-walk. For each kmer with matches, iterate queries
         * that have it and targets that have it, emit counts. */
        int64_t i = 0;
        while (i < n_entries) {
            int64_t kmer = kbuf[i] >> 20;
            int64_t run_end = i + 1;
            while (run_end < n_entries && (kbuf[run_end] >> 20) == kmer) run_end++;

            int64_t t_start = kmer_offsets[kmer];
            int64_t t_end   = kmer_offsets[kmer + 1];
            int64_t n_t = t_end - t_start;
            int64_t n_q = run_end - i;

            if (n_t > 0 && n_q > 0) {
                for (int64_t qi_idx = i; qi_idx < run_end; qi_idx++) {
                    int32_t qc = (int32_t)(kbuf[qi_idx] & 0xFFFFF);
                    int16_t* qrow = counts + (int64_t)qc * num_db;
                    for (int64_t h = t_start; h < t_end; h++) {
                        int32_t t_id = (int32_t)(kmer_entries[h] >> 32);
                        if (qrow[t_id] < 32767) qrow[t_id]++;
                    }
                }
            }
            i = run_end;
        }

        /* Step 4: extract candidates per query */
        for (int32_t qc = 0; qc < chunk_sz; qc++) {
            int32_t qi = chunk_start + qc;
            int32_t* row = out_targets + (int64_t)qi * max_cands;
            int16_t* qrow = counts + (int64_t)qc * num_db;
            int32_t nc = 0;
            for (int32_t t_id = 0; t_id < num_db && nc < max_cands; t_id++) {
                if (qrow[t_id] >= min_total_hits) row[nc++] = t_id;
            }
            out_counts[qi] = nc;
        }
    }

    for (int32_t t = 0; t < n_threads_actual; t++) {
        free(kmer_buffers[t]);
        free(counts_bufs[t]);
    }
    free(kmer_buffers);
    free(counts_bufs);
}
