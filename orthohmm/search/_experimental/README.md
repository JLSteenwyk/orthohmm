# Experimental kernels

This directory holds implementations that were built and measured during
development but are **not** in the production path. They are kept here
because they document the engineering iterations and negative results
discussed in the accompanying paper. The production runtime never loads
them.

Do not edit these files expecting them to affect installed OrthoHMM
behavior. `setup.py` deliberately skips this directory.

## Kernels

| File | Purpose | Outcome |
|------|---------|---------|
| `viterbi_cuda_numba.py` | Naive 1-thread-per-pair CUDA via `numba.cuda`, workspace in global memory | ~10× slower than CPU due to 500-cycle global-memory latency per DP cell. Replaced by `csrc/hmm_viterbi_cuda.cu` (warp-collaborative, shared memory). |
| `csrc/hmm_viterbi_striped.c` | Farrar-striped SIMD Viterbi (AVX2 int32), with lazy-D wraparound correction | Unbanded: 100% parity, 2.8–3.6× over full-DP scalar. Banded: 0.8–0.9× — per-lane band masks cost more than they save. The shipped "multipair AVX2" kernel (2.7× over banded scalar) wins on realistic workloads. |
| `csrc/kmer_merge.c` | MMseqs2-style batched sorted-merge Phase A prefilter: cross-query k-mer batching + merge-join against target posting lists | 100% parity but 7.7× slower than the shipped CSR prefilter. Merge pattern produces a 10 MB counts matrix per chunk (doesn't fit in L2); CSR's per-query 40 KB counts array is already cache-optimal at our problem sizes. |
| `csrc/hmm_forward.c` | Log-sum-exp Forward algorithm (3-state HMM, same interface as Viterbi) | 100% correct, 25× slower than Viterbi, no F-score improvement in the Forward-vs-Viterbi sweep. Kept for methods completeness. |

See `docs/DEVELOPMENT_HISTORY.md` in the repo root for the full story of
each iteration, with numbers.
