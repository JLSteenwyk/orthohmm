# Production Performance Optimization

Date: 2026-08-28

This work improves the production OrthoHMM CLI while preserving its computed
graph and orthogroups. QfO and OrthoBench are the primary accuracy benchmarks;
Three Kingdoms remains a secondary stress test.

## Reproducible Harness

Use `benchmark_tools/benchmark_production.py`. It invokes `python -m
orthohmm`, refuses to reuse an output directory, and records:

- the exact command, Git commit and dirty state, and source checksums;
- input and output checksums;
- stage wall and CPU time;
- sampled summed Linux process-tree RSS;
- candidates, significant hits, edges, genes, and orthogroup counts.

The historical `bacterial_scaling/run_orthohmm.py` is an experimental second
implementation. It previously normalized built-in search scores a second time.
That duplicate normalization was removed locally, and the bacterial runner now
uses the production harness. New performance claims must not use the legacy
driver.

## Accepted Changes

1. Search hits remain local `int32` query/target indices with numeric score and
   E-value arrays. Gene names are stored once per species.
2. Significant hits are filtered before multiprocessing IPC.
3. `--cpu` is a total budget. The default is four workers with eight native
   threads each on a 32-CPU run.
4. The adaptive coverage probe and candidate prefilter reuse one prepared
   target k-mer index within each species-pair task.
5. Profile emission lookup uses `numpy.take(..., out=...)`; candidate query
   indices use `numpy.repeat` instead of a Python pair-building loop.
6. Reciprocal-best-hit thresholds and graph edges remain integer indexed.
   Stable packed integer keys provide deterministic maximum-weight edge
   deduplication.
7. Leiden consumes the compact graph directly and uses seed 0. ABC is still
   written as a checkpoint and for MCL compatibility.
8. Orthogroup materialization uses one gene-to-row index instead of scanning
   the complete gene table for every orthogroup.
9. Linux process-tree RSS and detailed stage metrics are opt-in through the
   benchmark harness, so normal CLI runs do not pay monitoring overhead.

## Controlled Results

The `n005` bacterial subset contains 15,932 proteins from five proteomes.
Every accepted run produced identical edge, pre-singleton cluster, and final
orthogroup SHA-256 checksums. Counts remained 225,859 candidates, 67,337
significant hits, 5,416 edges, and 12,995 orthogroups.

| production state | wall | peak process-tree RSS | edge thresholds | graph edges | orthogroup materialization |
| --- | ---: | ---: | ---: | ---: | ---: |
| instrumented baseline | 15.77 s | 2.56 GiB | name-based | 0.24 s | included in combined stage |
| compact graph, 4 workers x 8 threads | 13.11 s | 1.36-1.38 GiB | name-based | 0.022 s | 5.65 s |
| final optimized path | 6.98-7.08 s | 1.35-1.37 GiB | 0.052 s | 0.018 s | 0.44 s |

The worker/thread comparison used a fixed 32-CPU budget:

| workers x threads | search wall | measured stage total | peak RSS |
| --- | ---: | ---: | ---: |
| 16 x 2 | 8.54 s | 15.86 s | 4.45 GiB |
| 8 x 4 | 5.99 s | 13.25 s | 2.51 GiB |
| 4 x 8 | 5.83 s | 13.16 s | 1.39 GiB |
| 2 x 16 | 8.47 s | 15.62 s | 0.77 GiB |

Four workers by eight threads is the best measured speed/memory balance and
is the recommended scheduler default. Users can select a different balance
with the hidden benchmark option `--threads_per_worker`.

The full 12-proteome OrthoBench input contains 251,378 proteins. Both runs
used the production CLI with a fixed 32-CPU budget and produced the same
45,440,595 candidates, 9,743,225 significant hits, 1,273,172 edges, and
80,347 orthogroups.

| production state | wall | peak process-tree RSS | search | edge thresholds | graph edges | clustering + refinement + materialization |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| instrumented baseline | 2,145.66 s | 2.81 GiB | 956.30 s | 130.07 s | 2.35 s | 1,056.74 s |
| final optimized path | 978.63 s | 2.83 GiB | 954.17 s | 1.95 s | 2.32 s | 20.00 s |

The optimized path reduced end-to-end wall time by 54.4% (2.19x) while the
full-input peak RSS was effectively unchanged. The large-input gain comes
from compact numeric thresholding and removal of repeated full-table scans
after search, not from changing the search result.

## Accuracy Status

The fresh pre-optimization production OrthoBench run scored F=63.8,
precision=85.2, recall=51.0, with 12 exact RefOGs. It produced 45,440,595
candidates, 9,743,225 significant hits, 1,273,172 edges, and 80,347 final
orthogroups. The final optimized production run reproduced F=63.8,
precision=85.2, recall=51.0, and 12 exact RefOGs. Its three scientific
checkpoints are byte-identical to the pre-optimization baseline:

| checkpoint | SHA-256 |
| --- | --- |
| edge list | `e8038ffd18a88a69a3ba618789671aaa1a27e801964f59d6aa9c0a671423f7f9` |
| pre-singleton clusters | `c989cdea5fd927d4318d4ef235e3560feea7f9068e8d11e1d3f1deeaba50760e` |
| final orthogroups | `16444a7fc1818a0bc2921cd6a390ed10f71147a8d2ccae402acb86e9c8c064d8` |

The high-sensitivity path now productionizes the bounded k=4 search,
multi-pass graph inference, and strict profile expansion. A committed-source
replay of the historical fixed-`mc100` checkpoint reproduced F=70.4,
precision=78.9, recall=63.5, with 13 exact RefOGs. The first fresh production
run scored F=69.6, precision=79.1, recall=62.1, with 12 exact RefOGs. It used
2,203.51 seconds and 11.87 GiB peak process-tree RSS on 32 CPUs. The fresh run
therefore improved substantially over the 63.8 production baseline but did not
initially reproduce the cached result.

The fresh run produced 298,209,630 prefilter candidates, 18,486,508
significant hits, 1,698,762 RBNH edges, and 1,865,672 final graph edges. The
historical checkpoint contained 18,235,373 significant hits and 1,803,122 RBNH
edges. Pair-level diagnostics showed that the biological hit sets agreed; the
global significant-hit difference was one same-sequence self-hit per protein.
Those self-hits had been excluded from the historical RBNH calculation but
were retained while production estimated same-species reciprocal-best
thresholds. They raised the threshold and removed valid paralog edges.

Production now excludes self-hits before estimating reciprocal-best
thresholds. A fresh run from the corrected source scored F=70.4,
precision=78.9, recall=63.5, with 13 exact RefOGs. It produced 1,803,122 RBNH
edges and 1,935,903 final graph edges. The refined cluster SHA-256 is
`8ee100f4370fef6e3e6a3e0fd36112adcc3aefeb60197331115667c519ca121f`,
exactly matching the historical checkpoint. This clears the first production
milestone and improves F-score by 6.6 points over the prior default, but it
remains 2.3 points below OrthoFinder 3.1.5 at F=72.7.

The corrected run used 2,185.48 seconds and 11.90 GiB peak process-tree RSS.
It was co-scheduled with QfO and Three Kingdoms, so that wall time is not a
controlled runtime comparison. The independent first fresh run remains the
best estimate of high-sensitivity cost: 2.25x the optimized standard runtime
and 4.2x its peak RSS. High sensitivity therefore remains opt-in while its QfO
validation is in progress.

Cached QfO evidence remains approximately 0.674 for
`orthohmm_split_s150_c10`. No QfO search result is attributed to the new
production harness until a full production run and official assessment are
completed.

## Rejected Experiments

- Target-species batching retained fewer concurrent indices but made search
  slower. A batch limit of 16 changed search from about 5.99 to 9.00 seconds;
  a limit of 2 took 6.52 seconds. The scheduling change was removed.
- Advanced-index profile vectorization created large temporary arrays and
  regressed `n005` to 8.33 seconds search and 3.06 GiB RSS. The bounded
  `numpy.take(..., out=...)` implementation was retained.
- Independent proteome chunking was rejected because orthology is a global
  graph problem; separate chunks change the scientific result.
- A dense NumPy/Numba similarity matrix was rejected as both quadratic and
  scientifically inappropriate for variable-length unaligned proteins.
- A pure-Python MCL implementation was rejected in favor of native Leiden and
  the existing native MCL option.
- Async local FASTA I/O was rejected because this pipeline is compute bound
  and consumes files sequentially.
- The proposed generic NumPy and fixed-size Cython HMMs were rejected. They do
  not implement the existing profile model safely, and OrthoHMM already uses
  C/OpenMP, AVX2, Numba, and optional CUDA kernels.
- Fixed wall-clock assertions were rejected as hardware dependent. The harness
  records relative measurements while unit tests enforce numerical and
  representation parity.
- Relaxing strict profile self-thresholds or sequence anchors reduced the best
  comparable OrthoBench F-score from 70.4 to 69.7 or lower.
- Multiplying RBNH thresholds by 0.90-1.05 did not beat the unmodified
  reciprocal-best boundary. The factor-1.00 screen was best at F=69.8 before
  strict profile expansion.
- A global BLOSUM45 search scored F=65.9 after refinement. A wider BLOSUM62
  extension band scored F=67.7. Cross-matrix score mixing is not defensible
  without matrix-specific calibration.
- A second strict-profile iteration doubled candidates and profile builds,
  increased replay wall time from about 319 to 721 seconds, and reduced the
  final score from F=70.4 to F=70.2.
- A reduced-alphabet prefilter lost 10,168 exact-alphabet significant pairs in
  the screened species pair while adding 2,073. Running it only for exact
  no-hit queries rescued 110 of 14,300 queries and 122 significant pairs, too
  little evidence for its cost and specificity risk.

No unmeasured 15-25x claim is supported. Only the controlled measurements
above should be reported.

## Recommendation

Adopt the compact representations, 4 x 8 scheduler on 32 CPUs, deterministic
direct Leiden path, numeric thresholds, and indexed orthogroup materialization
as production defaults. Retain the corrected high-sensitivity workflow as an
opt-in profile: its fresh OrthoBench checksum and accuracy gates pass, but it
has a substantial compute and memory cost and does not yet beat OrthoFinder.
Do not adopt multi-matrix consensus, relaxed profile boundaries, iterative
profiles, or reduced-alphabet rescue from the present evidence. Complete the
fresh QfO assessment before considering high sensitivity as the default.
