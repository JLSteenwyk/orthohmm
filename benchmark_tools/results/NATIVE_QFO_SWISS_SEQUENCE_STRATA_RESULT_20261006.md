# Native SwissTrees Frozen Sequence Strata

## Actual Finding

The [generated native table](native_qfo_swiss_sequence_strata_20261006_v1/TABLE.md)
projects admitted P0/C0/R0 and P0/C0/R1 onto the unchanged input-only sequence
bins frozen on 2026-09-18. These are native ablation cells, not relabeled
high-sensitivity/selected-default competitor results. All 18 families and 563
represented proteins match the frozen descriptors. The original feature
inventory, frozen native plan and both admitted conversions record identical
78-FASTA identities. This is reuse of bound input evidence, not a fresh rehash
of the entire proteome collection.

The [report](native_qfo_swiss_sequence_strata_20261006_v1/report.json) and
[TSV](native_qfo_swiss_sequence_strata_20261006_v1/scores.tsv) preserve 22 native
cell/bin rows, 11 descriptive differences and 36 individual family rows.
All retained bins are represented; empty bins have null/NA statistics, not
zero. Primary entropy bins contain nine families each. Native R1-minus-R0
F1 is +19.974 percentage points in the higher-entropy bin and +0.420 in the
lower-entropy bin. Precision differences are +37.818/+23.224 points;
recall differences are -1.749/-11.318 points, respectively.

For the seven short-relative families, F1 changes by +7.731 points and
recall by -8.753; for the other 11 families, F1 changes by +11.381 and
recall by -5.121. Full-family F1 remains 68.918%/78.957%, precision
64.394%/94.915% and recall 74.127%/67.593%. All machine values remain at
full precision; percentages are presentation only. F1 is the harmonic mean
of macro precision and recall, not pooled-pair F1 or mean individual-family F1.

These observations localize heterogeneous trade-offs in prespecified bins.
They do not establish a composition-driven mechanism, subgroup significance,
new uncertainty or superiority of selected defaults. Global entropy is not
local low complexity or evolutionary divergence. Relative shortness is not
a fragment diagnosis; no literal fragment text occurs here, and its absence
does not prove completeness. Bins overlap and can differ in taxa, architecture
and duplication history. Development exposure and small-family scope remain.
Both cells retain initial HMM search; R changes group-clique versus resolved
native-pair semantics, not only final group splitting.

## Verification

The primary recomputes bins from the frozen retained descriptors and requires
exact equality with the original memberships/cutoff/features. It binds the
existing scientific snapshot, ordinary/recovered count audits, corresponding
accuracy admissions, original feature/extractor/protocol sources and original
native plan. It re-reads both admitted compressed raw files, validating actual
per-family counts, represented members and identical reference truth.
Audited family statistics and native aggregates reproduce; original rounded
native endpoint values are preserved.

The [independent stdlib readback](native_qfo_swiss_sequence_strata_readback_20261006.json)
runs under Python `-I -S -B`, imports no export aggregation/count helper,
and independently enumerates all 21,530 raw relation rows. It uses the
equivalent doubled-count prior and inverse-reciprocal harmonic mean to check
all 36 family rows, all 22 projections, all 11 differences, empty-bin
missingness, pair semantics, frozen membership lists and the complete TSV.
Direct input/source/report/output digests are checked before and after export
and readback. This is not a repeated transitive native admission or biological
replication.

352 joined tests pass in 12.42s without failures, errors or skips. New modules
contribute 46 tests, including arithmetic oracles distinguishing pooled and
mean-family F1, unchanged membership, failed-timing preservation, complete
synthetic source-bound export/readback refusals and actual retained-result
readback. Synthetic original reviews/FASTA records are explicitly fixtures,
not production admission. Preserve the earlier 23-case and 45-case XMLs.
Prior search/reconciliation/transition/count/uncertainty/export contracts are
joined; earlier functional-pair diagnostics are outside this joined suite.

The actual export uses original Python 3.10.13 and records 0.43s/56,832 KiB
maximum RSS. Independent stdlib readback records 0.10s/19,968 KiB; both exit
0 with zero swaps. These are shared-host postprocessing observations, not
native inference timing or an efficiency comparison. Contention effects are
unknown and potentially tool-dependent. A fresh capacity check found about
769 GiB available RAM and 9.5 TiB disk; full host swap occupancy remains
disclosed. R1's failed native timing remains ineligible.

Report: 31,585 bytes, SHA256
`4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3`.
Generated table: 2,843 bytes, SHA256
`4debe245da9b0769fdfc65546ef785aeb032e55c0ee6fac7d16ac613845cc423`.
TSV: 2,125 bytes, SHA256
`876b2f97a4d775b61d1c83e9a28f5197b1a0276c59e264b33bf7d5f40184baf6`.
Independent readback: 8,718 bytes, SHA256
`8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a`.

## Reproduction And Remaining Work

Use the [prospective descriptive projection protocol](NATIVE_QFO_SWISS_SEQUENCE_STRATA_PROTOCOL_20261006.md)
and fresh output paths from the repository root:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/export_native_qfo_swiss_sequence_strata.py \
  --binding benchmark_tools/results/native_qfo_swiss_uncertainty_binding_22449_20261006.json \
  --binding-sha256 85673a114b7f7c6e05de5da189c8cbe8d77ec5d1378604dbae7e5902dfd996cc \
  --strata benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json \
  --output benchmark_tools/results/native_qfo_swiss_sequence_strata_replay

env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/readback_native_qfo_swiss_sequence_strata.py \
  --report benchmark_tools/results/native_qfo_swiss_sequence_strata_20261006_v1/report.json \
  --report-sha256 4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3 \
  --output benchmark_tools/results/native_qfo_swiss_sequence_strata_readback_replay.json
```

The second command independently checks the retained result, not the new
replay. Original locally available artifacts are required; this is not
whole-study archive restoration. Frozen method/scorer/defaults and historical
manuscript/PDF/archive bytes remain unchanged. Incorporate this companion
evidence at a complete manuscript assembly; archived components do not already
contain it. No inference, alignment/tree/reconciliation, endpoint scorer,
bootstrap or original job is rerun. No unfinished index8 output is read and
no successor launch is authorized.

Remaining native cells, matched-search contribution, graph/error/tree analyses,
broader uncertainty, independent generalization, original TreeFam mapping,
provenance, full reproducibility and release requirements remain open. The
full publication goal remains active and completion unproven.
