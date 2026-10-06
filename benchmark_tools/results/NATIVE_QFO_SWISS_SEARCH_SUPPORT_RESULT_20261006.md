# Native SwissTrees Search Support Independently Checked

## Actual Finding

The [complete direct-search trace](native_qfo_swiss_search_support_20261006_v1.json)
examines all 2,023 removed SwissTrees reference relations from the preceding
[reconciliation localization](NATIVE_QFO_SWISS_RECONCILIATION_TRACE_RESULT_20261006.md).
It traces 361 selected genes through the original significant-search checkpoints
of native P0/C0/R0 and P0/C0/R1. Both checkpoints contain 984,137 genes and
90,687,327 directed hit records. All six checkpoint files in each view were
inventoried by the original scientific output validation, unlike the newly
observed tree/node artifacts in the preceding localization.

For the 334 removed true positives, each view has 60 pairs with no direct hit,
6 with one direction and 268 with both directions. For the 1,689 removed false
positives, each view has 479 pairs with no direct hit, 20 with one direction and
1,190 with both directions. Thus 274 removed TP and 1,210 removed FP have at
least one retained direct hit. The [generated case ledger](native_qfo_swiss_search_support_20261006_v1.tsv)
retains every relation, including absent directions and directional counts.

Every selected pair has exactly equal directional score multisets between the
two views, preserving multiplicity and ignoring execution-dependent hit-row
ordering. Each view has 2,942 selected directed records. The different
whole-array hashes remain recorded: selected-pair agreement is not whole-hit-set
equality or proof of identical upstream search histories. The generic trace
preserves duplicate instances; the actual selected evidence has no duplicate
row instances.

These checkpoints store significant homology-search hits before RB-NH graph
construction. They are not graph edges, cluster-derived orthology pairs or raw
prefilter candidates. Stored scores are not calibrated orthology probabilities;
per-hit E-values are not present and are not reconstructed. Absence of a direct
hit does not identify whether a prefilter omitted a candidate, scoring was
insignificant or grouping used indirect connectivity.

Together with the preceding complete exclusion trace, this distinguishes
direct search support from downstream duplication-rule orthology exclusions.
Both true and false positives can have bidirectional homology support before
the later exclusions. It does not establish true duplication history, correct
inferred trees, the biological appropriateness of an exclusion or calibrated
specificity. Initial HMM search is on in both cells; this is not a total-HMM
ablation or a sensitivity-matched non-HMM comparison.

## Independent Check

The [independent readback](native_qfo_swiss_search_support_code_readback_20261006.json)
does not import the primary lookup. The primary uses a selected-endpoint mask
and Python directed-pair membership; the readback uses bounded int64 directed
keys and sorted-key binary lookup. Both read-only mmap scans validate array
shape/dtype, all endpoint bounds and finite scores without constructing an NxN
matrix or materializing all hit rows as Python objects.

The actual independent scan checks all 181,374,654 checkpoint records,
all 2,023 selected pairs, all 361 target genes, every 5,884 selected directed
record, and complete absent directions. It reproduces every selected row
offset and unmodified score, all 12 summary cells, exact directional multiset
equality and the entire TSV inventory. Original validation/manifest/file pins
and directly supplied artifacts are checked before and after the scan. These
are direct diagnostic checks, not a new transitive scientific admission or
independent biological confirmation.

The primary JSON is 1,999,028 bytes with SHA256
`787bde01ffacd60547fb290bdbb1d04242230ed05ce77d25d270227310252d5a`.
The case TSV is 139,368 bytes with SHA256
`2854fed4328874e0dcf55fe7dbc4e5d93c309e6de805b1210100467eef42b74b`.
The independent readback is 8,585 bytes with SHA256
`8ffc58fe192fa060de4e068161e72b8ab626066446608e5c0011925e62c67ede`.

306 joined tests pass in 11.94s, with no failures, errors or skips. The new
modules contribute 35 primary and 26 independent tests, covering exhaustive
tiny-array oracles, duplicates/absence, invalid arrays/controls, original-review
binding and corrupted full-readback fixtures. Those fixtures use synthetic
checkpoint inventories, not production admission. One retained-receipt test
checks the actual report bindings/results without repeating the large scans.
The joined suite also covers prior reconciliation, transitions, native counts,
uncertainty binding and score export; it does not include the earlier functional
pair diagnostics. Preserve both the earlier 305-case XML and final 306-case
XML rather than overwriting test history.

Actual analysis uses original Python 3.10.13/NumPy 2.2.6; tests use Python
3.12.3/NumPy 2.2.6. GNU time records 8.21s/1,464,672 KiB maximum RSS for the
primary and 13.60s/1,499,152 KiB for the independent scan, both exit 0 with
zero swaps. These are shared-host postprocessing observations, not inference
timings or a speed ranking. Competing CPU, memory-bandwidth and I/O demands
may distort them by an unknown, potentially tool-dependent amount. The
pre-diagnostic capacity check found about 771 GiB available RAM and 9.5 TiB
available disk; full swap occupancy was observed and disclosed.

## Reproduction And Limits

Follow the [retrospective protocol](NATIVE_QFO_SWISS_SEARCH_SUPPORT_PROTOCOL_20261006.md)
and retain locally available original native/checkpoint artifacts. From the
repository root, select fresh output paths:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/trace_native_qfo_swiss_search_support.py \
  --localization benchmark_tools/results/native_qfo_swiss_reconciliation_trace_20261006_v1.json \
  --localization-sha256 53c692cf7dab9affac7b0240cf5191d2530bd4369b66f4165e8acb811391dacb \
  --output benchmark_tools/results/native_qfo_swiss_search_support_replay.json \
  --ledger benchmark_tools/results/native_qfo_swiss_search_support_replay.tsv

env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/readback_native_qfo_swiss_search_support.py \
  --report benchmark_tools/results/native_qfo_swiss_search_support_20261006_v1.json \
  --report-sha256 787bde01ffacd60547fb290bdbb1d04242230ed05ce77d25d270227310252d5a \
  --output benchmark_tools/results/native_qfo_swiss_search_support_readback_replay.json
```

The second command checks the retained result, not the newly created replay.
These commands are diagnostic reproduction, not whole-study archive restoration.
The frozen core, original inference/scoring/defaults and historical manuscript,
PDF and archive remain unchanged. Incorporate this companion evidence during
the next complete manuscript assembly, without claiming the archived component
already includes it. Failed R1 timing remains failed/ineligible. No native
inference, graph construction, alignment, tree inference, reconciliation,
endpoint scoring, new confidence interval or accuracy admission is performed.
No unfinished index8 output is read, no original job is restarted and no
successor launch is authorized by this evidence.

Further graph-stage/error-strata/tree-error work, missing native cells, broader
uncertainty and independent generalization/provenance/release requirements
remain open. The full publication goal is active and completion unproven.
