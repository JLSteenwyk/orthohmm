# Native Candidate VGNC Pair Changes

## Actual Findings

The [new decomposition](native_qfo_candidate_vgnc_20261006_v1/report.json)
uses the independently admitted P0C1R0 output of assessment23894 and the
previously verified P0C0R0 complete pair table. Both configurations retain
initial HMM search and submit group-derived clique pairs with reconciliation
off. This is not a total-HMM ablation or a selected-default tool comparison.

| Configuration | TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P0C0R0 | 19,981 | 16,013 | 3,953 | 0.555120 | 0.834837 | 0.666834 |
| P0C1R0 | 20,143 | 18,146 | 3,791 | 0.526078 | 0.841606 | 0.647445 |

The complete [42,080-pair transition table](native_qfo_candidate_vgnc_20261006_v1/pair_transitions.tsv)
has exactly these five nonzero states:

| Baseline State | Candidate State | Pairs |
| --- | --- | ---: |
| FN | FN | 3,791 |
| FN | TP | 162 |
| FP | FP | 16,013 |
| TP | TP | 19,981 |
| not_scored | FP | 2,133 |

Candidate expansion recovers162 asserted TPs and adds2,133 scored FPs, retaining
all baseline scored TPs/FPs. No other scored transition is observed. Precision
falls2.904231 percentage points, recall rises0.676861 points and F1 falls
1.938806 points. This accounts for the pooled-score arithmetic; it is not a
confidence interval, causal biological mechanism or universal orthology claim.
`not_scored` is not a TN or proven biological non-orthology state. Native VGNC
eligibility and its asserted reference, not all submitted predictions, define
the scored denominator. Prediction edges were not independently rescored.

The [complete candidate sparse table](native_qfo_candidate_vgnc_20261006_v1/candidate_blocks.tsv)
contains25,507 nonzero unordered block cells,8,663 crossing reference blocks.
There are18,144 cross-block and two within-block FPs, versus16,011 and two
in the verified baseline. All2,133 added scored FPs therefore cross unchanged
reference overlap blocks. Shared-label reconstruction remains16,863 labels,
36,986 proteins,23,934 asserted pairs and16,844 blocks; these blocks do not
establish independent biological resampling units. No VGNC CI is admitted.

## Execution And Independent Checks

[Protocol](NATIVE_QFO_CANDIDATE_VGNC_PROTOCOL_20261006.md)61856584 is pushed
before raw-row/transition inspection. Aggregate outcomes were already known;
the protocol explicitly does not claim outcome-blind selection. Tested code
37a0c67f is pushed before the actual export. Initial119 tests pass2.56s;
expanded125 pass3.45s, zero failures/errors/skips, including46 new cases and
79 existing decomposition/mapping cases. New cases include complete synthetic
export/readback, aliases/shared labels/category overlap, rebound semantic
failures, database mutation, source/input/output corruption and failed-timing
relabeling. Original primary/independent kernels and14 allocation sources
remain unchanged. No package installation or scientific environment change.

The final [joined suite](native_qfo_candidate_vgnc_joined_tests_20261006_v1.xml)
passes212 cases in8.98s, zero failures/errors/skips, including the existing
three-cell figure/SwissTrees readback contracts and two new actual-artifact/
manuscript-table checks. These checks verify derived/source identity and
presentation correspondence, not another raw audit or scientific admission.
Final primary source SHA256aee858025a4528fe59edb9c33e3367c50580f5627563bb81529945aac3a1bdcc;
reader SHA25667255b6c4993c3f31eede03ae91fca8a799a9c340725f9ef11cf64ce32c7b2ad.

Primary and [independent readback](native_qfo_candidate_vgnc_readback_20261006_v1.json)
each run once successfully in original Python3.10.13. The reader uses-I-S-B,
imports no primary exporter and reuses the unchanged independent stdlib
reference/union-find, SQLite and raw-count kernels. It separately parses the
baseline CSV and all42,080 candidate raw rows, verifies every candidate sparse
cell and transition pair, and uses Fraction arithmetic for pooled ratios and
differences. All79 alias rows/no identical extra rows and the original mapping
digest match. Both routes fully hash the new prediction database before and
after reading; direct admission/execution/source/aggregate/raw inputs are bound.

The original baseline raw/database audit and its independent readback are
reused, not repeated or relabeled as newly raw-validated. Current-byte checks
do not establish uninterrupted historical integrity, complete transitive
scientific admission or completeness of omitted FPs. No inference, search,
clustering, reconciliation, scoring, FAS sampling or bootstrap is repeated.
Failed recovered R1 timing stays null/ineligible; no next identity is released.

GNU-time receipts in benchmarks/results retain primary5.01s/101,668KiB and
reader3.82s/134,008KiB, both exit0 and zero process swaps. These are shared-host
postprocessing observations, not native timing or comparative speed evidence.
Available host memory650346384KiB and near-full swap were observed beforehand.
CPU, memory-bandwidth and I/O contention have unknown potentially tool-dependent
effects. Unrelated workloads are untouched. Native23902 remains RUNNING and
review23910 remains dependency-pending; no terminal outcome is inferred.

## Identities And Reproduction

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| report.json | 15,598 | 0c83ab917417ac4e97570b8a68b4bbb0bdfd72eeb301b930bf24446e1a623aea |
| Independent readback | 1,782 | f75cc5775eb63f5a4aa9c72d844782f285d39b00481a380dca6450669c109c57 |
| candidate_blocks.tsv | 478,085 | 3a247896b9c318df4e1022ba9ed391b3f3dd77d40cc23d684c07c3990318aa58 |
| pair_transitions.tsv | 1,452,962 | 4a130a7dcbbd916bd79ffbc9ccc67631e51902e7cfa76bb7c601c9b01e0a6cf7 |

With original local inputs, run to fresh destinations. Retained outputs need
not be regenerated on continuation. The report binds reused baseline report,
readback and pair table, current snapshot and original candidate artifacts.

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.export_native_qfo_candidate_vgnc \
  --output benchmark_tools/results/native_qfo_candidate_vgnc_replay

env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/readback_native_qfo_candidate_vgnc.py \
  --report benchmark_tools/results/native_qfo_candidate_vgnc_20261006_v1/report.json \
  --report-sha256 0c83ab917417ac4e97570b8a68b4bbb0bdfd72eeb301b930bf24446e1a623aea \
  --output benchmark_tools/results/native_qfo_candidate_vgnc_readback_replay.json
```

The reader command checks the retained report. Supply the replay report and
its actual digest to check a new export. Absolute local source dependencies
remain; this is not portable full-study restoration or public deposition.
Manuscript/claims incorporate the bounded diagnostic; old figures/supplement,
rc4 archive and printed PDFs remain unchanged. Four fresh QfO scores, wider
uncertainty, original TreeFam provenance and full publication requirements
remain open. The full goal is active and incomplete.
