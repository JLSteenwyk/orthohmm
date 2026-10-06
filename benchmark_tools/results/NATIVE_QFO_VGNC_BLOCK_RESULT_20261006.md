# Native VGNC Decomposition And Paired Changes

## Findings

The [actual native report](native_qfo_vgnc_blocks_20261006_v1/report.json)
decomposes the two completed P0/C0 cells using unchanged reference-defined
overlap blocks. These are native ablation configurations with initial HMM
search on, not the historical high-sensitivity/satellite_v2 selected defaults.
All counts reproduce their admitted VGNC precision/recall/F1 within 1e-12.
The harmonic-mean admission and direct count formula differ by only 2e-16
for R1 F1; no endpoint or score was changed.

| Native Cell | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| P0/C0/R0, group-derived clique pairs | 19,981 | 16,013 | 3,953 | 0.555120 | 0.834837 | 0.666834 |
| P0/C0/R1, phylogenetically inferred pairs | 19,518 | 9 | 4,416 | 0.999539 | 0.815493 | 0.898185 |

The [complete pair transition table](native_qfo_vgnc_blocks_20261006_v1/pair_transitions.tsv)
retains all 39,947 distinct scored pairs across both cells:

| R0 Status | R1 Status | Pairs |
|---|---|---:|
| TP | TP | 19,518 |
| TP | FN | 463 |
| FN | FN | 3,953 |
| FP | FP | 9 |
| FP | not_scored | 16,004 |

No new scored pair, FN-to-TP recovery or category overlap was observed.
R1 excludes 16,004 scored FP pairs while losing 463 asserted TP pairs.
`not_scored` is deliberately not relabeled as TN or biological non-orthology.
Observed differences are precision +44.442 percentage points, recall -1.934
points and F1 +23.135 points. These are development-exposed descriptive
differences, not a confidence interval, causal evolutionary explanation or
general superiority claim.

R0 has 24,529 nonzero unordered block cells, including 7,685 cross-block
cells: 16,011 cross-block and two within-block FPs. R1 has 16,850 nonzero
cells, six cross-block cells with nine FPs and no within-block FPs. Full
[R0](native_qfo_vgnc_blocks_20261006_v1/p0_c0_r0.tsv) and
[R1](native_qfo_vgnc_blocks_20261006_v1/p0_c0_r1.tsv) sparse tables preserve
every TP/FP/FN contribution; omitted cells contribute zero counts, not a
proven lack of biological eligibility.

## Evidence And Verification

[Protocol](NATIVE_QFO_VGNC_BLOCK_PROTOCOL_20261006.md) was committed/pushed
at `e3c40104` before inspecting selected raw outcomes. Implementation and
79 passing focused tests were committed/pushed at `b6f4ac70` before actual
export. The final focused suite has 69 new tests and 10 existing mapping
tests, with no failures/errors/skips. Earlier 66-test and final 79-test XMLs
remain locally retained. Tests cover complete synthetic exports/readbacks,
shared labels and transitive overlap, aliases, within/cross-block FPs,
cross-category overlap, duplicate/truth/eligibility rejection, source and
artifact corruption, failed-timing relabeling and database mutation during
mapping. Synthetic fixtures are not scientific admission.

The original reference has 16,863 labels, 36,986 proteins and 23,934 asserted
pairs. Eleven proteins have multiple labels, resulting in eleven merged
label groups and 16,844 reference blocks. Standard-library traversal exactly
reconstructs the historical summary and the full
[reference inventory](native_qfo_vgnc_blocks_20261006_v1/reference_blocks.tsv).
Its SHA256 remains `4f9f196d61b0547e1c4f166d0528e96129cd6d882a7ef7eee245c24dc9155dc2`.
The unchanged selected mapping digest is
`caf76f99fcc6bf1207a48aee90a4ae8e34f1924f91b6471e7889d4116597f8ab`.
Both cells have 79 alias rows and no identical extra rows. Native last-label
and last-rowid alias reconstruction matches retained raw annotations; native
SQL has no general explicit alias-order guarantee.

The [independent stdlib readback](native_qfo_vgnc_blocks_readback_20261006_v1.json)
imports no exporter, reference/aggregation helper or scientific library.
It separately reads all original reference rows, selected SQLite mappings
and all 63,890 raw category rows, reconstructs blocks with union-find and
checks native eligibility, the complete asserted TP/FN partition, each sparse
cell, every reference-inventory row and every transition pair/count. Exact
`Fraction` arithmetic reproduces the pooled ratios and differences.
Primary and reader fully hash both originally inventoried prediction databases
before and after their reads, unlike the historical block analysis's partial
database checks. Original admissions, execution inventories, aggregates/raw,
snapshot/source/plan and helper identities are directly bound and checked.
This is current-byte consistency with original records, not uninterrupted
integrity proof or a new transitive scientific admission. Prediction edges
were not queried; omitted-FP completeness was not independently rescored.

Actual export/readback each completed once with exit 0 and zero process swaps
under original Python 3.10.13; reader uses `-I -S -B`. GNU-time receipts retain
5.07s/101,836 KiB and 4.99s/131,332 KiB, respectively, in
`benchmarks/results/native_qfo_vgnc_blocks_{primary,reader}_20261006_v1.time.txt`.
These are shared-host postprocessing observations, not native inference timing
or a speed comparison. Other analyses competed for CPU, memory bandwidth and
I/O, with unknown potentially tool-dependent effects. About 627 GiB RAM and
9.5 TiB disk were available; near-full host swap was disclosed. No workload
was stopped, reniced or re-affinitized. Failed R1 inference timing remains
ineligible. No inference, search, clustering, tree/scoring, bootstrap,
native-admission or completed diagnostic was repeated.

Source SHA256s:
primary `6d7b925066796c8e999e33ae87b2e5aaedbe33ac3c1129568e3494ea4d94f3bd`;
reader `0f8fc093af6ea027c739e798160508bdfdf91316d1abbb9b441d6afbb1660fe5`.
Report: 15,722 bytes,
`23cb9fd4b1c2219239ff354c97cd7dab9758667700667f88e8424414184e4c1f`.
Readback: 1,638 bytes,
`6eadf57f4f3032bbc3b787c6771824f73f47baae9c643d524b7422709ee94720`.
Pair table: 1,490,320 bytes,
`7825130427760a10fe7457ca7fd4966dffe0f3b9c07245bbcc10c6bf31a49c2d`.
Original databases and compressed raw score files remain local, not committed.
The committed tables are derived decompositions, not full prediction datasets.

## Reproduction And Remaining Work

From the repository root with the original bound local artifacts, use fresh
output paths. Default arguments pin the existing snapshot and block report.

```bash
env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE -u LD_PRELOAD \
  -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONNOUSERSITE=1 \
  PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 \
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/export_native_qfo_vgnc_blocks.py \
  --output benchmark_tools/results/native_qfo_vgnc_blocks_replay

env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE -u LD_PRELOAD \
  -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONNOUSERSITE=1 \
  PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/readback_native_qfo_vgnc_blocks.py \
  --report benchmark_tools/results/native_qfo_vgnc_blocks_20261006_v1/report.json \
  --report-sha256 23cb9fd4b1c2219239ff354c97cd7dab9758667700667f88e8424414184e4c1f \
  --output benchmark_tools/results/native_qfo_vgnc_blocks_readback_replay.json
```

The reader command checks the retained report, not the new replay. Supply
the new report's actual digest to check a replay instead. Absolute local
evidence paths remain requirements; this is not a portable full-study archive.

Neither reference overlap blocks nor method-dependent prediction components
establish independent biological sampling units. No arbitrary bootstrap/CI
is admitted; earlier failed rare-error/shared-clade uncertainty validation
remains unresolved. Native decomposition supplies correct count foundations
without resolving the sampling law. Original manuscript/PDF/archive remain
unchanged and do not yet include this companion. Remaining native cells,
matched-search/interactions, tree/error strata, wider valid uncertainty,
independent generalization, original TreeFam files, provenance and full
manuscript/release/archive scope remain open. No publication-ready claim.

Original 22444 remained RUNNING at 8:06:03; 22445/22450/22451/22452 remained
dependency-pending. Their unfinished outputs were not read, original handles
were not restarted and no next native identity was released. Frozen plan and
historical goal bytes were rechecked unchanged. Full goal remains active.
