# Native Duplication-Annotation Precision-Recall Trade-Off

[Protocol](NATIVE_QFO_SWISS_DUPLICATION_STRATA_PROTOCOL_20261006.md) pushed
at91ca7385 before outcome projection; first implementation1b893235. Initial
actual run fails before scoring/output by comparing tree-leaf aliases directly
with native accessions. [Correction](NATIVE_QFO_SWISS_DUPLICATION_ALIAS_AMENDMENT_20261006.md)
and final source/tests10822eea pushed before successful computation. Failed
receipt and pre-correction tests remain retained.

[Generated table](native_qfo_swiss_duplication_strata_20261006_v1/TABLE.md),
[report](native_qfo_swiss_duplication_strata_20261006_v1/report.json) and
[independent readback](native_qfo_swiss_duplication_strata_readback_20261006_v1.json)
extend native P0C0, not selected-default comparator scores. Initial HMM search
on, profiles/candidate expansion off. R0 group-clique pairs, R1 resolved native
pairs; failed R1 timing remains ineligible.

| Fixed Annotation Bin | Families | R0 F1 (%) | R1 F1 (%) | F1 Change (pp) | Precision Change (pp) | Recall Change (pp) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| All | 18 | 68.918 | 78.957 | +10.039 | +30.521 | -6.533 |
| Lower, at/below7/48 | 9 | 69.876 | 84.167 | +14.291 | +32.615 | -1.598 |
| Upper, above7/48 | 9 | 67.957 | 73.577 | +5.620 | +28.427 | -11.468 |
| Missing | 0 | NA | NA | NA | NA | NA |

Upper-bin descriptive recall penalty is larger and F1 gain smaller, despite
precision increasing in both bins. No significant interaction, causal mechanism
or true evolutionary-history claim. Fractions are reference-derived and related
to the same orthology labels; size/composition confounding and development
exposure remain. Default-S is not explicit speciation. No new uncertainty,
bootstrap draws or transferred historical intervals.

## Checks And Execution

Primary Fraction/median logic and independent integer-cross-product bin decisions
reproduce unchanged September23 median/ties/memberships. Both re-read original
40,019,520-byte compressed identifier mapping and require exact reference-entry
set equality for18 families/563 native genes. Retained Ensembl/other leaf aliases
are not compared literally with accessions; HOX aliases deduplicated by entry ID.
Informative nodes, not unique proteins minus one, remain the denominator.

Both re-read21,530 native raw rows, identical truth/complete memberships, all36
family records, eight projections/four differences and every TSV cell, including
null/NA missing-bin values. Alternative macro precision/recall and harmonic F1
agree, not pooled pairs or mean family F1. Direct checks repeated after reads.
Eight other original feature-source records inherited, not rechecked. Original
Darwin traversal/acquisition and whole-study admission not repeated; same-data
implementation consistency is not biological replication.

Final [joined JUnit](native_qfo_swiss_duplication_strata_20261006_v1/native_qfo_swiss_duplication_strata_joined_tests_20261006_v2.xml)
reports131 tests, zero failures/errors/skips (1.830s XML; displayed1.88s),
including50 new cases. Exact odd/even medians/ties, missing/empty bins, alias
denominators/entry joins, wrong IDs despite equal counts, invalid vectors,
row coverage/semantics/statistics and output refusal tested. Earlier40/121-case
receipts are pre-correction history. Initial scientificPython3.10 test invocation
fails because pytest is absent; tests use retained testPython3.12, no install.

Actual scientificPython3.10 runs remove PYTHONPATH/PYTHONHOME/PYTHONUSERBASE/
LD_PRELOAD/LD_LIBRARY_PATH/LD_AUDIT; set PYTHONNOUSERSITE=1,
PYTHONDONTWRITEBYTECODE=1, PYTHONHASHSEED=0, OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1. Independent reader uses-I-S-B/stdlib only.
[Primary receipt](native_qfo_swiss_duplication_strata_20261006_v1/native_qfo_swiss_duplication_strata_export_20261006_v2.time.txt):
7.10s/1,573,724KiB maximum process RSS.
[Reader receipt](native_qfo_swiss_duplication_strata_20261006_v1/native_qfo_swiss_duplication_strata_readback_20261006_v1.time.txt):
6.52s/1,550,024KiB. Both exit0/zero swaps. Shared-host postprocessing, not
inference costs or isolated speed; CPU/memory-bandwidth/I/O contention effects
unknown, potentially tool-dependent. Pre-run available RAM649,647,836KiB;
near-full swap25,316KiB free. No unrelated-work modifications/contention retry.

```bash
benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.export_native_qfo_swiss_duplication_strata \
  --repo . --output /absolute/fresh/native-duplication-projection
benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/readback_native_qfo_swiss_duplication_strata.py \
  --report /absolute/fresh/native-duplication-projection/report.json \
  --report-sha256 ACTUAL_REPORT_SHA256 --output /absolute/fresh/readback.json
```

Apply the sanitized environment above. Pinned original paths must be available;
not a portable final archive or full transitive workflow certification.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Exporter | 13728 | `8b2e0ab7e85e6d33fc8cee439245ec8e13c0698c9d41ebf9041c012774757e0f` |
| Independent reader | 11383 | `1a362d6a406e0bd29e865741d5363b58fd17ba69b923cd8a47ed010b380a1800` |
| Report | 41397 | `e0ecbc9d26738d0d3d573ca31772085b28d145ba2deaa9383c8cff75bd0d67b5` |
| Readback | 12310 | `84eb3d9fe02b5984cecaa06355d1cb1d68c54aec09ab588eef8b3ada42f1957b` |

## Remaining Scope

Advances4.3, not full validated divergence/duplication-history strata, broader
uncertainty, matched-search/HMM interactions, generalization or final package.
Five native QfO score cells unavailable. Original22444.0 completes0:0 at11:18:09;
22444 batch RUNNING11:20:18, original22445/22450/22451/22452 pending. Native-step
completion does not prove whole-run success. No unfinished output read or next
identity released before actual terminal/review gates. No original job restarted;
historical manuscript/supplement/archive/frozen bytes unchanged. Full goal active.
