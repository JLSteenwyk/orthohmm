# Native Candidate Group Trace and Failed Identifier Join

The [prespecified complete grouping trace](NATIVE_QFO_CANDIDATE_GROUP_TRACE_PROTOCOL_20261006.md)
is implemented and tested, but its selected export fails at accession
localization. No successful pair-path ledger/report exists. The original
[failure](native_qfo_candidate_group_trace_20261006_v1/failure.json) is retained,
not overwritten or automatically retried. Tested primary/independent reader
83e87c8d is pushed before that execution. The original114-test suite passes
9.72s; its earlier synthetic test failure exposed harmless repeated inventory
entries, not a scientific failure. Separate conflicting entries are rejected.

The [follow-up protocol](NATIVE_QFO_CANDIDATE_GROUP_ID_DIAGNOSIS_PROTOCOL_20261006.md)
and original failure are pushed as1d468e27 before inspection of unmapped counts.
Tested diagnostic496aeaa7 is pushed before its single actual execution. It
uses the separate stdlib partition parser/undirected graph traversal, not the
primary tracer or its accepted-union kernel. The81-test diagnostic/export
suite passes1.41s before execution. Both original native output inventories,
conversion/terminal/admission links, prior decomposition/readback and source
identities are checked before/after reading. No transitive scientific
admission, inference, conversion, endpoint scoring or raw database audit is
repeated.

Final127joined tests pass9.99s with zero failures/errors/skips, including
actual output/source identities and manuscript/claims/result limitations.

## Actual Diagnosis

| Quantity | Retained Result |
| --- | ---: |
| Complete baseline/candidate gene universe | 984,137 |
| Baseline groups | 394,328 |
| Candidate groups | 353,638 |
| Accepted unions | 40,690 |
| Scored transition union rows | 42,080 |
| Changed scored pairs | 2,295 |
| Changed pairs without missing native accessions | 2,287 |
| Changed pairs with unmapped accessions | 8 |
| Affected candidate TP / FP pairs | 0 / 8 |

The independent graph reconstruction exactly matches the WHOLE candidate
partition, including untouched groups and all genes;394328-353638=40690.
Every transition row is unordered-pair canonical and unique. The eight join
failures concern two accessions,`Q17QN5_BOVIN` and`Q1RMT5_BOVIN`, absent from
the injectively normalized complete native universe. All eight are
`not_scored` to`FP` transitions. The entire affected-pair table is retained,
not a favorable subset. No duplicate/order defect explains this failure.

This diagnoses a join boundary, not the reason for biological false positives.
It does not establish why those identifiers differ, identify a valid alias
replacement, rescore the eight pairs, change the original162TP/2133FP
arithmetic or establish identical upstream search histories/numeric seed IDs.
The2,287 other changed pairs are NOT admitted as a complete localization.
No suffix stripping, guessed alias, approximate join or default change is
used. Pair-path localization remains explicitly incomplete. A future mapping
investigation must use original inventoried mapping evidence and gene/protein
identities, not infer equivalence from name resemblance or tune on scores.

## Identities

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Original trace failure | 1,076 | `bfad1eefd8f084adc6ab66d5da385c9ce8a1f17440948c59330cf1a761efb695` |
| Independent identifier diagnosis report | 14,723 | `3bc6646e1dc7cb6488ceb6aaaf763e3a4e3cdcbd41268f2066ef2d8008f5be74` |
| Complete unmapped-pair table | 454 | `6ec0ed60e8004f99ed8822d09ad87da665ddc25883bcf9ffd1b8c073400073b7` |

[Report](native_qfo_candidate_group_id_diagnosis_20261006_v1/report.json)
and [all affected pairs](native_qfo_candidate_group_id_diagnosis_20261006_v1/unmapped_pairs.tsv)
retain original inventories, exact source/input identities and false
localization/scoring/CI/readiness flags. Original failure source remains
unchanged. The primary18.30s/1,227,516KiB/exit1 and diagnosis10.85s/
1,256,952KiB/exit0 observations are local shared-host postprocessing, not
native timing, speed comparisons or isolated efficiency evidence. Each
records zero swaps; available memory660104220KiB was observed before launch.
Full swap is also retained as a host observation, not a quietness certificate.

## Reproduction

Run tests in the existing test environment; the original selected export
namespace remains immutable. Reproduce ONLY the separate diagnosis in a new
directory, not by overwriting the original failure or invoking a native retry:

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/diagnose_native_qfo_candidate_group_ids.py \
  --failure benchmark_tools/results/native_qfo_candidate_group_trace_20261006_v1/failure.json \
  --failure-sha256 bfad1eefd8f084adc6ab66d5da385c9ce8a1f17440948c59330cf1a761efb695 \
  --output NEW_FRESH_DIAGNOSIS_DIRECTORY
```

Original scientific/allocation sources and results remain unchanged. Failed
R1 timing remains ineligible. Existing23902/23910 continue; this read-only
failure diagnosis supplies no new score, admission, native retry, CI or
publication-readiness claim. All full-goal requirements remain in force.
