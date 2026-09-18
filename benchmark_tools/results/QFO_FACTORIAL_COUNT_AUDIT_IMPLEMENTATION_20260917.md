# QfO Factorial SwissTrees Count Collector

Implemented `benchmark_tools/audit_qfo_factorial_swiss.py` to bridge the eight
independently admitted factorial assessments and the frozen 42-endpoint
bootstrap. No empirical eight-cell count report has been generated yet.

## Required Evidence

Supply a JSON manifest containing `cells`, exactly in binary P,C,R order
(`p0_c0_r0` through `p1_c1_r1`). Each entry contains `cell` and `admission`;
the latter is a file record with absolute `path`, `bytes` and `sha256`.
Pass the manifest hash explicitly. Assemble and commit this inventory only
after all relevant admissions exist; missing or failed cells are not imputed.

Indices 0 and 4 use their `admitted_reused_assessment` reports. The collector
requires their original frozen admission hash, exact original admitted-stage
record, participant identity, unchanged partition, pairs and coverage, and
the explicit no-rerun marker. Other indices require
`fresh_factorial_assessment_admitted`, successful saved scheduler evidence,
and the expected factorial participant. This collector consumes admissions;
it does not replace or independently repeat their scheduler/provenance gates.

All nested file records in the admissions and their execution reports are
hash-checked before and after assembly, deduplicated by path; conflicting
identities fail. Raw SwissTrees files must occur uniquely in the original
execution output inventory. Uninventoried raw files are not accepted.

The existing four-stage count audit is pinned to SHA-256
`546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183`.
It provides the previously audited native reference orientation, family
inventory and raw reference universe. Every cell must reproduce exactly the
same canonical reference pair identities, truth labels and represented genes,
as well as per-family relation counts and mapped protein counts. Matching
aggregate counts alone is insufficient.

Native counts are reconstructed as raw relation counts / 2 + 1. Every
family precision/recall and macro aggregate must reproduce the native metric
records within 5e-8. The final eight-cell structure must also pass the frozen
bootstrap input validator. F1 is harmonic macro precision/recall, not mean
family F1. The audit adds no intervals or biological-independence claims.

## Validation and Use

Eighteen new tests cover synthetic eight-cell counts, actual saved baseline
binding, rejected failures and wrong participants, changed family coverage,
duplicate pairs, membership changes, wrong native scores, truth-label swaps
with unchanged counts, full synthetic file assembly, manifest hash mismatch,
and raw-file tampering after admission. Combined with the count-reader and
factorial-bootstrap tests, 42 tests pass.

After all admissions are available, from the repository root:

```bash
python benchmark_tools/audit_qfo_factorial_swiss.py --manifest ADMISSION_INVENTORY.json --manifest-sha256 INVENTORY_SHA256 --baseline-counts benchmark_tools/results/qfo_swiss_counts_20260917.json --output NEW_COUNTS.json
python benchmark_tools/bootstrap_qfo_factorial.py --counts NEW_COUNTS.json --counts-sha256 COUNTS_SHA256 --protocol benchmark_tools/results/QFO_FACTORIAL_PROTOCOL_20260917.md --output NEW_INTERVALS.json --markdown NEW_INTERVALS.md
```

The placeholders denote future audited inputs, not existing empirical results.
Both tools require fresh output paths. Preserve the inventory, count report,
interval report and all original evidence. This SwissTrees workflow does not
provide uncertainty for GO, EC, FAS, VGNC, TreeFam-A or the secondary mean.
