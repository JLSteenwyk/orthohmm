# FastOMA Native Pair Uniqueness Audit

The complete historical QfO FastOMA native output passed strict row validation
against all 976,504 input accessions across 78 species. The disk-backed converter
and an independent GNU sort of the retained prediction files produced identical
raw and QfO-filtered pair sets.

| Quantity | Count |
| --- | ---: |
| Native rows | 15,320,615 |
| Distinct cross-species pairs | 15,320,615 |
| Duplicate relations, including reversed duplicates | 0 |
| Retained by historical QfO mapping | 15,277,489 |
| Removed for identifiers absent from mapping | 43,126 |

The native gzip was read to completion, including integrity verification. Every
row had exactly two known, distinct, cross-species accessions. Input files and
source records were checked before and after conversion. The old pair files
were not edited. Distinct native pairs have SHA-256
`efeb607b108c162541656b12676df252706ed5902b5dce662e4afa2c587db16b`;
filtered pairs have SHA-256
`4cb9c0f21fc5c46e15b891696d902bd4efe6af7249fcea74a33ef7162185e217`.
These match the independently sorted retained files byte for byte.

Machine-readable evidence: `fastoma_distinct_pair_audit_20260918.json`, SHA-256
`52c39c666d0d66cfe8886e2f48277244d1190ec668b989eea2f1024fe902edb5`.
Large generated pairs remain in
`benchmarks/work/fastoma_distinct_pair_audit_20260918/`, outside the commit.

Reproduce with fresh output paths:

```bash
python benchmark_tools/audit_fastoma_distinct_pairs.py \
  --root . --work benchmarks/work/fastoma_distinct_pair_audit_fresh \
  --output benchmarks/work/fastoma_distinct_pair_audit_fresh.json
```

This closes the historical converter's uniqueness limitation, not native
workflow-completeness or biological-correctness questions. No new accuracy
scores were computed, and this is not a corrected-release result. The mapping
exclusions remain explicit; the audit does not establish their absence of bias.
Shared-host execution was not a matched efficiency experiment. The corrected
FastOMA tree, fresh launch, native admission and scoring remain unfinished.

Validation: 41 focused tests passed, including the original streaming converter,
disk-backed conversion and independent sorting comparison. The complete real-data
audit also passed; no prediction-set changes or benchmark endpoint changes were
needed.
