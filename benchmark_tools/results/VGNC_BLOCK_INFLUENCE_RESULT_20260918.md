# VGNC Fixed-Table Block Influence Result

## Execution And Validation

Executed the [exploratory protocol](VGNC_BLOCK_INFLUENCE_PROTOCOL_20260918.md)
with implementation713b3f7, committed and pushed before inspecting results.
All four historical scored tables passed the existing source/dependency
audit and input hashes were checked again after execution. No inference,
native scorer eligibility or scientific parameters were changed.

The complete table contains67,376 deletions across16,844 reference blocks
and four stages. All remaining precision/recall/F1 ratios were defined.
The report is [vgnc_block_influence_20260918.json](vgnc_block_influence_20260918.json),
SHA-256`2d2233c880fc2688978a4cf4f163aa127d5055421fe39bd02bf3d892d97f33f7`.
Its7,570,440-byte table remains at
`benchmarks/work/vgnc_block_influence_20260918/all_block_deletions.tsv`,
SHA-256`341a828244ed3ae7556535adb76a4ce92cf325712f38235fed98fe576eff6301`.
Protocol SHA-256`f2a623f686cf27720a4a3aa99efafd797f95e9524a5acdc96171ecdf8674c3a9`.

Seventeen focused tests pass. In addition to synthetic edge cases, the
retained-result test recomputes every TSV metric and contrast range and
directly excludes raw rows for the ten largest absolute effects per stage,
independently of the production incident-count accumulator.

## Observations

Differences below are F1 percentage points. The deletion ranges are NOT
confidence intervals.

| Contrast | Full-table difference | Single-block deletion range | Sign retained |
|---|---:|---:|---:|
| multipass_refined - multipass | +39.2002 | +39.0496 to +39.2274 | 16844/16844 |
| strict_profiles_refined - strict_profiles | +39.2637 | +39.1141 to +39.2909 | 16844/16844 |
| strict_profiles - multipass | -0.1149 | -0.1739 to -0.0567 | 16844/16844 |
| strict_profiles_refined - multipass_refined | -0.0514 | -0.0738 to -0.0335 | 16844/16844 |

The largest absolute individual-stage influence is the OR10A5 block:
deleting its982FP and6TP increases multipass F1 by0.1573percentage points;
deleting its984FP and6TP increases strict_profiles F1 by0.1563points.
For each refined stage the largest effect is HSPA4, removing30FP and6TP
and increasing F1 by approximately0.0269points. These are influence
diagnostics, not biological validation of these assignments.

No single block deletion reverses the four observed differences. This does
not imply statistical significance, robustness to deleting multiple blocks,
generalization, or superiority over competitors. Shared errors make these
deletions dependent. The native eligibility rules were not recomputed;
therefore these are not benchmark results for reduced reference datasets.
The historical stages are not substitutes for pending corrected QfO runs.
VGNC paired uncertainty remains unresolved.
