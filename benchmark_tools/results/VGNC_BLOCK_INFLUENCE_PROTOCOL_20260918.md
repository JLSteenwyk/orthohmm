# VGNC Fixed-Table Block Influence

## Scope

Exploratory diagnostic specified after inspecting historical aggregate scores
and dependency structure, but before inspecting these deletion results.
This is not confirmatory inference, a bootstrap, or a jackknife variance.

Use the four frozen historical scored tables already checked by
`audit_vgnc_family_dependencies.py`, in order: multipass,
multipass_refined, strict_profiles, strict_profiles_refined. Construct
reference blocks by merging family labels sharing a reference protein.
Preserve every block, including those with no incident scored rows.

For each stage and each block separately, remove every scored TP, FP or FN
row touching that block. Remove a within-block row once. A cross-block FP
is removed once in each endpoint's separate deletion, never twice in one
deletion. Recompute precision, recall and F1 from remaining integer counts.
Fail on duplicate/conflicting pairs, negative counts or undefined ratios.
Retain the complete deletion table and summarize the ten largest absolute
F1 changes per stage with lexicographic tie-breaking.

Apply the same deletion to both stages for four contrasts: 1-0, 3-2, 2-0,
3-1. Report full-table differences, minimum/maximum deleted differences,
and positive/zero/negative counts over all blocks. Do not select contrasts
or exclusions after inspection. Hash input files before and after execution.

## Interpretation

This deletes rows from fixed scored tables. It does not rebuild the reference,
rerun inference or recompute the native scorer's eligibility rules. Deletions
are dependent and their ranges are not confidence intervals. Stability to a
single deletion does not establish stability to joint deletions, new datasets
or alternate references. It cannot establish a causal effect of refinement.
These historical stages do not replace corrected-release or competitor runs.

## Execution

```bash
python benchmark_tools/diagnose_vgnc_block_influence.py \
  --output benchmarks/work/vgnc_block_influence_20260918
```

Output must be fresh. The full TSV stays in the work directory; retain its
SHA-256 and the JSON summary in the publication evidence inventory.
