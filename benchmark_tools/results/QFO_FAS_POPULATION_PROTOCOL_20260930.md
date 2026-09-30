# Retained FAS Population Recount Protocol

Prospective diagnostic, frozen before examining recounted strata or population
score sums. No scientific default, submitted prediction, native benchmark
score, primary endpoint or prior result is changed.

## Question

Do the corrected eight-method FAS logs' precomputed, missing and unannotated
counts reproduce from the retained prediction databases? What deterministic
range contains a hypothetical complete eligible-pair mean if every currently
uncomputed score is in [0,1]?

The existing [saved-sample audit](QFO_CORRECTED_FAS_SAMPLE_AUDIT_20260928.md)
and [sample-attrition bounds](QFO_FAS_SAMPLE_ATTRITION_20260928.md) do not recount
the eligible populations. The
[lookup exposure check](FAS_SAVED_COMPLEXITY_EXPOSURE_20260928.md) validates
only scores relevant to saved samples. This new pass has a distinct purpose:
all canonical precomputed scores and all retained native-query populations.

## Frozen Inputs And Rules

- Use all eight admitted methods in corrected v7 manifest SHA256
  `042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc`.
  Keep the published method order; do not omit the large sequence-only database.
- Reuse the retained attrition logs/commands/raw samples and complete annotation
  panel. Validate their recorded files and counts; do not re-run annotations.
- Use the frozen precomputed lookup: 590,454,327 compressed bytes, SHA256
  `6904475bcc4b623305c4fb3105b5e20b0a0959363000abd1e43f2fbf26bb005f`.
- Reproduce the native DISTINCT accession query literally. No new taxon,
  gene/family, confidence, pair-count or score filter. Skip accession aliases
  with underscores exactly as the scorer does. Precomputed membership takes
  precedence over annotation eligibility.
- Stream the complete lookup once. Use compact uint64 accession-pair keys;
  stable canonical overwrites retain the last valid entry, as in the native
  loader. Reject duplicate exact JSON keys or unsupported shapes, nonfinite
  scores and out-of-range scores rather than silently changing their treatment.
- Read each database through SQLite `mode=ro`, `query_only=ON`, file-backed
  temporary sorting and a 64 MiB page cache. Reject live journal/WAL sidecars.
  No database or reference is modified.
- Newly hash database bytes before/after scanning. These are current retained
  identities, **not independently established historical database checksum
  pins**. Task commands link their paths; exact recount/log agreement strengthens
  consistency but does not by itself prove upstream conversion or predictions.
- Recount all DISTINCT query pairs, alias skips, precomputed pairs, annotated
  uncomputed pairs and unannotated pairs. Require precomputed/missing/no-feature
  counts to equal native logs, and their eligible sum to equal the immutable
  corrected comparison's FAS assessed-relations count. Check all saved sample
  lookup classifications and precomputed values as an additional cross-check.

## Bounds And Interpretation

For the fixed eligible population let P be the precomputed count, M the
uncomputed count, and S the sum of its precomputed scores. The hypothetical
completed mean is in `[S/(P+M), (S+M)/(P+M)]` if every remaining score lies in
[0,1]. Summation uses `math.fsum` within batches and across batch sums; this is
double-precision numerical evaluation, not symbolic exact arithmetic.

These are deterministic conditional sensitivity bounds, **not confidence
intervals**, new native FAS scores, evidence that missing scores are defined,
family-independent validation or a generalization estimate. There is no
pair-IID resampling or missing-at-random assumption. Precomputed score validity,
annotation accuracy, historical RNG state, native sampling-law verification,
shared-pair covariance and other QfO uncertainty requirements remain separate.
Keep all eight bounds irrespective of whether they preserve any ranking.

## Execution And Validation

One local analysis allocation: one CPU, 16 GiB RAM, 90-minute maximum, no GPU,
exclusive-host claim, requeue, retry, new inference or scientific timing result.
Use the existing Python environment with NumPy 2.2.6 and ijson 3.5.0; do not
install or upgrade shared packages. Set numerical-library thread counts to one.
Only this job may be polled/managed; do not signal unrelated work or access DGX.

```sh
/home/bizon/anaconda3/bin/python -B -m benchmark_tools.audit_fas_population \
  --output benchmarks/work/qfo_fas_population_20260930
```

The directory must be fresh. Per-method outputs preserve partial evidence;
a final report requires all eight methods, all recount checks and stable
input/parser source hashes. On discrepancy/failure, stop and retain the
allocation/log/partial artifacts; do not retry to obtain a passing score.

73 focused tests pass, including 26 new tests. They cover collision-free pair
encoding, range checks, literal native-query equality, read-only DISTINCT alias
behavior, exact native-loader agreement on canonical overwrite/invalid-value
fixtures, duplicate-key rejection, empty/no-eligible populations, annotation
precedence, score-domain rejection, and batch-invariant counts/bounds. These
are component/synthetic tests, not the actual eight-method population result.

The helper and protocol must be committed and pushed before allocation launch.
After terminal completion, inspect scheduler exit/accounting, retain input/
output pins, independently reconstruct counts/bound arithmetic, and report
any discrepancy. No output may be promoted to full QfO uncertainty admission
or publication readiness.
