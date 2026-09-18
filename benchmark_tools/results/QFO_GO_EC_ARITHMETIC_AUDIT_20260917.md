# GO and EC Arithmetic and Interval Audit

All 24 retained results (eight historical comparators and four recovered
stages, for each endpoint) pass the count, mean and native uncertainty-field
checks. The machine-readable report is
`qfo_go_ec_arithmetic_audit_20260917.json`. No scores were changed.

## What the Native Scorers Calculate

For each eligible predicted relation, GO/EC `avg Schlicker` selects the
maximum term similarity for each participating annotation on either side,
then averages these annotation-specific maxima across both sides. This is
not an unweighted average of the two directional averages when the numbers
of participating annotations differ. It then averages the resulting protein
pair scores equally. GO requires annotated terms with an eligible shared
ontology; EC requires EC annotations. `NR_ORTHOLOGS` is the number of scored
relations, unlike the FAS eligible-pair denominator. Neither is all predicted
relations. Native traversal is one direction by entry number, not alphabetical
accession order. The audit canonicalizes accession pairs only to detect
duplicate undirected relations; no duplicate was found.

The pinned Darwin image's `/darwin/lib/Stat` defines `StdErr` as
`InverseStudent_t(0.95, n-1) * sqrt(sample_variance / n)`: a 95% confidence
half-width under the native independent-observation model, not one SEM.
The same file implements the critical-value function. The audit executes
that function inside the pinned image for all 24 actual sample sizes and
records the program, output and source hashes. Loading `Stat` first is
necessary to make its critical-value helper available; the audit rejects
missing or malformed native output.

**The generic QfO `stderr` field has different meanings across endpoints:**
GO and EC use this Student-t half-width; the inspected FAS scorer uses a
single SEM. These fields must not be pooled or uniformly relabeled.

## Numerical Checks

Raw GO/EC scores are serialized to six decimals, whereas aggregate means
and intervals are computed before that rounding. The audit uses online
sample variance and a rounding-derived tolerance rather than demanding
impossible bitwise equality. For per-score rounding error epsilon=5e-7,
the SEM perturbation is bounded by epsilon/sqrt(n-1); multiply by the
native critical value for the half-width bound. Small additional allowances
cover the aggregate JSON precision and numerical accumulation.

Observed maximum mean discrepancy was 1.264e-8; maximum half-width discrepancy
was 1.360e-10. All assessed counts exactly matched raw rows, ranging from
81,001 to 2,022,806. Every result includes at least 22,398 proteins appearing
in multiple scored pairs. These are descriptive dependency checks, not
independent resampling units.

## Limits

Arithmetic agreement does not validate the underlying ontology, annotation
evidence, frequency estimation, similarity calculation or completeness of
predicted relations. Native intervals do not account for shared proteins or
homologous families and do not provide paired method-comparison uncertainty.
Full-precision raw scores cannot be recovered from the rounded files.
Historical execution environments are not reconstructed by matching their
aggregates to formulas in the retained pinned image. The four recovered
stages remain distinct from historical comparator runs.

Reproduce using the retained raw outputs and pinned image:

```bash
python benchmark_tools/audit_qfo_go_ec.py --output NEW_AUDIT.json
python -m pytest -q tests/unit/test_audit_qfo_go_ec.py
```

Fifteen tests cover both endpoints, direction handling, duplicate and invalid
records, serialization bounds, wrong counts/means and incorrectly interpreting
the native uncertainty field as a single SEM.
