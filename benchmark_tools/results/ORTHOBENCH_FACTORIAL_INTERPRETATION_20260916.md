# OrthoBench Component Results

[Generated results](ORTHOBENCH_FACTORIAL_RESULTS_20260916.md) contain all eight
prespecified cells, 70 RefOGs, 20,000 paired draws (seed 20260918) and all
36 multiplicity-adjusted contrast/metric endpoints. Source JSON is
`orthobench_factorial_results_20260916.json`, SHA256
`6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141`.
P denotes profile expansion, C candidate expansion and R reconciliation.
This remains development-exposed evidence, not independent confirmation.

## Findings

- Reconciliation raises F1 in every matched setting: +2.942/+2.952 percentage
  points without candidate expansion and +6.763/+6.899 with it. All four
  adjusted F1 intervals are above zero. Precision increases in all four;
  recall decreases, with adjusted recall intervals below zero for expanded
  candidates. R also changes output semantics from candidate co-membership
  to root HOGs and includes the configured membership filter where applicable.
- Candidate expansion alone reduces observed F1 by about 3.1 points. With
  reconciliation it raises observed F1 by 0.697/0.795 points, but all four
  adjusted candidate-expansion F1 intervals include zero. Recall increases
  and precision decreases in all four settings, with adjusted intervals
  excluding zero. Larger candidates are not unconditionally better.
- Profile expansion raises observed F1 by 0.568-0.704 points across the four
  settings. All adjusted F1/P/R intervals include zero. Profile-off still
  retains the initial HMM search; this is not a total HMM-versus-sequence
  comparison. Most RefOGs tie on family F1 (59-61 of 70).
- The candidate-expansion/reconciliation interaction is descriptively
  +3.821/+3.946 F1 points across the two profile settings. No inferential
  interaction test was prespecified; these numbers are not significance claims.
- The full P1/C1/R1 configuration achieves 74.106074% F1, reproducing the
  historical satellite_v2 point estimate. This factorial contains no new
  OrthoFinder contrast and does not establish superiority over OrthoFinder.

## Validation And Costs

All eight F1/P/R scores match the official OrthoBench CLI to its one-decimal
printed precision; exact-family counts also match. All cells partition the
same 251,378 input genes, retaining singletons. Multispecies-group membership
and group counts are reported separately from accuracy.

All four native processes succeeded but their original batches failed the
known cwd-dependent package-inventory postflight check. Independent recovery
rehashes all outputs and verifies frozen sources, inputs, commands, tools,
native summaries, inferred trees and membership policies. Original FAILED
statuses remain in the JSON; the table's "Failed Cells: None" means no
scientific cell was excluded after recovery, not failure-free batch execution.

Reconciliation wall times are 1,515-1,964 seconds with sampled summed
process-tree RSS peaks of 1.467-1.577 GiB on a shared node. These are cached
incremental costs, not end-to-end runtimes or matched scalability evidence.
Candidate-only cells have no independently measured per-cell runtime.

QfO factorial evaluation, an unconstrained membership-filter diagnostic, a
matched sequence-search control and independent curated validation remain
required. Do not tune the frozen YGOB method in response to these outcomes.
