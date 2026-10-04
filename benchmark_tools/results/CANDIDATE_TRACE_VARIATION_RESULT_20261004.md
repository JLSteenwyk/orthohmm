# Candidate Cap And Tie Sensitivity

The [retained-trace report](candidate_trace_variation_20261004/diagnostic.json)
and [summary](candidate_trace_variation_20261004/diagnostic.md) compare every
accepted merge in the original OrthoBench satellite arm and three completed
12-proteome native repeats. This extends the previous
[cost/partition linkage](FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md)
with a tested failure mechanism, not another benchmark run or method change.

## Actual Trace Evidence

Each trace contains 8,440 accepted merges. Match records by round and sorted
named source/target memberships, not numeric group labels. The three native
traces share 8,428, 8,439 and 8,426 semantic records with the original. Original-
only records split across rounds 0/1 as 4/8, 0/1 and 4/10; native-only counts
are identical. Two repeats additionally have two original-only and two native-
only anchor identities, consistent with propagation after first-round changes.

Among common accepted records, all source and target cluster indices match.
Forward/reverse hit counts also match exactly. Maximum common support
difference is `2.1316282072803006e-14`; other feature deltas are retained rather
than assumed zero. This is not a comparison of every raw hit or rejected merge.

For common anchors whose selected source memberships differ, there are 2, 1
and 4 anchors respectively. Every one has four accepted attachments in both
traces, reaching the frozen per-round attachment cap. The changed selections
include exactly tied supports and spreads of only a few floating-point units.
The complete named memberships, accepted selections, numerical differences
and direct trace identities are preserved in the report.

Frozen source sorts competing eligible selections by decreasing support,
decreasing margin, source index and target index. The cap makes this ordering
consequential. However, accepted records alone do not reveal all candidates
at the cutoff or prove the origin of the score-bit differences. Unchanged
common cluster indices contradict attributing these observed comparisons
to wholesale relabeling. Do not infer that competing host workloads caused
the assignment differences.

## Controlled Frozen-Engine Experiment

Call the unchanged frozen `merge_supported_satellite_candidate_clusters`
function with the original satellite parameters, including four attachments
per anchor per round and two rounds. Use nine equally supported singleton
satellites plus a ten-gene anchor: 19 genes and 180 directed hits.

| Change | Unattached Satellite | Accepted Merges | Rounds |
| --- | ---: | ---: | ---: |
| Baseline | 8 | 8 | 2 |
| Reverse only the initial cluster order | 0 | 8 | 2 |
| Change only satellite 8's hit scores by one floating-point unit | 7 | 8 | 2 |

The score perturbation is `2.220446049250313e-16`. No evidence gate, source,
gene/species membership or biological label changes. The order-only case
changes no score at all. This demonstrates sufficient mechanisms for different
assignments at equal group/merge counts; it does not causally reconstruct every
historical change or establish an accuracy advantage.

Run the diagnostic with reporting Python 3.12.3 and the exact private Python
3.10.13 interpreter used for native timing. Both use NumPy 2.2.6. All report
fields, including numerical trace projections and fixture outcomes, agree
except the recorded interpreter version. The
[native-interpreter report](candidate_trace_variation_native_runtime_20261004/diagnostic.json)
is retained separately. These are tiny candidate-stage fixtures and trace
readbacks, not full native inference, cross-platform certification or timing.

## Validation And Consequences

Initial combined tests pass 76 cases. Final coverage adds the actual second-
interpreter comparison: 77 cases pass, including full trace projections,
membership validation, infinite-margin handling, direct input/source hashes,
frozen fixture replay and the previous resource-table checks. Both diagnostic
collections exit zero on their first invocation. The
[execution receipt](candidate_trace_variation_execution_20261004.json)
pins source/helpers/tests, both reports and all three JUnit files. The final
portability-guard readback also passes all 77 cases here; new raw-asset tests
skip only when local retained assets are absent, not on checksum mismatch.

The prior structural readback found no changed candidate or final group
intersecting any of the 70 OrthoBench reference families; this is inherited
evidence, not new scoring. Keep all native repeats, difference lists, failures
and shared-host timing limitations. Do not select matching outputs or replace
the frozen method with score rounding chosen after seeing these benchmarks.

A future determinism change needs tests separating reproducible aggregation,
stable semantic tie ordering, cap interactions and legitimately ambiguous
assignments. Stable index tie-breaking alone cannot remove score-bit effects.
If memberships change, freeze a new method and obtain new independent
confirmation. This turn neither implements that new method nor modifies
production defaults, the manuscript PDF or the rc3 archive.

Next audit the genuinely remaining per-configuration cost and broader
publication requirements before authorizing selective new measurements.
Scientific uncertainty, development-family inventory, biological-stratum and
distribution requirements remain active. This mechanism result does not
establish publication readiness or complete the overall goal.
