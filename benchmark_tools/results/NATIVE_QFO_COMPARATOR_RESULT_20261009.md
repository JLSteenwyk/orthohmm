# Native-Cell SwissTrees Comparator Sensitivity

The [frozen retrospective protocol](NATIVE_QFO_COMPARATOR_UNCERTAINTY_PROTOCOL_20261009.md)
and tested consumer were committed before the calculation. Production executed
once at source revision `fc9ff1282be85d5b1d6632392fc0c8cee4c76c61` with exit 0.
The [machine-readable result](native_qfo_comparator_uncertainty_20261009_v1/report.json)
is 50,060 bytes, SHA256
`aea49a83bcce9a80237fd1c8f04ada43efeb9bd41f6e17ba0aa0c7f911018c02`.
The [execution and independent readback](native_qfo_comparator_execution_20261009_v1.json)
retains both actual commands, executed code, outputs and zero exits. No native
inference, conversion, scoring, timing or failed-attempt retry occurred.

## Methods

Four admitted native cells and both retained OrthoFinder 3.1.5 modes share 18
ordered SwissTrees families, 563 disjoint represented genes and 10,765 reference
relations. The official count conversion halves raw relation counts before a
unit prior: PPV is `(TP+2)/(TP+FP+4)` and TPR is `(TP+2)/(TP+FN+4)` in terms of
the stored raw counts. F1 is the harmonic mean of macro-family PPV and TPR,
not mean family F1. Raw serialized native endpoints remain distinct from these
full-precision reconstructed statistics.

The fixed calculation uses 100,000 shared multinomial family draws with PCG64
seed 20260920. Each replicate recomputes the actual macro statistic. Candidate
minus comparator differences use raw 0-to-1 units; nominal percentile intervals
use quantiles 0.025/0.975, and Bonferroni intervals use 0.05/96 and 1-0.05/96
with linear interpolation. All 48 planned endpoints remain in the correction:
eight cells, two comparators and three metrics. There are 24 estimated and 24
unavailable endpoints, not a correction restricted to successful cells.

Independent readback imported no production helpers. It reconstructed every
family statistic from counts, checked member lists/truth totals and macro
arithmetic, replayed the same deterministic weights, verified all 24 nominal
and adjusted intervals and family differences/outcomes, and checked every TSV
and Markdown row. This replay is numerical verification, not another dataset
or new independent scientific confirmation. Seventeen source/helper/protocol/
input pins were checked before and after readback. The final focused/adjacent
source suite passed 72 tests; this does not itself prove bootstrap coverage.

## F1 Results

All differences and intervals below are raw 0-to-1 values. Full-precision
precision/recall differences, nominal intervals and all missing rows are in the
[48-row TSV](native_qfo_comparator_uncertainty_20261009_v1/intervals.tsv) and
[complete table](native_qfo_comparator_uncertainty_20261009_v1/TABLE.md).

| Native Cell | Comparator | F1 Difference | 48-Endpoint Adjusted Interval | Family Wins/Ties/Losses |
| --- | --- | ---: | --- | --- |
| P0/C0/R0 | Full OrthoFinder | -0.159229 | [-0.272661, -0.007023] | 1/0/17 |
| P0/C0/R1 | Full OrthoFinder | -0.058839 | [-0.164040, 0.018543] | 3/4/11 |
| P0/C1/R0 | Full OrthoFinder | -0.162704 | [-0.298476, 0.013904] | 1/0/17 |
| P1/C0/R1 | Full OrthoFinder | -0.062041 | [-0.168638, 0.018129] | 3/4/11 |
| P0/C0/R0 | Sequence-Only OrthoFinder | -0.001331 | [-0.094152, 0.098224] | 7/5/6 |
| P0/C0/R1 | Sequence-Only OrthoFinder | +0.099060 | [-0.117354, 0.302180] | 12/0/6 |
| P0/C1/R0 | Sequence-Only OrthoFinder | -0.004805 | [-0.081246, 0.067536] | 5/8/5 |
| P1/C0/R1 | Sequence-Only OrthoFinder | +0.095857 | [-0.120221, 0.300133] | 12/0/6 |

All four F1 point estimates are below full OrthoFinder. Only the P0/C0/R0 F1
adjusted interval excludes zero in that direction; the other intervals include
zero. Both unreconciled cells have lower precision than full OrthoFinder with
adjusted intervals below zero. Reconciled precision differences relative to
full OrthoFinder are small and their intervals include zero.

Relative to the sequence-only checkpoint, P0/C0/R1 and P1/C0/R1 each have a
positive precision difference with adjusted interval above zero and a negative
recall difference with adjusted interval below zero. Their F1 point differences
are positive, but their F1 adjusted intervals include zero. The P0/C0/R0 recall
interval is also below zero. Overall, eight of 24 estimated adjusted intervals
exclude zero. This is a conditional descriptive accounting, not a global
significance count or whole-study correction.

## Claim-To-Evidence Limits

| Claim | Assessment And Supporting Evidence |
| --- | --- |
| Native reconciled cells outperform full OrthoFinder on SwissTrees F1 | Not established; both differences are negative and adjusted intervals include zero. See the F1 rows and complete result above |
| Native reconciled cells show a precision-recall trade-off versus the sequence-only checkpoint | Supported conditionally in these 18 exposed families: precision intervals are positive and recall intervals negative in both cells. All values and family outcomes are retained in the TSV |
| An interval including zero proves equivalence | Not supported; no equivalence margin or equivalence test was specified |
| These intervals validate a new default or initial-HMM benefit | Not supported; all cells retain initial HMM search, and P toggles downstream profile refinement. The protocol forbids tuning or promotion based on these exposed results |
| These intervals cover other QfO endpoints, the six-metric mean or resource rankings | Not supported; they apply only to F1/PPV/TPR for SwissTrees and do not transfer to VGNC, TreeFam-A, GO, EC, FAS or timings |

Sequence-only OrthoFinder is the retained MCL-checkpoint diagnostic, not a
separately timed full method. R changes group-clique to resolved native-pair
semantics as well as reconciliation. P0/C0/R1 has recovered admitted accuracy
but failed/ineligible timing; no timing was recovered by this analysis.
P1/C0/R0 is absent from the supplied fresh snapshot, not globally missing:
its historical high-sensitivity method remains separate. P0/C1/R1 failed before
native inference, P1/C1/R0 failed assessment with OOM, and P1/C1/R1 failed native
inference with SIGSEGV. None supplies an imputed comparison.

This is retrospective, development-exposed and conditional, not selection-
adjusted inference or independent generalization. Eighteen disjoint reference
families do not prove exchangeability or remove merged-prediction correlations.
Finite bootstrap-tail resolution and approximate percentile coverage remain
limits. Earlier 24-endpoint comparator and 42-endpoint factorial analyses are
unchanged. New results do not narrow their original scope or constitute a
whole-study multiplicity correction.

## Integration Status

This closes the executable native-cell versus OrthoFinder SwissTrees interval
calculation/readback, not the seven-part publication goal. The older manuscript,
figures, rc5/rc6 archives and frozen requirement reconciliation remain unchanged.
This result is an explicit scientific addendum pending integration into a
separately generated manuscript and comparison figure. The
[terminal resource consolidation](NATIVE_FACTORIAL_TERMINAL_RESOURCE_RESULT_20261009.md)
also closes the older reconciliation's stale resource-join next action; it does
not supply isolated performance estimates. Other endpoint uncertainty,
biological-descriptor truth and historical provenance limitations remain.
