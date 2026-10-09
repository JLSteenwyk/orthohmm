### Retrospective Native-Cell Comparator Sensitivity

Native minus OrthoFinder SwissTrees F1 differences and adjusted intervals
below use raw 0-to-1 units. Family outcomes are descriptive, not independent
pair counts. All planned missing comparisons remain unavailable.

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

All four F1 point estimates are below full OrthoFinder. Only P0/C0/R0
has an adjusted F1 interval excluding zero in that direction. Both
unreconciled cells have lower precision with adjusted intervals below zero.
Reconciled precision differences versus full OrthoFinder have intervals
including zero. Relative to its sequence-only checkpoint, P0/C0/R1 and
P1/C0/R1 have positive precision intervals and negative recall intervals;
their positive F1 point differences have intervals including zero.
P0/C0/R0 also has lower recall than this checkpoint. Eight of 24 estimated
adjusted intervals exclude zero, conditionally within this separate family.
This is not a whole-study significance count, equivalence test or general
superiority finding. No default selection or initial-HMM causal effect follows.

**Native Comparator Figure.** The [six-panel figure](native_qfo_comparator_integration_20261009_v1/native_comparator_sensitivity.pdf)
shows F1, precision and recall against full OrthoFinder (top) and its
sequence-only checkpoint (bottom). Dots: point differences. Thick lines:
nominal 95% intervals. Thin lines: 48-endpoint adjusted intervals.
The figure alone uses percentage points; the table/TSV retain raw units.
All eight cells appear in both rows; unavailable points are not zeros.
[Complete 48-row table](native_qfo_comparator_uncertainty_20261009_v1/TABLE.md),
[full-precision machine result](native_qfo_comparator_uncertainty_20261009_v1/report.json),
[claim-to-evidence addendum](NATIVE_QFO_COMPARATOR_RESULT_20261009.md).

P1/C0/R0 is absent from this supplied fresh-cell snapshot, not globally
missing: the historical high-sensitivity method remains separate.
The other absent cells failed before inference, during assessment or
during native inference, respectively; no score is imputed. The recovered
P0/C0/R1 timing remains failed/ineligible. R changes pair-output semantics
as well as reconciliation; initial HMM search remains on in every cell.
These intervals do not cover other QfO endpoints, the secondary mean,
resource timings or missing outcomes. The terminal failure section below
retains its earlier scope; its no-new-draw statement concerns that failure
reporting update, not this separately specified retrospective calculation.

The [terminal resource consolidation](NATIVE_FACTORIAL_TERMINAL_RESOURCE_RESULT_20261009.md)
joins all 13 observed ablation attempts with explicit accuracy/failure and
resource-scope flags. It does not repair failed measurements or establish
isolated performance. Original shared-host contention limits remain.

