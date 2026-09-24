# QfO Parameter Uncertainty: Low CPM Added

The frozen parameter analysis now includes the control, four threshold variants,
and the admitted low-CPM arm. High CPM remains unavailable, not zero. The earlier
threshold-only analysis remains retained as a historical snapshot.

## Method

Unchanged `run_qfo_parameter_uncertainty.py` reconstructed all 18 SwissTrees
family count tables from retained raw predictions and checked their truth and
membership against the frozen reference. It used 100,000 paired family draws,
PCG64 seed 20260925, linear quantiles, and Bonferroni adjustment across all 18
planned endpoints, including the unavailable high-CPM contrasts. Each draw
recomputes native macro precision and recall followed by their harmonic mean;
this is not mean family F1 or independent-pair resampling.

## Low-CPM Results

| Metric | Control | Low CPM | Difference, percentage points | Adjusted interval, percentage points | Family wins/ties/losses |
| --- | ---: | ---: | ---: | --- | --- |
| F1 | 0.833513218 | 0.791013733 | -4.2499 | [-21.9395, 0.7887] | 3/13/2 |
| Precision | 0.955176749 | 0.927710548 | -2.7466 | [-13.7247, 0.4512] | 2/13/3 |
| Recall | 0.739341256 | 0.689427572 | -4.9914 | [-25.3053, 1.2660] | 2/14/2 |

Unadjusted paired F1 interval: [-13.3411, 0.4600] percentage points. All
intervals include zero; neither equivalence nor a general decrease is established.
The four threshold contrasts are unchanged from the earlier analysis.

Per-family F1 differences are positive for BAR (+0.027377), CASP (+0.008040),
and HOX (+0.028665), negative for NOX (-0.035549) and VATB (-0.944437), and zero
for the remaining thirteen families. VATB is a priority for an explicitly
post-hoc pipeline trace, not proof of any search, grouping, or reconciliation
mechanism. Family wins do not weight the magnitude of changes, and their count
does not determine the native aggregate.

## Evidence And Limits

- [Input admission inventory](qfo_parameter_uncertainty_partial_inventory_20260923.json).
- [Reconstructed counts and intervals](qfo_parameter_uncertainty_partial_20260923.json).
- [Independent numerical reproduction](qfo_parameter_uncertainty_partial_reproduction_20260923.json).
- [All low-CPM QfO endpoints](qfo_cpm_low_endpoints_22096_0.json).

The result SHA-256 is
`1d6bd4f68a9d6728e0180bdde31f55ca0da7b4fc22e09bcc797dcb8db6cec628`.
The separate reproduction uses the same random generator and quantile library
but reconstructs the statistic independently from admitted counts. Source and
input identities are checked before and after computation. All 15 available
endpoint contrasts reproduced within an absolute tolerance of 1e-12.

Five of six contrasts are estimable; the full panel is incomplete. This is
development-exposed sensitivity evidence, not independent confirmation and not
selection-adjusted evidence for new defaults. Eighteen families provide limited
resolution and may share evolutionary history. No corresponding family-level
uncertainty is established here for TreeFam-A, the other QfO endpoints, or the
project-defined six-score mean. Shared-host times are not controlled scaling
evidence. No defaults change.
