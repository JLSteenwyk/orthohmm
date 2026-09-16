# Corrected Variable-Length Simulation Results

[Results table](SIMULATION_VARIABLE_NATIVE_RESULTS_20260916.md) is generated
from `simulation_variable_native_results_20260916.json`. Pinned assembler
b66225d separately verified corrected OrthoHMM array 21142 and original
OrthoFinder array 21010, using their own executors and provenance. These are
the first inspected accuracy results for this frozen heterogeneous-length
panel. Scientific settings were not revised after viewing them.

## Accuracy And Uncertainty

Full OrthoFinder has higher paired mean F1 than each OrthoHMM mode in all seven
conditions. All seven high-sensitivity differences have Bonferroni-14 intervals
below zero. For satellite_v2, four adjusted intervals are below zero; turnover,
missing20 and uneven_taxa intervals include zero. These approximate bootstrap
intervals are conditional on successful paired seeds, not failure-adjusted
population comparisons or proof of universal superiority.

Satellite_v2 reduces the observed gap relative to high sensitivity. Baseline
mean F1 is 99.46% versus 99.94% for full OrthoFinder; turnover is 98.84% versus
99.36%. The largest deficits occur under divergence: satellite_v2's paired
F1 differences are -11.74 points for divergent (five paired seeds) and -12.28
for divergent_turnover (eight paired seeds). The available-case table means
use different seed sets when methods fail and must not be subtracted to
reconstruct those paired estimates. High precision with lower recall is the
dominant observed OrthoHMM tradeoff here; causal stage attribution needs the
planned error tracing, not speculation from aggregate scores.

## Failures And Scope

- High sensitivity passes native admission for all 70 datasets; satellite_v2
  for 67; full OrthoFinder and its parent-gated sequence checkpoint for 65 each.
- Satellite_v2 fails species-tree inference for divergent seed 20261109 and
  divergent_turnover seeds 20261109/20261110. Native logs report no connected
  single-copy family coverage for n11/n12 (seed 9) or n2 (seed 10).
- Full OrthoFinder has nonfinite graph weights in divergent seeds 20261101,
  20261102, 20261107, 20261108 and 20261109. Both full and checkpoint results
  are excluded for those runs. The exact cause of these remaining numerical
  failures requires further diagnosis; the fixed-length explanation must not
  be assumed to account for every heterogeneous-length failure.
- All 140 OrthoHMM records show profile construction (62-193 profiles), but
  zero added profile edges. This panel therefore does not demonstrate an
  accuracy contribution from multi-sequence profile expansion. It retains the
  HMM-based initial search and is not an HMM-free control.

Keep this panel separate from the fixed-length stress panel. It has synthetic,
family-specific lengths without within-family indels or realistic domain
architecture. It does not establish generalization to arbitrary proteomes.
Record failures, limited seed counts and shared-machine timing caveats. The
results do not support an overall OrthoHMM-superiority claim. Independent
curated validation, ablations and the remaining publication work are required.
