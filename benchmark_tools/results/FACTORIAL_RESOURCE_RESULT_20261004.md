# Retained Factorial Costs Consolidated

## Evidence Recovered

The [16-cell table](factorial_retained_resources_20261004/resources.md),
[JSON](factorial_retained_resources_20261004/resources.json) and
[TSV](factorial_retained_resources_20261004/resources.tsv) now consolidate
the original OrthoBench and corrected-QfO factorial stage measurements.
P denotes profile expansion, C candidate expansion and R reconciliation;
P-off still retains initial HMM search. This supplies missing corrected-QfO
resource reporting without repeating inference or scoring.

All eight candidate-preparation arms have retained elapsed-time observations.
Each arm is shared by its R-off and R-on cells; its repeated table value is
not two independent runs. All eight R-on cells have recorded incremental
wall/user-CPU/system-CPU and sampled process-tree RSS observations. R-off
reconciliation is not applicable, not a zero-cost full pipeline. All 16
full-pipeline entries within these original cached factorial executions remain
unavailable; cached stages cannot be added and called an observed full run.

For corrected QfO, recorded reconciliation wall times are 4,445.429 / 6,861.549 /
4,416.658 / 6,822.365 seconds for P0C0 / P0C1 / P1C0 / P1C1, respectively.
Their sampled summed tree-RSS peaks are 5.925 / 7.332 / 5.921 / 7.308 GiB.
OrthoBench's four retained reconciliation observations span 1,514.539 to
1,963.790 seconds and 1.467 to 1.577 GiB. These are descriptive stage records,
not causal component overhead or isolated speed comparisons.

## Scope And Verification

The collector uses the corrected QfO candidate manifest and four separately
admitted native pair outputs; superseded original-release QfO runs are not
substituted. It checks candidate and membership-constraint bindings, all 15
direct source/admission/metrics file identities, exact cell/factor inventory,
recorded CPU budget, finite costs and the historical RSS convention.
OrthoBench's four original FAILED batch states remain FAILED. Separate
scientific output recovery does not relabel them as clean successful batches.

Actual collection with isolated standard-library Python 3.12.3 (`-I -S -B`)
exits zero. Twenty-three fixture cases pass; adding an actual-artifact readback
brings the focused total to **24 passing**. The independent readback checks
all 16 identities, eight preparation observations and all 32 reconciliation
cost values directly against retained metrics, plus original scheduler states
and deliberately absent full-pipeline fields. The
[execution receipt](factorial_resource_collection_execution_20261004.json)
pins the collector, outputs, tests and both successful JUnit receipts.

The [OrthoBench accuracy/coverage results](ORTHOBENCH_FACTORIAL_RESULTS_20260916.md)
and [corrected-QfO accuracy/prediction-volume table](qfo_corrected_factorial_complete_20260919/scores/scores.md)
remain the companion evidence, using the same eight cell labels. No score,
coverage statistic, method/default, uncertainty endpoint or completed timing
panel changes. The current 34-page manuscript review and rc3 archive remain
unchanged snapshots, not retroactive containers of this new addendum.

## Remaining Cost Requirement

This closes stage-cost consolidation, not all per-configuration cost evidence.
The scoped collector makes no claim that full-inference observations are
absent from every production or scaling record. Those observations remain
separate: assess exact input bytes, frozen settings, runtime and measurement
scope before associating any with a factorial configuration. In particular,
inspect the final 12-proteome resource panel before commissioning new costs;
do not assume its subsets, sequence order or defaults match the factorial.
Only genuinely missing costs justify new runs. Existing successful native,
statistical and scaling measurements must not be repeated merely to fill a
new version or seek faster times.

Historical summed RSS can double-count shared pages and miss short-lived or
between-sample peaks. It is not the newer panel's lifetime cgroup peak; no
later measurement-validity pass or observed contention series is imputed.
Preparation excludes initial search, cache loading and profile/seed building.
Shared-host contention effects are unknown and may differ among configurations;
no corrected isolated estimate or causal timing claim follows. Appropriate
QfO uncertainty, development-family inventory, biological strata and wider
native/archive/deposition requirements remain. No DGX or quiet-window gate is
reintroduced, and no unrelated job/service is modified.
