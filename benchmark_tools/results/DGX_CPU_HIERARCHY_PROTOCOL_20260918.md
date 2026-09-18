# DGX CPU Hierarchy Controls

Prospective engineering controls only; do not alter scientific inclusion
thresholds or rerun the27-method timing panel based on these results alone.

Run one exclusive spark-7ff0 allocation,2task CPUs and2GiB,10minute limit,
no requeue. Sequential fixed order: quiet(1second sleep), completed batch
burst(one0.75CPU-second worker), sustained batch load(four sequential
0.75CPU-second workers). A separate one-CPU native step sleeps throughout.
Do not stop unrelated jobs or change service/kernel settings.

Capture host counters before/after each observation point, job-parent CPU
before/after its immediate step-child CPU reads, all immediate step children,
and exact monotonic read brackets. Require stable same-job native/batch
identities, disjoint immediate children, stable hierarchy, nondecreasing
counters, closed read brackets and stable source hashes. Preserve raw values.

For each loaded trial, test whether batch CPU increment is at least0.5seconds
per intended worker and native sleeping-step increment is less than0.25seconds.
These are coarse load-localization checks, not calibrated bounds or a new
scientific screen. Quiet trial is descriptive, with no new pass threshold.
Report every trial and unmet expectation; no selective retry.

Report host busy, each step increment, job outer/inner increments and signed
job-minus-step-sum/host-minus-job arithmetic. Parent and child reads are not
atomic: differences include read padding, direct parent activity and delayed
accounting. They must not be labeled foreign CPU or subtracted from runtimes.
This cannot establish native compute overhead, non-CPU isolation or general
accounting-error rates. Complete-command integration and scientific inclusion
policy remain separate future work.
