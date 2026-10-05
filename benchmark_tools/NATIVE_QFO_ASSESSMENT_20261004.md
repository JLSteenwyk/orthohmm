# Full-Native QfO Assessment

This handoff extends the [full-native pair conversion](NATIVE_QFO_CONVERSION_20261004.md)
without changing any source bound to the live native inference panel. It does
not accept older cached-candidate admissions or relabel their executor revisions.
Seven full-native QfO identities, indices 6 through 12, retain their distinct
`ohmm_qfo_full_native_CELL` participant namespace.

## Admission Before Execution

`run_native_factorial_qfo_assessment.py` requires a successful two-CPU
conversion on bizon, its exact manifest digest, and a bound successful native
terminal review. Recheck the native scheduler, request, plan, input inventory,
independent output/ownership admission, converted counts and protein coverage.
R-on count must equal the independently validated resolved-native pair count;
R-off retains the audited cluster-derived semantics. Empty predictions remain
explicit and are not replaced by fabricated scores.

Recheck the frozen QfO environment manifest and every declared pipeline,
reference, Java, container, executable and Singularity-support record, plus the
native source bindings and current assessment/conversion/admission helpers.
These are file-identity checks, not a new hermetic host-runtime certification.
Conversion and scoring must use the same reference mapping. Source, mapping,
input, output or scheduler changes are refused; no failed native result is scored.

Read-only `--check-only` performs the same admission/path checks without
creating output directories or executing Nextflow. Execution requires its own
eight-CPU allocation. Before an actual launch, use the existing safe-capacity
and scheduler workflow; this adapter does not reserve resources or displace jobs.

## Execution

Use the frozen Nextflow command builder for the 2020 GO, EC, VGNC, SwissTrees,
TreeFam-A and FAS challenges. Never add `-resume`. The new work directories
`qfo_benchmark/w/nq06` through `nq12` retain Darwin's validated path-length
limit. New results are `qfo_benchmark/scoring/full_native_06` through `12`;
execution receipts live in `benchmarks/results/full_native_qfo_assessment_v1/CELL`.
All namespaces are create-once and separate from historical scoring/inference.

```bash
python -B benchmark_tools/run_native_factorial_qfo_assessment.py \
  --root /ABS/REPO --pairs /ABS/CONVERSION/results.json \
  --pairs-sha256 EXACT_CONVERSION_SHA --conversion-job CONVERSION_JOB \
  --check-only
```

Remove `--check-only` only within the scheduled assessment allocation.
Capture immutable preflight/final receipts, command, helper/runtime/reference
pins, scoring log and complete observed output inventory. Postflight rechecks
all records. Failures, interrupts and partial outputs are retained with no
admitted accuracy, silent retry or zero imputation. A zero process exit is
only `process_succeeded_pending_independent_admission`.

Recorded monotonic interval covers endpoint execution and postflight, not
preflight/conversion/inference. It is not CPU accounting, a peak-memory
measurement or an isolated-performance estimate. Shared-host contention has
an unknown, potentially tool-dependent timing effect.

## Independent Endpoint Review

`admit_native_factorial_qfo_assessment.py` rechecks conversion/native/scorer
bindings and actual successful eight-CPU assessment completion. Compare
execution and preflight identities/command, exact output inventory and all
15 fresh completed Nextflow tasks; reject cached or incomplete task traces.
Validate every native metric identity, finite value, standard error and exact
aggregate agreement against the frozen endpoint templates and SwissTrees
family inventory. Require metric outputs to be in the execution inventory.

```bash
python -B benchmark_tools/admit_native_factorial_qfo_assessment.py \
  --root /ABS/REPO --pairs /ABS/CONVERSION/results.json \
  --pairs-sha256 EXACT_CONVERSION_SHA --conversion-job CONVERSION_JOB \
  --assessment-job ASSESSMENT_JOB \
  --output-directory /ABS/REPO/benchmarks/results/full_native_qfo_admission_v1/NEW_NAME
```

The admission directory is separate and create-once. Failed reviews retain a
failure receipt; only fully validated evidence gets
`full_native_factorial_qfo_assessment_admitted`. This does not authorize another
inference identity, establish independent biological validation or complete
the publication goal.

## FAS And Score Semantics

Preserve the pinned native FAS protocol: eligible precomputed relations plus
missing relations with both annotations, not all predictions; newly computed
pairs capped at 9000; unseeded native shuffling and ratio-preserving precomputed
subsampling; uncomputed/failed missing scores omitted. Do not change sampling
or missing-score conventions after observing outcomes.

Independently read the retained raw FAS sample, require canonical unique
accession pairs, verify every sampled pair belongs to the submitted prediction
file, and reproduce the aggregate mean and native pair-IID SEM. Report sample
size, eligible count and sampling fraction separately from submitted pair
volume and input-protein relation coverage. This does not certify the full
eligible-population denominator, annotation completeness, score implementation,
representative sampling, pair independence or paired-method uncertainty.
Those publication-level audits remain required rather than being hidden by
successful endpoint admission.

Only VGNC, SwissTrees and TreeFam-A summaries use native TPR/PPV harmonic
means. GO/EC similarity and FAS are not F1. The six-endpoint mean remains a
project-defined secondary summary. Development exposure, native standard-error
limitations and group-versus-resolved-pair semantics remain explicit.

## Verification Scope

Synthetic tests cover all seven identities including empty predictions,
native/conversion scheduler and source/input/output/mapping refusals, fresh
execution and immutable paths, retained failures/interrupts, and joined native
metric/task/FAS admission. Only the actual scientific runs and their admitted
outputs can supply new benchmark scores; fixture results are not evidence of
method accuracy or speed. Keep initial failing test reports as diagnostics.
