# Prospective Native Dual-Bracket Diagnostic

Freeze before collecting native outcomes with the new collector. Purpose:
test complete-command paired CPU brackets under actual native workloads,
following the successful single-core controls in job 21911. This is not an
overhead comparison, scientific scaling panel or acceptance-policy revision.

## Fixed Runs

Run three fresh, sequential tasks in this order: OrthoHMM high sensitivity,
OrthoHMM satellite_v2, full OrthoFinder. Use the existing four-proteome
scaling input (73,266 proteins, 36,474,860 amino acids) and frozen native
commands/runtime manifests from pressure panel v2. The command templates
are tasks 1, 3 and 8 of that plan, SHA-256
`b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
Selection is one instance of each prespecified method, not selection by its
previous duration or outcome. Relocate only run-specific output/cache paths
into `dual_bracket_native_v1`; leave dataset and installed runtime paths
unchanged. Retain a fresh pinned plan and transitive execution recipe before
submission. Do not reuse historical native outputs as fresh measurements.

Use spark-7ff0, 20 task CPUs, 96 GiB, exclusive allocation, concurrency one,
900-second native timeout and one-hour scheduler limit per task. No GPU.
Verify idle scheduled state before submission. No SSH inspection of the
active workload; controller-only status monitoring is permitted. Preserve
terminal scheduler records. Do not stop or modify unrelated host services.

## Measurements

Use `measure_native_dual_bracket_step.measure`. It preserves worker startup,
ready/go/completion/release handshakes and records a final sample after
command completion. Collect at one-second intervals, requiring all native
pressure and frontier identity checks. Save every raw point and step-memory
peak, with the usual native exit/timeout record and GNU-time companion.
Preparation and verification remain outside native command boundaries.

Retain original screening unchanged and add narrow per-interval results.
Preserve both whole-command original screening and the separately labeled
dual observation-window summary; these scopes are not interchangeable.
Use the existing CPU thresholds without modification. Never clip signed
residuals, subtract overhead, sum overlapping outer residuals, or attribute
residuals directly to competing work.

## Evaluation

After all three tasks terminate, independently replay complete original and
narrow screens, native pressure and frontier identity checks. Verify the
recipe, runtime, input order/content, native output semantics and complete
native artifact inventory. Compare canonical outputs with the same-method,
same-input previous periodic runs for equivalence, reporting any difference.
Output differences or wrapper failures remain failures, not exclusions.

Report every interval flag and retained pressure total, command duration,
memory and output comparison. A diagnostic is fully observed only if native
execution succeeds and all required measurements validate. For each method,
report whether all narrow intervals pass as a separate descriptive outcome;
do not discard a method or extend a run to dilute a flag. No selective
reruns, alternative species subsets, parameter changes or method-order changes.

Failures or remaining flags require a mechanistic investigation and a new
prospective experiment if necessary. Even three clean runs cannot establish
overhead budgets, repeatability, larger-input behavior or absence of non-CPU
interference. A complete repeated overhead panel and frozen scientific
inclusion policy remain necessary before the 27 scaling runs are admitted.
