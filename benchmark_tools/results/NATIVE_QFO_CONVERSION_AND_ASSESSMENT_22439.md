# Full-Native QfO Conversion Completed And Assessment Started

Resume from pushed `4c561b6d` and verify original native `22437` remains
live. Completed native `22435`, index6/P0C0R0, has a successful original
terminal review. Do not rerun inference or reuse older cached QfO scores.
Initial sensitive HMM search remains on; downstream profiles, candidates
and reconciliation remain off for this configuration.

## Retained Pre-Execution Correction

The [original held submission22438](native_qfo_conversion_submission_22438.json)
incorrectly used the path returned by `record()` for the Python invocation.
That helper resolves symlinks: the prepared review environment's interpreter
has the same binary as the native base environment, but different installed
packages. An actual sanitized-environment import probe fails at `Bio` with
the base path and succeeds with the prepared environment path. Equal binary
checksums do not establish equal Python environments.

Cancel this owned held job before release. Retain the
[cancellation and environment probes](native_qfo_conversion_cancelled_22438.json)
and [actual original batch](convert06_22438.sh). Scheduler accounting records
CANCELLED, zero elapsed, no start and no steps; at cancellation no output
directory or execution logs existed. This is a retained erroneous submission, not an
executed converter failure, successful run or native inference retry.

The first cancellation observer incorrectly used the native-specific
`verify_terminal()` on a conversion allocation. Its identity check raises
`ValueError` after cancellation had executed. The retained readback records
that error and checks the actual converter controller/accounting directly;
no second cancellation or duplicate submission is hidden.

## Corrected Conversion

[Submit22439 explicitly](native_qfo_conversion_submission_22439.json),
preserving the prepared review environment's invocation path separately from
its binary identity. Actual environment import succeeds: Python3.10.13,
Biopython1.87, NumPy2.2.6 and psutil7.2.2. The converter source is unchanged:
SHA `472f26a140cf89aed36f8212a31fef7d8f5061204dac2e0ce80c3c7e6746630e`.

The [independent held gate](native_qfo_conversion_pre_release_22439.json)
checks owner, exact two-CPU/32-GiB/six-hour/no-requeue envelope, request
Comment, generated batch/parsed argv, original request and successful
native review, source/interpreter/environment, fresh destination, all920
helper pins and safe RAM683,144,740,864bytes. Both shell syntax checks pass.
Retain [one actual release](native_qfo_conversion_released_22439.json) and
[scheduler-generated batch](convert06_22439.sh).

Job22439 and its batch complete0:0 in1:16. The
[byte-identical conversion result](native_qfo_conversion_result_22439.json)
retains all source/output bindings and reports:

| Quantity | Value |
| --- | ---: |
| Cross-species group-derived clique pairs | 9,009,082 |
| Pairs lost to reference mapping | 0 |
| All inference input accessions | 984,137 |
| Accessions present in at least one pair | 546,450 |
| Fraction of all inputs present in a pair | 0.555258058583307 |

These are R-off group cliques, not phylogenetically resolved pair predictions.
Coverage is not precision, recall, F1 or reference-family coverage. Singleton
and no-relation proteins stay in the all-input denominator. Original and
filtered pair files are byte-identical,138,866,068bytes each. Raw pairs stay
outside Git.

The [independent readback](native_qfo_conversion_readback_22439.json)
streams all pairs, checks canonical row format/count and unique relation
accessions, validates the existing assessment preparation/binding gates,
and retains the successful converter scheduler. It is not an independent
reimplementation of clique conversion or another raw-resource replay.
An exact-process scan occurs after this short job has already completed;
zero candidates are not a job failure. No converter startup PID is claimed.

## Six-Endpoint Assessment Now Running

Fresh scoring namespaces are `full_native_qfo_assessment_v1/p0_c0_r0`,
`qfo_benchmark/w/nq06` and `qfo_benchmark/scoring/full_native_06`.
[Held assessment22441](native_qfo_assessment_submission_22441.json) uses
unchanged `run_native_factorial_qfo_assessment.py`, eight CPUs/64GiB/24hours,
no requeue and the exact prepared Python invocation. The
[independent held gate](native_qfo_assessment_pre_release_22441.json)
passes owner/envelope/Comment/batch/fresh paths/pins and safe RAM;
[release](native_qfo_assessment_released_22441.json) occurs once.
Retain the [actual batch](assess06_22441.sh).

[Actual startup](native_qfo_assessment_start_22441.json) verifies driver
PID283857, exact argv, matching job cgroup/eight-CPU environment and a64-GiB
ancestor memory cap. The job-tagged preflight binds the original prepared
command, conversion and native index. Nextflow has submitted validation,
conversion and metric scheduling tasks. Startup is not final success,
endpoint admission, continuous integrity or hermetic runtime certification.

Evaluate VGNC/SwissTrees/TreeFam-A F1, GO/EC similarity and FAS separately;
the six-metric mean is only a project-defined secondary summary. Preserve
the frozen FAS eligible-population, unseeded sampling and missing-score
limitations; native pair-IID SEM is not a paired-family confidence interval.
After whole assessment termination, inspect actual execution, all native
tasks/endpoints and FAS sample/provenance with the existing independent
admission CLI. No score or admission is currently claimed.

The [milestone readback](native_qfo_conversion_assessment_milestone_20261005.json)
checks all14 retained JSON receipts, four byte-identical published batches,
all920 helper pins and the actual original native PID/scheduler handles.
Its first observer passed a path to a record-only reader and failed before
writing a report; the corrected observer retains that scope error. No
scientific output, job or helper source changed to make the readback pass.

Conversion and scoring are separate from matched inference timing and may
overlap original22437; its host monitor records competing workloads. CPU,
memory-bandwidth and I/O contention have unknown, potentially tool-dependent
effects. These are shared-host observations, not isolated performance or
corrected speed estimates. No unrelated job/service is modified. No code,
defaults or frozen scientific/helper settings change, and no new test-suite
result is claimed. The full publication goal remains active and unproven.
