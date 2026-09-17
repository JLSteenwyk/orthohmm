# Matched Scaling Protocol

Prospective first panel, frozen before new timing outcomes. This addresses the
practical-efficiency requirement, not accuracy tuning or independent validation.
Existing shared-node timings remain descriptive and are not substituted here.

## Inputs and Runs

Use the twelve checksum-pinned complete OrthoBench proteomes from
`orthobench_factorial_prepared_20260916.json` (SHA256
5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382).
Order basenames by SHA256 of `orthohmm-publication-scaling-20260916-v1`, a newline,
and the basename. Use nested prefixes of 4, 8 and 12 proteomes, without removing
proteins or rewriting FASTA content. Record protein and sequence-character counts,
membership and file hashes. Ordering uses neither accuracy nor runtime outcomes.

Run frozen OrthoHMM high-sensitivity, frozen satellite_v2 with inferred phylogeny,
and OrthoFinder 3.1.5 full inference three times each at every size: 27 planned
runs. Rotate the three-method order by repeat plus size index. A run requires a
fresh output directory; do not share inferred search results between methods or
repeats. Record OS page-cache state as uncontrolled; do not flush global caches
on this shared machine. Method settings and native dependency hashes must match
the admitted baseline or be explicitly recorded as a different method version.

OrthoFinder's pre-phylogenetic checkpoint is not a separate complete run. A
checkpoint duration may be reported only if its timestamp and included stages
are observed directly; never subtract unrelated historical timings. Other tools'
historical resource records remain visible with their comparability limitations.

## Resource and Execution Gates

Use identical 32-CPU, 128-GiB limits and a 24-hour per-run wall limit on the same
machine, with no overlapping panel runs. Check requested and actual affinity,
thread settings, CPU model and native runtime. Schedule an exclusive allocation
where available, but also inspect non-Slurm workloads: scheduler exclusivity alone
does not establish a quiet host. Do not terminate unrelated user processes.

Before execution, freeze exact CLI commands, stage boundaries and a tested
resource collector. Record wall time, total CPU time and simultaneous process-tree
memory, with sampling interval and missed-short-process limitations. Retain GNU
time and scheduler records as additional evidence; neither an individual-process
RSS peak nor tracemalloc is the simultaneous total of all native child processes.
If cgroup accounting is available, document its scope and distinguish page cache
from process RSS. Monitor unrelated host activity throughout each run. Classify
contaminated runs explicitly and retain their measurements; replacements, if
necessary, need a documented rule before inspecting comparative speed.

Separate input preparation, tool inference, output conversion and scoring.
End-to-end inference starts before tool database construction and ends only after
native final output completion. Cached/reused stages must not be mislabeled as
end-to-end. Verify complete valid output coverage separately from process exit.
Record failures and timeouts with consumed resources; do not treat them as fast
successful runs or silently omit them from the planned inventory.

## Reporting and Limits

Show every run, median and range for the three repeats, input sizes, failures and
workload flags. Report wall time and memory separately; do not select each tool's
fastest repeat or fit a universal complexity law to three dataset sizes. Paired
comparisons require comparable completed runs, with exclusions made explicit.

This is one nested series: taxon identity and input size co-vary, and full
proteomes preserve their actual duplication and sequence-length distributions.
Results cannot alone establish general scalability or extrapolation to QfO-sized
inputs. Broader dataset coverage and existing resource limitations remain part of
the publication completion audit. Prepared inputs or passing collector tests are
not evidence that the scaling experiments have completed.

## Accounting Feasibility Check (17 September)

The read-only `slurm_resource_snapshot.py` now verifies a process's requested
Slurm job/step subtree, reads cgroup-v2 CPU/memory counters and inherited job
memory limits, and samples RSS for its descendant processes. It neither changes
limits nor resets peaks nor moves a process out of Slurm. The
[QfO live snapshot](qfo_checked_full_resource_snapshot_20260917.json) confirms
these interfaces are readable on this machine; this is not a scaling result.

The kernel's `cpu.stat` counters include descendants. Its memory accounting
includes charged file-cache and kernel allocations; `memory.peak` reflects the
peak since creation or reset, not necessarily the inference interval. These
quantities must remain distinct from sampled summed process RSS, which is
non-atomic and can count shared pages repeatedly.
[Linux cgroup-v2 documentation](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html).

A bounded time-series collector is now implemented in `monitor_slurm_resources.py`.
Its [five-sample live smoke test](qfo_resource_series_smoke_20260917.json) retained
[raw observations](qfo_resource_series_samples_20260917.jsonl), verified stable
anchor identity/scope and nondecreasing cumulative CPU counters, and measured an
8.003-second observation span. The CPU difference was8,000,914 microseconds;
the maximum sampled summed RSS was6,861,684,736 bytes. These measurements only
describe that window, which was not isolated from other workloads or unit tests.
Sampling failures are retained and do not imply that the Slurm job terminated.

`measure_slurm_command.py` now wraps command launch through exit in a dedicated
Slurm task subtree, retaining zero/nonzero exit, timeout, spawn failure and
measurement failure separately. The [Slurm lifecycle smoke](resource_command_smoke_20260917.json)
completed as job21330 with two allocated CPUs and a one-GiB inherited limit.
Thirteen observations captured wrapper-only before/after and wrapper plus parent
and child during the command. Real subprocess tests also cover failures and
timeout cleanup, including a TERM-resistant child. The smoke was not a scientific
scaling run or an isolated performance comparison. Native output correctness is
not inferred from exit zero, and a sampling failure does not stop native execution
before its own exit or declared timeout.

Scientific inference/stage command boundaries, unrelated-workload attribution,
whole-run workload checks and native-output validation are still required before the
27 timing runs. An anchor's task subtree does not necessarily contain every
Slurm step in a job; record the scope rather than labeling it whole-job memory.

## Host Competition Probe (17 September)

`observe_host_competition.py` supplies read-only CPU competition evidence from
two process snapshots. It matches PID plus creation time, verifies cgroup
membership, excludes the measured task subtree and monitor, and retains sampling
errors, unmatched processes and moved memberships. A fixed 0.25 average-core
threshold detects persistent foreign CPU work; a negative result is not proof of
host exclusivity. No unrelated processes are signaled or modified.

A three-second probe anchored to QfO job 21333 observed approximately 40.406
foreign average cores, one sampling error and two unmatched foreign processes.
These observations include any work outside the measured task subtree, not only
other projects. The raw local report is
`benchmarks/results/host_competition_qfo_21333_v1.json`, SHA256
5cc9ea4a693735e57fc82bbb5a69d9b1e9a35c9c8baa820068c2b64c2216ce3b.
The whole-host process inventory is not copied into the publication repository.
Observer source SHA256:
8b0f20b50aaa0efee5fb3f7a4b924323d771e39b791be5aaac15caac6c132113.

Twelve focused tests cover identity, scope, counter and uncertainty handling.
The probe does not observe short-lived work between snapshots, I/O, GPU or
memory-bandwidth contention. Whole-run integration and a controlled execution
window remain required. None of the 27 scientific scaling runs has begun.

## Command-Lifetime Workload Collection

The command wrapper now supports `--monitor-host`, using
`command_host_monitor.py` to stream local raw process snapshots before command
launch, during timeout-controlled waits, and after command exit. Only the last
host snapshot is held for interval comparisons. The report records whether
samples bracket the entire command, interval classifications, observation errors,
maximum observed competing CPU use and the largest between-snapshot gap.

Missing workload observations never become quiet evidence. Detected competition
takes precedence in classification; otherwise sampling errors, process churn or
incomplete coverage yield an inconclusive result. Workload errors do not change
the native command exit status or stop its execution. The workload result is
separate from the cgroup resource-measurement result. Every classification retains
`controlled_workload_verified: false`: visibility between samples, storage/GPU
contention and host exclusivity remain unproven. Collection overhead is included
in measured resources. Raw whole-host process inventories remain local rather
than entering published result tables.

Prepared a dedicated two-CPU, one-GiB parent/child Slurm smoke test of the complete
collector lifecycle. This is not one of the 27 scientific scaling runs. A
controlled execution window, frozen scientific commands and native-output
admission remain required before those comparisons.
