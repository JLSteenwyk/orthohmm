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

The [lifecycle smoke report](resource_host_command_smoke_20260917.json) completed
as job21458, exit0, scheduler elapsed11seconds, using frozen9629903 executor.
Four host observations bracketed the native command and allthree intervals
detected competing CPU work (maximum43.690 observed foreign average cores).
No workload observation exceptions occurred. The native command exited0, but
neither this smoke nor its resource measurements establish isolated efficiency.

Importantly, each whole-host snapshot took2.081-2.177seconds; total observed
cgroup CPU was8.937seconds over9.667seconds. Thus a requested0.5-second wait is
not a0.5-second sampling cadence, and the monitor itself materially burdens this
short smoke. Before scientific timing, reduce or amortize collection cost with
a prespecified cadence and retain the resulting missed-work limitations. Do not
treat functional monitoring success as acceptable measurement overhead.

## Lower-Cost Collection and Separate Host Cadence

Prospective update before any scientific scaling run: the observer uses
psutil's public `oneshot()` context for shared process reads and its fresh
`is_running()` PID-reuse check outside that context, followed by an uncached
cgroup-membership check. No private psutil API or custom proc-stat parser is
used. The installed7.2.2 implementation compares fresh monotonic process-start
identity on Linux; no per-process second epoch conversion is needed.

Three alternating old/new live snapshot pairs observed1,946-1,948 processes,
zero sampling errors, and consistent observer PID/start identity. Old elapsed
times were1.769263,1.768316,1.767543seconds; new times were1.060566,1.055259,
1.059101seconds. This is a local collector diagnostic on a busy host, not a
controlled inference speed comparison. Old source SHA256:
8b0f20b50aaa0efee5fb3f7a4b924323d771e39b791be5aaac15caac6c132113;
new source SHA256:
d8b66001bdabf4c436aab8185f218dac019dd45b3075c71521513235a1c07e0a.

The wrapper adds `--host-interval` (default30seconds) independently of the
cgroup `--interval`. It always attempts pre-launch/post-exit scans; intervening
scans become due30seconds after the preceding host scan finishes and occur on
the next command-wait timeout. Thus30seconds is not an exact start-to-start
cadence, and actual gaps remain reported. Short-lived competition can be
missed, and no quiet-host certification is introduced. Observer overhead is
still part of cgroup accounting. At this host size, one approximately1-second
scan per30seconds remains nontrivial; acceptable overhead for scientific
measurements still needs explicit evaluation.

A prepared35-second parent/child Slurm smoke uses2CPU/1GiB,1-second cgroup
waits,30-second host intervals and a90-second command timeout. It exercises
the independent cadence rather than changing or repeating inference. Frozen
execution and recorded output verification precede any use in scaling runs.

The [cadence smoke report](resource_host_cadence_smoke_20260917.json) completed
as21622, exit0, scheduler elapsed39seconds, frozen executor65b8b04. All linked
source/evidence hashes were reread and verified. Three successful host scans
bracketed the35.2966-second native command;35 resource observations spanned
37.9928seconds. Largest gap between host scans was30.1346seconds. Host scans
took1.2480,1.3390 and1.3719seconds. Cgroup CPU use was4.5153seconds, including
wrapper/observer cost; this must not be labeled inference-only CPU cost.

Both host intervals detected competing CPU work, with maximum46.1897 observed
foreign average cores, zero observation exceptions and no resource sampling
errors. The collector operates at its prescribed separate cadence, but this
does not establish acceptable overhead for every command or a controlled timing
window. None of the27 scientific scaling runs has started. Raw process
inventories remain local. Snapshot SHA256:
130d2d4ea52f3731baa6ac4b835251c03d6533a3e7ca8ea642dd82ae89eb6e47.

## Native Command Boundaries Prepared

The [27-command manifest](publication_scaling_commands_20260917.json) is now
prepared, SHA256
25345e7a5d49e7474b09188498dc4760b298512266cd21510fba56ad6401e53d.
Its source, frozen inputs, baseline environment, native profile smoke,
tool-resolution and source hashes passed preflight. All27 unique inference
destinations are still absent; no scaling inference has run.

The baseline OrthoHMM harness executes native inference and then hashes inputs,
outputs and source files and counts output lines. Timing the entire harness
would mix inference with external reporting. The manifest therefore records
both the reference harness configuration and its equivalent direct native CLI
command. Tests intercept the actual harness subprocess launch and require exact
argument equality for high-sensitivity and satellite_v2. The tested harness
bytes match the frozen core harness SHA256
47afa3439fea5e0c6c9ac9d62a1c1f819394ba56ec65a55f7192ac9407f918dc.

Only total CPU allocation changes from the admitted simulation configuration:
OrthoHMM `-c32` retains `--threads_per_worker4`, and OrthoFinder uses `-t32 -a32`.
OrthoHMM uses the frozen core working directory/PYTHONPATH, native metrics output,
and unchanged scientific settings, including inferred satellite_v2 phylogeny.
The timer must encompass native CLI startup through exit and native file writing,
while input copying, source/output hashing, conversion and scoring remain
separate. OrthoFinder gets fresh per-run input copies and no restart flags.

Independent native validation remains required: the existing simulation
OrthoHMM validator expects harness-added provenance, while its OrthoFinder
input-copy check assumes `.fasta` basenames. Scaling retains real `.fa` inputs
and bypasses the external harness. Do not waive these checks, rename inputs
silently, or alter old frozen validators. Prepare a scaling-specific validator
for direct native outputs and caller-recorded before/after provenance before
authorizing runs. Controlled scheduling and overhead assessment also remain
outstanding. This manifest is preparation evidence, not timing results.
