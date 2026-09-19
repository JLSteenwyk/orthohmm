# Prospective Lineage Read-Crossing Control

## Question And Limits

Can the aggregate lineage reader retain valid scope identity and full-span
CPU evidence when an owned service starts and disappears between its root
counter read and its next ancestor read? This extends the completed
`LINEAGE_LIFECYCLE_RESULT_21989.md` control, which tested disappearance
between snapshots, not during a snapshot. It does not test overhead,
attribution specificity or comparative scientific timing eligibility.

The [upstream Linux v6.11 source](https://github.com/torvalds/linux/blob/v6.11/kernel/cgroup/rstat.c#L524-L589)
uses system-wide CPU categories for root-cgroup usage, whereas non-root
usage uses cgroup runtime accounting. Root usage includes user, nice,
system, IRQ, softirq and steal categories. These are not interchangeable
with an atomic sum of child-cgroup task runtimes. This is upstream source
evidence, not verification of the exact Ubuntu/NVIDIA kernel patch set.
The current DGX reports linux-image-6.11.0-1016-nvidia version
6.11.0-1016.16 (source linux-signed-nvidia-6.11). Current root membership
includes kthreadd and kworker threads. Neither present membership nor the
upstream implementation identifies historical flagged CPU consumption.

## Frozen Workload

Run exactly three sequential controls in an exclusive 20-CPU, 4-GiB DGX
allocation with an exact two-CPU-bound observer step and a ten-minute
outer job limit. Do not run native benchmarking simultaneously. Preserve
all three outcomes; no selective replacement of failed controls.

Reuse the existing finite lifecycle payload and its unique owned service
names, .75-process-CPU-second burn, RuntimeMaxSec=10s, MemoryMax=256M,
TasksMax=8, one allocated CPU affinity, --collect and bounded disappearance
poll. No unrelated service is modified or stopped. This is intentional
outside-observer engineering load, not an undeclared timing contaminant.

For each trial, take a baseline root-to-observer snapshot. Then take a
crossing snapshot with a control-only callback immediately after reading
the root cpu.stat. Within that callback run the existing complete lifecycle
trial, including membership/affinity/CPU checks, disappearance evidence and
its manager/observer snapshots. Record callback start/end monotonic times.
The outer crossing snapshot then reads the remaining ancestor counters.
Take a final root-to-observer snapshot after the crossing snapshot returns.
Do not inject any event into production native measurements.

## Required Evidence

Require unchanged boot, scope device/inode identity, observer membership
and source hashes. Require the successful owned-service CPU and membership
checks from the existing trial, observed service disappearance, and its
manager/full-span outside-target responses of at least .5 CPU seconds.
Require crossing root-read finish <= callback start < callback end <=
next-ancestor read start. Require full baseline-to-final root-minus-target
response of at least .5 CPU seconds. These fixed response thresholds are
engineering checks, not timing-inclusion thresholds or causal bounds.

Retain baseline, crossing and final raw snapshots, callback times, nested
lifecycle evidence, all signed adjacent complements for the full span and
both partial spans, process/kernel/Python identity, source checksums and
Slurm allocation/terminal records. The service may be absent from the root
counter sampled before its creation and present in a later sample; neither
partial-span response is required to exceed a threshold. Do not sum
overlapping ancestor totals or clamp negative complements.

Any exception, missing scope, migration, changed identity, wrong event
ordering or failed service check is a retained failure with partial
evidence. A passing exit alone is not acceptance: independently replay
snapshot validation/comparisons, event ordering and all three outcomes.

This control cannot identify current native residuals, prove that the
monitor detects every transient workload, validate timing overhead with
callbacks enabled, or establish non-CPU isolation. Root/user-slice and
native-only accounting controls and a prospective scientific inclusion
policy remain separate. Historical flags and timing eligibility do not change.
