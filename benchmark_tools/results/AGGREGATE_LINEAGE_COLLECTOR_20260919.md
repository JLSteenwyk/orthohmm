# Aggregate Lineage Collector Development

## Motivation And Implementation

The completed dual-overhead panel 21920 lost tasks 3 and 17 when transient
service cgroups appeared during frontier enumeration. The original failed
observations and all original screening thresholds remain unchanged. No
historical measurement is repaired or admitted by this work.

`probe_cgroup_lineage.py` adds a separate prospective reader of `cpu.stat`
for each scope on the root-to-target path. It does not enumerate siblings.
The [kernel cgroup-v2 documentation](https://cdn.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html#cpu-interface-files)
states that `usage_usec`, `user_usec` and `system_usec` include descendant
processes. Aggregate CPU counters therefore provide a candidate way to
observe activity without requiring each transient service directory to
remain present. This is a design inference, not DGX validation of deleted
descendant accounting or a claim of synchronized measurements.

The reader retains raw counters and per-read monotonic times, checks every
lineage scope's device/inode before and after, and checks boot identity.
It rejects changed, missing, symlinked or aliased lineage scopes. Read failures
raise an exception containing the partial evidence; callers must preserve it.
Interval comparison rejects identity changes, overlapping observation windows
and decreasing cumulative counters. It reports each ancestor delta separately
and signed adjacent differences. Ancestor totals overlap and must not be
summed. Negative differences are retained, not clipped or interpreted as
negative foreign work. Their telescoping sum is an algebraic identity, not
independent verification of resource attribution.

The existing frontier and dual collectors do not import this reader, so no
frozen executor or timing gate has changed. There is no automatic fallback
from a failed frontier observation to a lineage observation.

## Verification

```sh
python -m pytest -q tests/unit/test_probe_cgroup_lineage.py
```

Synthetic tests cover sibling creation/removal during reads, target
replacement, boot changes, malformed counters, read failures with partial
evidence, missing and reordered rows, aliased identities, changed interval
identities, overlapping windows and signed negative complements. Filesystem
fixtures do not model kernel CPU accounting or establish transient-load
detection sensitivity.
All 26 new tests pass; the combined lineage, existing frontier and dual
bracket suite passes all 54 tests.

A one-second read-only smoke on the busy local host successfully sampled
five scopes from root to the current session. One adjacent difference was
negative, as the implementation permits for non-atomic counters. This was
not a controlled experiment or a retained publication timing result. No
native process was started, stopped or migrated by the probe.

## Required Next Validation

Before any replacement overhead or scaling panel, freeze a DGX control
protocol that places a known finite CPU workload in an owned sibling cgroup,
then removes that cgroup before the second observation. Record workload CPU,
membership, parent identities, raw counters and explicit cleanup outcomes.
Also test sibling creation/removal during sampling and intentional lineage
replacement as an invalid observation. Do not modify unrelated services or
infer success from synthetic tests. Use existing Slurm/delegation facilities;
if they cannot provide the required ownership and lifecycle control, report
the missing evidence rather than substituting an unverified control.

After that control, integration still requires native membership validation,
the existing narrow/outer host brackets, pressure and memory accounting,
raw failure retention and independent replay tests. The new reader must not
silently reuse the old frontier report schema. A complete fresh overhead
panel and a prospectively fixed scientific inclusion policy remain necessary.
No threshold is proposed from this local smoke. Root accounting differences,
process migration, accounting delay, overhead and non-CPU interference remain
limitations. This development does not close the controlled-timing requirement.
