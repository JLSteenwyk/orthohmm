# Remaining Native CPU Flags

## Scope

This is a post-outcome descriptive join of the already audited three-run
diagnostic, not a fresh raw-counter replay or causal analysis. All 1,992
intervals are retained, including all 22 narrow-screen flags. The complete
audit is SHA-256-pinned; each retained measurement report and OrthoHMM metrics
file is checked against that audit's inventory before and after reading.
No collector, threshold, output or scientific eligibility changed.

[Machine-readable description](dual_native_flag_description_21912_20260919.json),
SHA-256 `1df35c1c1b97260f83ed04cb30c5d8eedc5d4b1448d930a1c0db2ac312a7d391`.
The report contains every interval, all/flagged/unflagged summaries and
within-stage summaries, not just selected examples.

## Phase Localization

The frozen core at `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` establishes
serial stage order. The recorded stage durations lack start timestamps, so
all uninstrumented native-command time is allowed to occur before any stage.
Bounds also allow microsecond duration-rounding uncertainty. An interval is
labeled only when it lies inside a stage for every permitted gap placement.
Boundary intervals remain unresolved; no OrthoFinder stage is assigned.

All 21 satellite_v2 flags lie within the phylogeny stage under these bounds:
two around562-564seconds after command start, the others around703-744seconds.
There are 255 certainly phylogenetic intervals: 21 flagged and 234 unflagged.
The possible phylogeny start lies within559.945946-560.682319seconds and its
possible end within816.905326-817.641699seconds. No search/profile-expansion
interval is flagged in this run. High sensitivity has no narrow flags.

## Within-Phylogeny Comparison

Entries are medians of retained intervals, not independent samples or paired
causal effects. Host creation/context-switch counts use the host bracket;
cgroup and pressure observations use different windows.

| Recorded quantity | Flagged21 | Unflagged234 |
| --- | ---: | ---: |
| Host process creations | 5,695 | 1,520.5 |
| Host context switches | 35,840 | 14,992 |
| Native average CPU cores | 19.3954 | 19.4897 |
| Native CPU PSI some, microseconds | 198,706 | 77,972.5 |
| Outside-target named frontier CPU-seconds | 0.001232 | 0.001441 |
| Batch-step CPU-seconds | 0.052563 | 0.052460 |
| Signed root-minus-frontier CPU-seconds | 0.234512 | 0.106492 |
| Narrow read overhang, milliseconds | 1.380534 | 1.053092 |

Named outside-job cgroup CPU is small and is not elevated at the median in
the flagged group. Host process creation, context switching and native pressure
are higher. These associations are consistent with a workload-dependent
accounting or kernel-work contribution, but do not identify it. Root cgroups
contain direct tasks not assigned to the named frontier, and reads are
non-atomic. Do not subtract these differently bracketed quantities or treat
host process creations as authenticated native subprocess counts.

## OrthoFinder Counterexample

The OrthoFinder flag at index529 spans approximately528.992-529.988seconds.
It has only72 host process creations, 2,197 context switches, native average
19.7086cores, batch-step0.104649CPU-seconds, outside-frontier0.001458CPU-seconds,
root-minus-frontier0.217497CPU-seconds and1.828024ms narrow overhang.
Thus the high-creation pattern does not explain every retained flag. A general
claim that process creation causes all flags would exceed this evidence.

## Next Experiment

The successful earlier single-core controls did not characterize full-node
native execution with rapid process creation. A useful next prospective
experiment should compare bounded steady full-node work with bounded
process-creation work, alongside a known outside-target CPU positive control.
Record achieved native CPU/creation activity rather than assume matching
utilization. This can test whether native workload alone reproduces flags
and whether real outside work remains detectable. It cannot retroactively
identify the root tasks active in these runs.

Freeze conditions, repetitions, safety bounds and interpretation before
collecting those outcomes. Do not remove these flags, change thresholds after
seeing controls, or proceed directly to scientific scaling. The repeated
overhead budget and scientific inclusion policy remain separate requirements.

## Reproduction

```sh
python -m benchmark_tools.describe_dual_native_flags \
  --audit benchmark_tools/results/dual_native_audit_21912_20260919.json.gz \
  --repo . --output /new/path/dual_native_flags.json
```

The original archive paths referenced in the audit must be accessible.
Forty focused tests passed, covering conservative stage gaps, invalid durations,
proc-counter parsing, timestamp/flag mismatches, interval coverage and related
audit/reader behavior. Tests and this descriptive analysis do not admit
controlled timings or establish biological accuracy/publication readiness.
