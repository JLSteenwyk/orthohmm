# DGX Dual-Collector Overhead Panel

## Collection

All 18 tasks are terminal: 16 COMPLETED and two FAILED (`21920_3` and
`21920_17`, exit 1:0). Controller recorder `21922` completed exit 0:0.
The complete terminal gate was checked before any post-run SSH or native
output inspection. No jobs were retried or system services stopped.

The [controller replay](dual_overhead_scheduler_replay_21920.json)
checks all 2,505 immutable polls against the 18 first-terminal records;
there are zero observation errors and no missing tasks. Replay source
`replay_array_scheduler_capture.py` was committed/pushed at `1e1605a`.
All 31 collector/replay tests pass, including transient errors, expiry,
failed tasks, changed snapshots, missing polls and altered summaries.
Replay uses the collector's tested parser, not an independent parser.
It does not establish environmental or scientific timing validity.

```bash
python -B -m benchmark_tools.replay_array_scheduler_capture \
  --directory benchmarks/work/dual_overhead_scheduler_21920 \
  --array-id 21920 --tasks 18 \
  --output benchmarks/work/dual_overhead_scheduler_replay_21920.json
```

The original remote outputs remain intact. Archived native runs, recipe,
and four-proteome input are under
`benchmarks/work/dual_overhead_archive_21920/`, preserving paths relative
to `/home/jlsteenwyk/projects/orthohmm-publication`. Rsync transferred
648,063 regular files totaling 4,247,333,753 bytes, plus directory entries;
the 18 controller records and batch logs were collected separately.
Large raw outputs are not committed. SSH was performed only after the
quiet panel ended, over the dedicated Ethernet address `10.10.10.2`.

The original accounting download used `JobIDRaw`, following a preparation
note typo. The audit correctly rejected it before native metadata reads
because numeric allocation IDs do not identify array tasks. A separate
untouched [array accounting file](dual_overhead_array_accounting_21920.txt)
uses `JobID` and retains all tasks and steps. The preparation note is fixed;
no scheduler identity or validator rule was relaxed.

## Failed Measurement Evidence

Both native commands report exit 0 and `timed_out: false`. Both periodic
collectors preserved `failed_point.json` with status
`invalid_frontier_observation`, inner status `invalid_frontier_snapshot`,
and error `Cgroup frontier changed during sampling`. They produced no
complete `dual_bracket_report.json`.

| Task | Newly appearing cgroup between inventories | Native exit | Timing status |
| --- | --- | --- | --- |
| 3, satellite_v2 pair 0 periodic | `/system.slice/sysstat-collect.service` | 0 | Failed measurement |
| 17, satellite_v2 pair 2 periodic | `/system.slice/anacron.service` | 0 | Failed measurement |

There were no removed or identity-changed entries in either failed
snapshot, and ancestor direct-process counts were unchanged. This
localizes the collector failures; it does not quantify those services'
resource use or establish that they had no impact on native work.
Native exit success alone is not output or timing admission. Missing
interval observations cannot be reconstructed by deleting the invalid
point or substituting a later snapshot.

Both raw failed points are retained:
[task 3](dual_overhead_failed_point_21920_3.json),
[task 17](dual_overhead_failed_point_21920_17.json).
Their respective SHA-256 values are
`2e576016b7120e77b3f8f2dcbe975b2f62387e656437df83be4121526adbbb79`
and `d4234c72550405c210cce8e71acf44d41dc529bf09f7822d86016ffec576f57e`.
All partial points, native logs and completion records remain in the archive.

## Post-Run Audit

The existing frozen-plan/native-provenance, output-identity and measurement
replay audit is run against the complete archive, not only successful pairs:

```bash
python -B -m benchmark_tools.audit_frontier_overhead \
  --panel dual_21920 --archive benchmarks/work/dual_overhead_archive_21920 \
  --results benchmark_tools/results \
  --accounting benchmarks/work/dual_overhead_array_accounting_21920.txt \
  --output benchmarks/work/dual_overhead_audit_21920.json
```

This is a smallest-input collector-overhead panel, not the requested
matched-resource scaling comparison itself. Scientific timing admission
and publication readiness remain separate, unmet gates.

## Audited Results

The audit completed successfully. All 16 scheduler-successful tasks passed
native-output/provenance and measurement-replay checks. Both failed tasks
remain failed, not inferred successes. There were no temporal-order or
boot-identity issues. The full audit is retained as
[`dual_overhead_audit_21920.json.gz`](dual_overhead_audit_21920.json.gz)
(2,830,790 bytes, deterministic gzip). Compressed SHA-256:
`b8ecd2ffbd630a8912bbecdaa7c18822ad94d83112d4ae2f234f38dab55107ed`.
Decompressed JSON: 83,818,587 bytes, SHA-256
`aa259ceac0d4dbe09f408f5e2fab56678bd355c1593ab302168dfab77c9939bf`.
The source and helper hashes are embedded in that report.

Signed differences below are `100 * (periodic / boundary - 1)` using
native monotonic wall time. They compare collector modes within a method,
not accuracy or runtime superiority between methods.

| Method | Pair 0 (%) | Pair 1 (%) | Pair 2 (%) | Median (%) |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high-sensitivity | 1.721034 | 0.606718 | -0.816960 | 0.606718 |
| OrthoHMM satellite_v2 | Missing: task 3 failed | 1.313434 | Missing: task 17 failed | Not estimable for complete panel |
| OrthoFinder full | -0.523596 | 1.078953 | 1.041131 | 1.041131 |

All seven available pairs exceeded the prespecified 60-second minimum,
had equal canonical output identities, and met the 10% per-pair numerical
budget. The two fully observed methods met the 5% median budget. The
complete-panel budget is **null**, not passed: satellite has only one of
three pairs. Canonical output identity does not prove identical internal
work. Negative ratios reflect variability, not negative collection cost.

Whole-command screens passed for all 16 validated tasks, but the periodic
screens did not uniformly pass:

| Periodic task | Observation points | Original flagged intervals | Narrow-bracket flagged intervals |
| --- | ---: | ---: | ---: |
| 1 | 554 | 371 | 0 |
| 5 | 616 | 166 | 2 |
| 6 | 822 | 622 | 20 |
| 8 | 625 | 171 | 0 |
| 10 | 543 | 377 | 0 |
| 13 | 621 | 185 | 1 |
| 15 | 547 | 394 | 0 |

All 23 narrow-bracket flags have reason `excess_unassigned_cpu`. Both sets
of flags remain in the report; narrowing the observation bracket does not
erase the original screen. Boundary arms have only two observations and
cannot establish interval-level quietness. Neither screen causally assigns
the unassigned CPU to a service or measures its effect on inference.

The failed tasks' native completion records are also retained:
[task 3](dual_overhead_native_done_21920_3.json) and
[task 17](dual_overhead_native_done_21920_17.json).
No failed pair was selectively retried, no interval was dropped, and no
runtime was corrected by subtracting measured or estimated overhead.

The next timing step requires a prospective measurement/isolation protocol
that handles transient service cgroups without losing resource coverage,
followed by a newly frozen complete panel. The present archive cannot
recover the missing intervals and does not justify relaxing thresholds
after observing results. The 27-run scientific scaling comparison remains
unadmitted; this panel does not close the controlled-resource requirement.

## Compact Manuscript Summary

The [compact summary](dual_overhead_summary_21920_20260919.json) is a direct
projection of this already-audited report, not another admission or replay.
It retains failed task indices, incomplete-panel/null budgets, per-method
available-pair counts, numerical medians, all 23 narrow flags in validated
periodic tasks, and false scientific/environmental admission. The manuscript
and claim checklist now describe this completed result rather than the
historical running state.

After verifying the compressed audit SHA-256 against the value above,
reproduce the projection with jq using a fresh destination:

```sh
gzip -dc benchmark_tools/results/dual_overhead_audit_21920.json.gz |
  jq --arg sha256 b8ecd2ffbd630a8912bbecdaa7c18822ad94d83112d4ae2f234f38dab55107ed \
    -f benchmark_tools/summarize_dual_overhead.jq > /new/path/summary.json
```

The projection accepts the supplied digest as metadata; it does not itself
authenticate input bytes. Regeneration was byte-compared with the retained
summary after separately verifying the archive checksum. Original raw audits,
counter records and native outputs are unchanged.
