# Native Root-Context Diagnostic 22021

## Outcome

All three frozen four-proteome native commands completed and passed raw
lineage/root-context replay, provenance and native-output validation. Their
canonical outputs equal the prior validated lineage diagnostics. No panel
order, boot or named-scope identity issue was detected. This establishes
native integration of the supplementary context collector, **not scientific
timing admission, collector overhead or a method-speed ranking**.

| Native method | Native wall seconds | Intervals | Original flags | Narrow flags | Output equivalent |
| --- | ---: | ---: | ---: | ---: | --- |
| OrthoHMM high sensitivity | 551.830510409 | 552 | 0 | 0 | Yes |
| OrthoHMM satellite_v2 | 809.677239777 | 810 | 51 | 10 | Yes |
| OrthoFinder 3.1.5 full | 613.159839346 | 614 | 1 | 1 | Yes |

All commands used the same 73,266 proteins and frozen method settings.
Native output counts were 35,314 orthogroups for high sensitivity;
35,627 orthogroups/root HOGs and 51,256 native pair rows for satellite_v2;
24,052 checkpoint groups and 90,490 native pair rows for full OrthoFinder.
These are different native output types, not comparable accuracy endpoints.

GNU-time user/system CPU seconds were 9452.29/5.80, 12329.42/758.43 and
4346.11/148.73, respectively. Native-step cgroup `memory.peak` values were
2,267,328,512, 2,462,818,304 and 6,936,137,728 bytes. All retained native-step
memory-event counters were zero. GNU-time maximum-process RSS is retained
separately and must not be conflated with simultaneous cgroup memory.

## Session And Environment

Slurm job 22021 completed `0:0` in 00:33:34, from
2026-09-19 19:53:58 to 20:27:32 America/New_York, in one exclusive
20CPU/96GiB allocation with a one-hour limit and no requeue/restart. The
bounded SSH `sbatch --wait` receipt passed independent comparison to the
terminal scheduler record. Remote/local bounds remained 3720/3750 seconds.
No extra DGX SSH inspections or transfers occurred during the panel.

The bounded user journal covers 23:53:50 UTC through 00:27:40 UTC. Its
printed timestamps use UTC-07:00, while the scheduler record uses UTC-04:00.
The manager startup is recorded at the job start, with no manager shutdown
reported in that window. The unrelated `samwise-daemon-samwise.service`
scheduled 383 restarts during the printed job window after failing because
its working directory was missing. Two additional restart entries occur
after the printed job end within the retained journal margin. No unrelated
service or persistent login setting was changed. This is not a background-free
environment; journal co-occurrence does not establish the cause of CPU flags.

## Descriptive Context

The report retains all 1,976 intervals, their 52 original and 11 narrow
flags, and original/narrow flagged/unflagged subsets. Empty subsets have
missing distributions, not zero-valued estimates. Root-membership changes
were observed in 15, 22 and 17 intervals, respectively.

For the root-minus-three-named-children CPU residual, all-interval medians
were 1,336.5, 1,529 and 303 microseconds. Narrow-flagged medians were missing,
228,474.5 and 178,814 microseconds. Signed minima remain negative and are
not clipped. These are descriptive counter differences over non-atomic,
distinct read windows, not root-task CPU attribution. Every host category
and scope window is reported separately; overlapping guest/user categories
are not summed, and dependent intervals are not treated as replicates.

## Provenance And Retention

Execution source: `9821c81`; recipe SHA-256:
`2190269735646d5996832f841c4439d38e0b7f05fc5b3402af000f8b2bfe73f4`.
The post-run reporter was added locally in `abdd117`, without changing the
frozen DGX execution tree. Prior output comparison uses the pinned lineage
audit SHA `0e3a77f785b645f44e0ca5cb8bbaeffbdb7106e2104350045e299c658de9050a`.

| Retained artifact | SHA-256 |
| --- | --- |
| `root_context_native_audit_22021_20260919.json.gz` | `9b2a243a13cd83ab431025a7f24d5b166672d50cbaa581697a2eba2557e54b23` |
| `root_context_native_description_22021_20260919.json.gz` | `7ddf1fda24df7077c50ff19cd69a2431f6bdedcb60df33a23328473cd19bd7b7` |
| `root_context_native_receipts_22021.tar.gz` | `838560c3e43dccc5dcba4782a521df1a43b0ccd7e9074b04876b263ecfee5805` |
| `root_context_native_session_audit_22021_20260919.json` | `2f3c0c7d454a90680391b98521439b8b6091838c761186cbbd52888e8cf380b0` |

Large raw archives remain outside Git in `benchmarks/work/` and on the DGX:
`root_context_native_archive_22021.tar.gz` (SHA
`f1c0094d76951e2e67c8df29a34d8cd4ac50f778e8177dbf5b03d147ac5eee8f`)
and `root_context_native_inputs_22021.tar.gz` (SHA
`91830e6740313a74214a447db0aff379fe3913fd5a3e025af9e6c200d786dc9f`).
The receipt tar includes the bounded manager journal (SHA
`bf7393dda20769816de7307d7c22f9e2729378a27a477431c312f2cb32e1c63c`).
No external archival deposit or input redistribution authorization is implied.

An initial local audit invocation preceded completion of archive extraction
and rejected the incomplete recipe inventory. After extraction finished,
the audit passed without changing native evidence or rerunning any method.

A fresh relocated archive also validates all three tasks without panel issues.
Its result projection (native timing/accounting, screening, root context,
output identity/counts and observation bounds) exactly equals the first audit;
only location-dependent evidence records are omitted from that comparison.
`root_context_native_relocation_22021_20260919.json` records both uncompressed
audit hashes. The relocated raw audit remains in `benchmarks/work/`.

## Reproduction

Extract both raw archives into a fresh archive directory. Retain the prior
lineage audit's referenced evidence, as it is independently rechecked.
From the repository root, run (with fresh output filenames):

```bash
python -B -m benchmark_tools.audit_root_context_native \
  --archive benchmarks/work/root_context_native_archive_22021 \
  --results benchmark_tools/results \
  --recipe benchmark_tools/results/dgx_root_context_native_recipe_20260919.json \
  --recipe-sha 2190269735646d5996832f841c4439d38e0b7f05fc5b3402af000f8b2bfe73f4 \
  --scheduler benchmarks/work/root_context_native_scheduler_22021/scheduler_22021.txt \
  --job 22021 \
  --prior-audit benchmark_tools/results/dgx_lineage_native_audit_21995_20260919.json.gz \
  --output benchmarks/work/root_context_native_reaudit_22021.json
python -B -m benchmark_tools.describe_root_context_native \
  --audit benchmark_tools/results/root_context_native_audit_22021_20260919.json.gz \
  --audit-sha 9b2a243a13cd83ab431025a7f24d5b166672d50cbaa581697a2eba2557e54b23 \
  --output benchmarks/work/root_context_native_redescription_22021.json
```

The original supplementary collector overhead remains unmeasured in a matched
native comparison. Non-CPU isolation and eligibility for prospective scaling
measurements also remain unresolved. No historical scientific score, method
default or CPU-screen threshold changes follow from this diagnostic.
