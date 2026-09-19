# Retained DGX Host-Bracket Replay

## Question and Method

The frontier collector preserves `hierarchy_host_after` before its frontier
and pressure reads, then replaces the final host snapshot with a later read.
The archived points therefore support a paired diagnostic: keep the exact
same native CPU counters, native timestamps, left host bracket and original
thresholds, but compare the two retained right host snapshots.

`diagnose_pressure_brackets.py` performs this comparison on all eight
validated periodic tasks from panel 21889. It first reproduces every outer
interval dictionary and flagged-index inventory exactly, then validates and
replays the earlier bracket. It checks native CPU differences and durations
are identical. Original files and measurement decisions are not modified.
Failed tasks 5 and 12 remain present without inferred interval results.

The source summary and eight reports are hash-checked before and after
reading. This is not a repeat of the full frontier/scheduler admission.

## Results

| Task | Method | Intervals | Original flags | Earlier-bracket flags | Earlier residual median (cores) |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | High sensitivity | 548 | 358 | 0 | 0.083773 |
| 3 | Satellite v2 | 815 | 527 | 14 | 0.083432 |
| 6 | Satellite v2 | 820 | 493 | 16 | 0.075573 |
| 8 | OrthoFinder full | 617 | 177 | 0 | 0.037028 |
| 10 | High sensitivity | 547 | 380 | 0 | 0.083952 |
| 13 | OrthoFinder full | 619 | 175 | 1 | 0.035474 |
| 15 | High sensitivity | 546 | 333 | 0 | 0.042587 |
| 17 | Satellite v2 | 790 | 554 | 15 | 0.084428 |

Across 5,302 intervals, the original screen flags 2,997 and the earlier
bracket flags 46. Exactly 2,951 change from flagged to unflagged; none
changes in the opposite direction. All 46 remaining flags are positive
`excess_unassigned_cpu` flags. Earlier-bracket median overhangs are
0.919-1.721 ms, compared with the wider collector's 9.467-18.884 ms.

The additional host CPU counted by the later endpoint has median
0.21-0.33 CPU-seconds per interval for OrthoHMM and 0.04 for OrthoFinder.
This establishes the arithmetic consequence of the bracket extension on
these archived counters. It does not establish who consumed that CPU,
causal interference, synchronized accounting, or a wall-time correction.
Overlapping host windows must not be summed as foreign CPU time.

## Consequences

The broad flags cannot be interpreted as direct evidence of competing CPU
use. However, this retrospective replay is not a replacement acceptance
policy: 46 residual flags remain, native/host accounting delays are still
possible, and two overhead tasks remain failed. No historical screen is
relabeled and no scientific scaling run is admitted.

A prospective collector can preserve the narrow hierarchy bracket for the
native CPU screen and keep frontier/pressure observations separately, each
with its own timestamps. Before timing admission, test that design with
known native-only and injected-load controls, including the remaining
burst behavior, and rerun the complete overhead panel under a frozen
policy. Do not selectively rerun only favorable or failed pairs.

## Validation and Reproduction

Eight new diagnostic tests plus the original frontier/interval tests pass:
36 tests total. Tests reject counter, chronology, boot, stored-result and
flag-inventory changes, ensure no input mutation, and isolate synthetic
added host ticks without changing native counters.

Retained result: `dgx_pressure_brackets_21889_20260919.json`, SHA-256
`53ad017d543d4b9bce2b610a6038fec14ce51b0c96adfd098cec9068c6e39ecd`.

```bash
python -m benchmark_tools.diagnose_pressure_brackets \
  --summary benchmark_tools/results/dgx_pressure_flags_21889_20260919_v2.json \
  --sha256 d560fa5d42458cd7914a7ce93dfc87f64ae2d19784a0f26d1c13a8a6ac255781 \
  --output /new/path/pressure_brackets.json
```

Requires the unchanged local raw-report archive. No DGX job was launched
and no unrelated process or service was stopped for this replay.
