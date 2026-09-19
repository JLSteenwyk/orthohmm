# Complete DGX Pressure Overhead Audit

## Outcome

Array `21889` is terminal: 16 completed tasks and two failed tasks (5 and
12). Recorder `21890` completed, preserving all 18 terminal scheduler
records across 2,489 polls, with zero observation errors or missing tasks.
All 18 recorded scheduler-file hashes were independently verified.
No DGX access occurred during the active panel. Collection began only
after local accounting confirmed all tasks terminal.

The existing `audit_frontier_overhead` auditor, using `--panel
pressure_21889`, validated all 16 successful tasks against the frozen
recipe, task and scheduler bindings, raw measurement replay, required
native-pressure observations and native-output fingerprints. Both failed
tasks remain failed. No temporal overlap or changed boot domain was found
among validated runs. Seven available paired comparisons had equivalent
canonical outputs and met the minimum 60-second duration requirement.

| Method | Available pairs | Median periodic/boundary change | Complete method numerical budget |
| --- | ---: | ---: | --- |
| OrthoHMM high-sensitivity | 3/3 | -0.710027% | Pass |
| OrthoHMM satellite-v2 phylogeny | 3/3 | +0.077954% | Pass |
| OrthoFinder full | 1/3 | Unavailable | Unavailable |

The one available OrthoFinder pair has change -0.892355%. Negative ratios
represent observed variability, not evidence of negative measurement cost.
There is no complete-panel numerical budget result. No timings were
corrected by subtracting overhead, and no failures were selectively rerun.

## Unresolved Validity

All validated whole-command CPU screens pass. However, all eight validated
periodic runs contain flagged intervals: tasks 1, 3, 6, 8, 10, 13, 15 and
17 have respectively 358, 527, 493, 177, 380, 175, 333 and 554 flags.
Boundary-only runs have no interval screening. Pressure observations are
diagnostics, not a replacement exclusion threshold or proof of contention.
These results do not establish environmental validity.

Tasks 5 (OrthoFinder periodic, pair 0) and 12 (OrthoFinder boundary, pair 2)
retain `verified_wrapper_failed` with error `Boot, target or frontier
identity changed`. Their native `done.json` records have exit code zero
and `timed_out: false`. Thus the scheduler failures must not be described
as demonstrated OrthoFinder inference failures. Determining the precise
identity change and the interpretation of interval residuals remains work
to do; the validation gates have not been relaxed.

**The 27 scientific scaling runs are not authorized or admitted by this
audit. Publication readiness remains unproven.**

## Evidence

- [Complete audit, gzip JSON](dgx_pressure_overhead_audit_21889_20260919.json.gz):
  gzip SHA-256 `d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006`;
  decompressed SHA-256 `d390edc3916e7ec5079e218f2343532328e96a0a2e4f1bf3ca9229c4ea0ff13b`.
- [Recorder capture](dgx_pressure_overhead_capture_21890_20260919.json):
  SHA-256 `fe9a0f74ae05b42b5788fc4cebeac949befccf98daec24a77ed215b09090e08a`.
- [Terminal accounting](dgx_pressure_overhead_accounting_21889_20260919.txt):
  SHA-256 `0ac50741994d5591e9d9c0f7753d6af021ef123cd541386eaf2c30c257babbea`.

Full local evidence resides at `benchmarks/work/pressure_overhead_archive_21889/`.
The main transfer copied 649,423 regular files (4,450,712,068 logical bytes)
from the output and recipe directories. Frozen four-proteome inputs,
scheduler records and all 18 batch logs were also collected. Large raw
outputs are not committed. Remote originals are unchanged.

Executed from the repository with the existing Python 3.10.13 environment:

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.audit_frontier_overhead \
  --panel pressure_21889 \
  --archive benchmarks/work/pressure_overhead_archive_21889 \
  --results benchmark_tools/results \
  --accounting benchmarks/work/pressure_overhead_archive_21889/accounting.txt \
  --output benchmarks/work/pressure_overhead_audit_21889_v1.json
```

The auditor exited zero. This is evidence validation, not a biological
accuracy assessment or an admission of controlled scientific runtimes.
