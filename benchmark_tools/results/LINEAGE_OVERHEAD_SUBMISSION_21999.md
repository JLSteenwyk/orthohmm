# Lineage Overhead Panel Submission

## Deployment

Launcher and tests were committed and pushed as `9632cc7`; all 99 focused
tests and batch syntax checks passed. Exported 496 regular files directly
from that commit, including benchmark Python modules, the batch script and
frozen plan/protocol. Remote archive SHA-256 matched locally; all extracted
file sizes and hashes matched the committed archive, with no symlinks.

- Archive: `benchmarks/work/lineage_overhead_recipe_9632cc7.tar`.
- Archive SHA-256: `071a0e2263380a0149c246adda2b735f8cb2582eb3de5bcbd4e1e99de7d67906`.
- Remote directory: `/home/jlsteenwyk/projects/orthohmm-publication/lineage_overhead_recipe_v1`.
- Recipe: `dgx_lineage_overhead_recipe_20260919.json`.
- Recipe SHA-256: `73ee9571228868989dd5e687f7c30af518961885a4fe4dfdc919a7bf1820fe75`.
- Authorization: `dgx_lineage_overhead_authorization_20260919.json`.
- Authorization SHA-256: `0847cb9ba6d09349f823b795e77ff01f3500a9d31f32b457a10180c82ddcf02f`.
- Batch SHA-256: `c5e52b94cba1846b3aa8ef08a19f33fd9eb21313abb304a56a192ca771bf52ab`.

Plan/protocol hashes remain those in `LINEAGE_OVERHEAD_LAUNCHER_20260919.md`.
Remote selection accepted all 18 tasks and confirmed the output root absent.
Both pinned runtime inventories passed verification: 26,673 application-tree
records and 10,066 system-tree records. The DGX queue was empty, and two
current one-second vmstat samples showed 100% idle, zero swap and steal.
These are preflight observations, not evidence of interval-level quietness.

## Submission And Recording

Submitted array 21999 on hold, all tasks 0-17, throttle 1, exclusive
20 CPUs/96 GiB per task, no requeue and one-hour scheduler limit. Native
timeout remains 900 seconds. Submission removed inherited LD_PRELOAD,
LD_LIBRARY_PATH, LD_AUDIT and PYTHONPATH and exported `TMPDIR=/tmp`.

Controller recorder 22000 runs on bizon, one CPU/512 MiB, 20-hour scheduler
limit and 71,400-second capture deadline. It uses the same committed source
archive extracted locally under
`benchmarks/work/lineage_overhead_controller_9632cc7/lineage_overhead_recipe_v1`.
It was confirmed RUNNING with successful polls while array 21999 was held.
Only then was the full array released. Latest check: task 0 RUNNING at 29
seconds on spark-7ff0; tasks 1-17 pending at the array throttle; recorder
RUNNING at 1:10. No native output has been inspected after release.

Capture directory: `benchmarks/work/lineage_overhead_scheduler_21999/`.
Recorder log: `benchmarks/work/lineage_overhead_capture_22000.log`.
Remote outputs: `lineage_collector_overhead_v1/run_00` through `run_17`.
Remote logs: `lineage_overhead_21999_INDEX.log` in the DGX project root.

## Follow-Through

No SSH or native-output inspection until all 18 tasks are terminal. Preserve
complete scheduler records, failed tasks and missing evidence. Then collect
the complete output/deployment archive and audit exact commands, runtime,
input order, resource identity, native outputs, raw collector replay and
within-pair canonical equivalence. Report all nine signed overhead pairs,
per-method medians and every original/narrow flag; no selected replacements.
The old frontier/dual replay is not a substitute for the lineage schemas.
Boundary interval coverage remains unavailable, not clean by default.

This engineering panel does not authorize scientific scaling or establish
environmental validity. The remaining CPU discrepancies, during-read
lifecycle behavior and prospective scientific inclusion policy remain open.
Publication readiness is not established.
