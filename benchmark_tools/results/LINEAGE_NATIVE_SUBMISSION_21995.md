# Native Lineage Diagnostic Submission

## Frozen Inputs And Deployment

Source/plan commit `d1462ce` derives the same three native method commands
and four-proteome inputs from pressure-panel v2 tasks 1, 3 and 8, with new
output/cache paths under `lineage_native_v1`. The collector is the separate
native lineage implementation introduced at `6166b83`. No inference setting,
dataset, timeout, CPU threshold or runtime manifest changed. All 74 focused
launcher/collector/replay tests and batch shell syntax checks passed.

Bindings:

- Protocol `LINEAGE_NATIVE_PROTOCOL_20260919.md`, SHA-256
  `7a20fe9e759e10c8e42c269d66775c28b1d3da763e8b73bdac9c18b4ee21d00b`.
- Plan `dgx_lineage_native_plan_20260919.json`, SHA-256
  `714ef04458904e1c01526d171d0b2fde658ac135bd416af49ab107dc5ed14bfe`.
- Recipe `dgx_lineage_native_recipe_20260919.json`, SHA-256
  `d05c1028f486d5e0a9788912fc8d3785a7a2ca504e9caafbe5e8094a7a7dc24b`.
- Source archive `benchmarks/work/lineage_native_recipe_v1.tar`, SHA-256
  `43a609173e7dfc44eef9289d7d2450c67e90abd7e9956a634fb4c6d592cf3129`.

The transferred archive hash matched remotely. All 491 deployed file hashes
were compared with the corresponding frozen Git blobs and matched. Remote
selection preflight accepted all three exact tasks and confirmed the output
root absent. Runtime/input checks still execute before and after each native
command; deployment checks are not evidence that those later checks pass.

## Preserved Launch Failures

Initial tasks 21990, 21991 and 21992 all failed at the loader-environment guard,
before native launch. The submission shell exported `LD_LIBRARY_PATH` with
CUDA library directories; no `/etc/ld.so.preload` was present on the DGX.
Each task's log identifies the same guard failure. Logs are retained as
`lineage_native_21990.log`, `lineage_native_21991.log` and
`lineage_native_21992.log`. The output root was confirmed still absent.

Recorder 21993 also failed because invoking `capture_job_scheduler.py` by
file path did not place the package on Python's import path. Its log is
retained. Corrected module invocation in recorder 21994 captured all three
terminal records after completion: one poll, zero observation errors, no
missing jobs. This is terminal capture, not continuous observation of the
failed launch panel. Its raw records are in `lineage_native_failed_scheduler_21990/`.

The complete three-task panel was resubmitted with the same source/plan/recipe
and still-unused outputs. No guard was relaxed. Each submission explicitly
removed loader and Python overrides from the exported environment:

```sh
env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT -u PYTHONPATH \
  sbatch [hold/dependency options] benchmark_tools/run_dgx_lineage_native.sh \
  d05c1028f486d5e0a9788912fc8d3785a7a2ca504e9caafbe5e8094a7a7dc24b INDEX
```

The launch correction was made before any native outcome existed. All methods
were retained; no settings or response thresholds were tuned from results.

## Active Panel

| Job | Method | Dependency |
| --- | --- | --- |
| 21995 | OrthoHMM high sensitivity | Initially held, then released |
| 21996 | OrthoHMM satellite_v2 | afterany:21995 |
| 21997 | OrthoFinder full | afterany:21996 |

All request exclusive spark-7ff0, 20 CPUs/96 GiB, no requeue, one-hour
scheduler limit and 900-second native timeout. The queue was empty before
the initial submission. All three corrected jobs were submitted before
release. Controller recorder 21998 runs on bizon with one CPU/128 MiB,
four-hour capture deadline and 4:10 scheduler limit, using:

```sh
/usr/bin/python3 -B -m benchmark_tools.capture_job_scheduler \
  --jobs 21995 21996 21997 --output benchmarks/work/lineage_native_scheduler_21995
```

Recorder 21998 was confirmed RUNNING and had successful polls before job
21995 was released. The first post-release check showed 21995 RUNNING and
21996/21997 dependency-pending. No SSH or native-output inspection is permitted
until all three become terminal. Do not infer native success from RUNNING.

After terminal capture, collect complete native outputs and raw measurement
archives, validate provenance and output semantics, and replay the new schema.
All flags, failures and missing evidence must remain visible. These runs do
not measure incremental collector overhead and cannot admit the scientific
scaling panel. The full publication goal remains incomplete.
