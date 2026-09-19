# Dual-Collector Overhead Submission

## Deployment

Launcher/batch/tests committed as `6e222a4`;39 focused tests and `bash -n`
passed. Exported452 regular files directly from committed Git objects. Every
remote recipe file's size and SHA-256 matches its source blob; no external
symlinks. Native workloads and frozen output/input order are unchanged.

- Recipe directory: `/home/jlsteenwyk/projects/orthohmm-publication/dual_overhead_recipe_v1`.
- Recipe manifest: `dgx_dual_overhead_recipe_20260919.json`.
- Recipe SHA-256: `ba18b0dca91698aa07cb451085c9340346828eaf86be90b7de53cfc2badceeb8`.
- Authorization: `dgx_dual_overhead_authorization_20260919.json`.
- Authorization SHA-256: `1179b6fc403b146d43aafb2363379c65c1a176afe442914046ec49d84ace4144`.
- Plan SHA-256: `1160a669a4033c66e7bdde5baddf429e83b9328d3821f99291fd29210cac8ab9`.
- Protocol SHA-256: `7df6f18476d20d9e96ab1a2f1fd69956e42347bf519ba93f35c0edeb1c4e761a`.
- Archive: `benchmarks/work/dual_overhead_recipe_6e222a4.tar`.
- Archive SHA-256: `fc438ab1d9c36aa52997e08dbde425e3deb0078c9b952a7bcad35e275432ca97`.
- Batch SHA-256: `2289d7050942c1f7ab81553c42fb6c83890349e9f646a646658e8bc1de7a8020`.

Both application and system runtime inventories verified on the DGX before
submission. Each task rechecks pinned inventories and recipe before/after
measurement. Queue was empty and both current one-second vmstat samples
reported100% idle, zero swap and steal. This is preflight evidence only.

## Live Panel

Submitted array21920, tasks0-17 with throttle1, exclusive20CPU/96GiB per
task, one-hour scheduler limit,900-second native timeout, no requeue.
Task0 (actual job21921) confirmed RUNNING on spark-7ff0 at00:00:13 with the
expected resources, working directory, throttle and zero restarts.

Submission explicitly clears inherited loader/Python overrides and supplies
remote TMPDIR:

```text
--export=ALL,LD_LIBRARY_PATH=,LD_PRELOAD=,LD_AUDIT=,PYTHONPATH=,TMPDIR=/tmp
```

The pinned batch specifies the remote working directory and interpreter,
disables bytecode writes, sets a fresh isolated cache prefix before Python
startup and restricts execution to the exact assigned array index. External
authorization covers all18 engineering tasks and no scientific execution.

Recorder21922 runs separately on bizon,1CPU/512MiB,20-hour scheduler limit,
71,400-second capture limit, using the same extracted committed recipe.
Evidence directory: `benchmarks/work/dual_overhead_scheduler_21920/`.
Log: `benchmarks/work/dual_overhead_capture_21922.log`.
No SSH inspection during this panel; controller polling only.

## Required Follow-Through

Wait for all assigned tasks and detailed terminal records, collect the complete
archive, then bind commands/runtime/resources and validate native outputs and
within-pair canonical equivalence. Retain all failures and numerical overhead
differences; no selective reruns. Summarize all original/narrow flags separately
from5% median/10% pair overhead budgets. Nothing in this submission establishes
environmental validity, scientific timing admission or publication readiness.
