# Lineage Overhead Launcher

`run_lineage_overhead_panel.py` binds the frozen plan
`90bcdb4b12655270dad2d69a3806174a4c4a63efb6e400b530e671b99938b1ed`
and protocol
`6107283f87892c801408a1717eaba42d3dc104650b4603875281d3d08ac35893`.
An external hash-pinned authorization must cover all 18 engineering tasks
and explicitly exclude scientific execution. The pinned deployment inventory
must contain every local benchmark Python module, the batch script, plan and
protocol. Existing completed scientific executors are unchanged.

The launcher follows the existing overhead workflow: exact DGX interpreter,
host, assigned array index and 20-CPU/96-GiB resource guards; isolated fresh
cache and disabled bytecode; no loader/Python overrides or system preload;
runtime checks before/after native execution; exact native input ordering;
and a task receipt retaining the plan, recipe and authorization hashes.
Periodic and boundary tasks select their corresponding lineage collectors
without the older frontier collector or optional pressure configuration.

`run_dgx_lineage_overhead.sh` defines an exclusive sequential 0-17 array,
no requeue, one-hour scheduler limit and the pinned remote working directory.
Submission must clear inherited loader overrides and use `TMPDIR=/tmp`, as
documented by earlier preserved launch failures. Start a controller recorder
before releasing held tasks, and do not inspect native outputs until every
assigned task is terminal. No submission has occurred at this milestone.

All 99 targeted tests pass:

```sh
python -m pytest -q tests/unit/test_run_lineage_overhead_panel.py \
  tests/unit/test_prepare_lineage_overhead_panel.py \
  tests/unit/test_run_dual_overhead_panel.py \
  tests/unit/test_measure_lineage_boundary_step.py \
  tests/unit/test_replay_lineage_boundary_measurement.py
bash -n benchmark_tools/run_dgx_lineage_overhead.sh
```

Tests cover all 18 selections, real checksum rejection, exact authorization
scope, changed source/protocol/inventory, both collectors and receipts, and
rejection before measurement of incorrect host/interpreter/resources/index,
existing cache, enabled bytecode, hash seed or loader override. Worker launch
tests are mocked, not live scheduler or kernel validation. Remote deployment
identity, runtime verification, full-panel execution and post-run audit remain
required. Overhead-budget passage would not establish clean scientific timing.
