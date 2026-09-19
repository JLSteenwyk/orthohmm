# Prospective Native Lineage Diagnostics

Run exactly three sequential native diagnostics on exclusive spark-7ff0
allocations, each 20 CPUs/96 GiB, one-hour scheduler limit, no requeue and
900-second native timeout. Use the existing frozen four-proteome input
(73,266 proteins), commands and runtime manifests from pressure-panel v2
tasks 1, 3 and 8, respectively OrthoHMM high sensitivity, OrthoHMM satellite_v2
and OrthoFinder full. Relocate only per-run output/cache paths and diagnostic
task metadata. Do not change method settings, inputs or inference code.

Use the new `measure_native_lineage_step` at one-second cadence, with its
separate lineage schema, both unchanged CPU screens, native pressure, native
timer and native-step cgroup memory accounting. Preserve each failed point,
command failure and timeout. These are live integration diagnostics, not
overhead measurements or admitted scientific scaling runs.

Freeze the protocol and derived plan before execution, deploy a new source
directory and hash its complete contents. Verify the pinned runtime and inputs
before and after each run through the existing `measure_run` path. Submit
all three jobs before releasing the first, using afterany dependencies to
retain all outcomes and enforce order even on failure. No selective retries.
Use fresh absent output roots and cache prefixes. Do not modify unrelated
jobs, system services or previous archives.

Inspect only scheduler state from the controller while any diagnostic is
active. No SSH/native-output inspection during the three-run panel. After
all jobs are terminal, retrieve complete raw measurements and native outputs;
replay the new schema and verify source/input/runtime, scheduler identity and
native output semantics before interpreting results. Preserve all original
and narrow flags. A native exit of zero is not sufficient evidence.

Report all three execution outcomes, native wall/CPU/memory/pressure,
measurement replay outcomes, output validation and unresolved environmental
uncertainty. Do not compare these one-off durations as method superiority,
subtract estimated collector overhead, or promote diagnostic success to
timing admission. During-read service churn, a complete fresh overhead panel,
non-CPU isolation and a prospective scientific inclusion policy remain open.
