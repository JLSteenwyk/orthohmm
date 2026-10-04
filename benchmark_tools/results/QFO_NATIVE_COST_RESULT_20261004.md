# Retained Full Native QfO Cost Verified

## Actual Result

The [association](qfo_native_configuration_cost_20261004/association.json)
and [table](qfo_native_configuration_cost_20261004/association.md) link the
completed corrected-QfO high-sensitivity job **21707** to `p1_c0_r0`.
This is one retained full native inference observation, not a new run,
another repeat, or the cost of the original cached factorial execution.

The [checker](../link_qfo_native_cost.py) is committed/pushed at `8cdaff17`
before actual collection. Collection succeeds once. Rehash 133 distinct
direct records, including all 78 corrected FASTAs and 40 frozen source
records; verify the recorded command, CPU budget, native outcome and
native/replay checkpoint chain. Parse all input identifiers and both raw
partitions: 984,137 genes occur exactly once in 391,908 groups with exactly
equal whole-partition membership. The selected score's candidate/replay
metadata uses this same arm. This is stronger than relying on the old
partition-equality flag alone, but does not rescore predictions or prove
search completeness or a historical process trace.

| Measurement | Retained value | Scope |
| --- | ---: | --- |
| Native wall | 70,886.239038 s | Internal search through group materialization |
| Native user CPU | 2,133,115.343629 s | Recorded native process/descendant CPU |
| Native system CPU | 2,146.661854 s | Same native interval |
| Native peak RSS | 18,245,820,416 bytes | Sampled sum of process-tree RSS |
| GNU-time wall | 70,888 s | Native child launch-to-exit, rounded |
| GNU-time peak RSS | 11,885,516 KiB | Maximum process RSS, not tree sum |
| Scheduler elapsed | 19:41:39 | Historical job/wrapper interval |

The native metrics include initial search, graph inference, profile
expansion, refinement and group materialization. Preparation, parent
validation, conversion and scoring are excluded. Keep all distinct interval
and memory scopes; do not add stage times into a synthetic observed total.
Sampled summed RSS can double-count shared pages and miss unsampled peaks.
Neither historical memory measure is the later panel's lifetime cgroup peak.

## Verification

Thirty-three preparation cases pass on their first invocation. After actual
collection, add an independent raw-artifact readback and run the new tests
plus the existing native OrthoBench and QfO stage tests: **101 pass**, no
failures/skips, in 6.29 seconds. Independent readback rehashes every direct
record and recomputes whole-partition equality from raw named/unnamed group
files without calling the checker's collection or partition helper. It
also checks every native resource/stage field, the separate GNU companion,
and all sixteen original unavailable full-cost entries. The
[execution receipt](qfo_native_cost_execution_20261004.json) records exact
source, outputs and local JUnit identities.

Reproduction requires the retained local absolute-path inputs/runtime, not
only the reporting JSON. From the repository root, select an absent output:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  benchmark_tools/link_qfo_native_cost.py --root . \
  --output-directory /tmp/orthohmm-qfo-native-cost-readback
```

The output directory is refused if it exists. This command rechecks retained
evidence; it does not submit inference or modify historical results.

## Remaining Work

Together with the separate [native OrthoBench association](FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md),
these reports cover three prescribed configurations: two OrthoBench arms
with three repeats each, and this QfO arm with one historical observation.
Thirteen other configurations remain without linked full native costs in
this combined inventory. This is not a claim that no relevant measurement
exists elsewhere. Check exact bindings before commissioning new runs.
All sixteen full-cost fields for the original cached factorial executions
remain unavailable and unchanged; these separate native measurements do
not retroactively fill them.

The historical shared host has unknown, potentially method-dependent
contention. Its 32-CPU/192-GiB requested allocation and original deployment
are not the newer private-runtime 32-core/128-GiB panel. No background series,
matching-memory pass, isolated performance or causal speedup is imputed.
Do not pool these observations across deployments or infer tool-speed ranks.

No scoring, frozen method/default, completed scaling run, older table,
manuscript/PDF/archive snapshot or unrelated job/service is changed. Rc4
does not contain this later addendum. Broader uncertainty, independent
validation scope, prespecified strata, native restoration and distribution
requirements remain; the full goal stays active and readiness unproven.
