# Native CPU Evidence Is Present

The [resource summary](threadripper_fixture_resource_scopes_20260929.json)
extends the retained private-runtime fixture memory summary with CPU counters
from the actual native-command read brackets. No job was rerun. These small
fixtures ran on a shared host and are not controlled efficiency measurements.

| Diagnostic job | Native path | Bracket CPU seconds | Command wall seconds |
| --- | --- | ---: | ---: |
| 22373 | OrthoHMM high sensitivity | 39.725260 | 12.807817 |
| 22374 | OrthoHMM phylogeny | 44.858304 | 16.324126 |
| 22375 | Full OrthoFinder | 22.036450 | 13.708840 |

CPU and wall intervals differ deliberately: CPU counter reads enclose the
command, so the difference includes wrapper work within that bracket. CPU
time can exceed wall time with parallel work. Do not interpret this table as
a speed ranking, compute a corrected runtime, or relabel the counters as
terminal whole-job CPU. The missing scheduler fields are still unavailable,
not measured zero. Native and job memory scopes remain unchanged.

The [kernel interface documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html)
defines `cpu.stat` usage/user/system counters in microseconds for the cgroup
and its descendants. `memory.peak` measures the cgroup/descendant peak since
creation or reset, not process RSS. These interface semantics do not by
themselves prove containment of a particular workload or absence of overhead.
The retained completion check establishes anchor-only boundary observations,
not continuous containment. No final-job accounting gate is waived here.

## Validation And Replay

The existing summary helper now requires exactly one retained `done.json`,
valid command/read ordering, unchanged host and cgroup identities, complete
nondecreasing usage/user/system counters, and agreement with the native
completion subtree. Evidence files are rehashed before and after derivation.
No CPU component is reconstructed by adding rounded user/system values.

Eighteen focused tests pass. All three actual retained fixture replays pass,
and separate direct parsing of the raw CPU strings reproduces all nine
reported counter differences. The first replay rejected a path-format
mismatch in the new helper; the correction maps absolute cgroup filesystem
paths to membership paths explicitly. No failure or raw record was discarded.

```bash
python -m pytest -q tests/unit/test_summarize_fixture_memory_scopes.py
python -B -m benchmark_tools.summarize_fixture_memory_scopes --repo . \
  --output /tmp/threadripper-resource-scopes-new.json
```

Use a fresh output path. The original memory-only receipt remains historical.
The prospective 27-run panel still needs a quiet window, validated complete
accounting and environmental handoff, source freeze and overhead validation.
No scheduler or service change was made.
