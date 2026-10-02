# Private Controller Startup

The preceding [runtime refresh](THREADRIPPER_RUNTIME_REFRESH_20261002.md)
validated native interpreter imports, not the controller's imports of the
current executor/environment helpers. Execute one bounded controller probe
using its unchanged binary and the existing lookup inspector's PROBE.
No source change, inference, allocation or host-contention sampling occurs.

## Actual Evidence

The pinned controller imports the current executor, environment worker,
lifecycle manager and native-overhead auditor. Child-only environment settings
match the launcher's bytecode/user-site/hash-seed settings; loader overrides and
PYTHONHOME are removed from that child, not the shared environment. This is a
standalone startup, not execution of the Slurm submission script.

| Check | Actual result |
| --- | --- |
| Child exit / stderr | 0 / empty |
| Imported file-backed modules | 410 |
| Project modules from the current benchmark_tools directory | 81 |
| Unique observed imported/mapped/executable files | 424 |
| File-backed mappings | 62 |
| Missing, changed or ambiguous pinned-file identities | 0 / 0 / 0 |
| Pure explicit-deployment path/digest/plan resolution | Passed |
| Isolated bytecode cache created | No |

The child reads only its own /proc/self/maps for the lookup observation, not
other processes. No measurement/observer function or executor execution path
is called. Pure deployment resolution is not selection of an authorized run.
The former full 57,958-record tree check remains its dated result and is not
repeated or claimed freshly executed here.

Separate stdlib readback parses captured process stdout, reconstructs the
unique file set and project origins, rehashes all 424 observed files, and
requires a single matching expected size/hash from the two pinned manifests.
This is independent count/hash/multimap checking, not independent loader
semantics, all-file-I/O tracing or continuous enforcement.

The [compact receipt](threadripper_controller_startup_20261002.json) pins nine
evidence files: 4,096 bytes, SHA256
`b5f3dcfbed3baa94bcce18f4b0b06a5246e214849cbffe4f01e4ede74d4f556f`.
The local complete observation is 126,628 bytes, SHA256
`efa22c377b47c74216d2c797a8c79436af5dd15d307967f00689c258732d1074`.
The first tool-wrapper input had a JavaScript syntax error before any nested
tool executed; no earlier startup attempt or filesystem change resulted.

## Remaining Work

This closes the missing observed controller-startup evidence, not all workload
branches, complete controller file I/O, the real executor/environment handoff,
causal observer slowdown, current readiness or a quiet window. No tests are
rerun because no implementation changes. Frozen scientific settings, prior
diagnostics, private runtimes and native lookup receipts stay unchanged.
No timing job, DGX access, service action, unrelated-job interruption, archive
rebuild or manuscript rerender occurs. All owned handles are terminal.
The original publication goal remains active and incomplete.
