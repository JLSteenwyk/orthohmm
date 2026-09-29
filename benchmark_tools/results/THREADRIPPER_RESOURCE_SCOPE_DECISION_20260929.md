# Timing Resource Scope Decision

Status: endpoint choice adopted prospectively under the goal's authorization
for autonomous engineering and phase-separated resource measurement. This
supersedes the earlier optional question; no user reply or approval is inferred.
The [machine-readable amendment](threadripper_resource_endpoints_20260929.json)
freezes the primary scopes and exact derivation sources before any production
measurement directory exists. It grants no execution permit or readiness pass.
The 27 identities, scientific settings, inputs, ordering and repeats remain unchanged.

## Adopted Endpoints

The goal asks for matched hardware, actual CPU use and consistent peak-memory
accounting, with preparation, inference, conversion and scoring separated.
The implementation additionally pursued terminal whole-job accounting,
including controller/reporting teardown. These are different endpoints.

Prospectively designate native-inference measurements as the
primary resource endpoints, with their exact instrumentation scopes disclosed.
Retain preparation/reporting measurements separately and final whole-job values
as unavailable. This could avoid requiring scheduler changes solely to obtain
teardown counters. It must not relabel partial job counters as complete job use.
The earlier optional scope question is resolved by this engineering decision.
Terminal whole-job teardown counters are supplementary unavailable fields;
their absence alone does not prevent measuring the selected phase endpoints.

| Quantity | Existing evidence | Required qualification |
| --- | --- | --- |
| Native wall time | Command start/end monotonic timestamps | Includes the entire native command, not installation or scoring |
| CPU seconds | Cgroup `cpu.stat` difference bracketing the command | Includes wrapper bracket work; preserve task-subtree identity and validate descendant accounting |
| Memory bytes | Kernel `memory.peak` for the native step, read after command completion while its anchor remains alive | Step-lifetime peak including launcher, not process RSS, command-only allocation or final whole-job peak |
| Preparation/reporting job peaks | Separate cumulative job-cgroup observations | End at their read boundaries; overlapping peaks cannot be added or subtracted |
| Terminal whole-job use | Not established | Remains missing; blank/zero Slurm fields are not measured zero |

The CPU task subtree and memory step are not identical scopes. A production
protocol must explicitly validate that the native command and its workers are
accounted for and disclose wrapper/step overhead. Boundary-only anchor checks
do not prove continuous containment. The three retained fixtures demonstrate
small-workload operation, not full-scale overhead or production eligibility.
See [CPU scope evidence](THREADRIPPER_NATIVE_CPU_SCOPES_20260929.md) and
[native fixture integration](THREADRIPPER_PRIVATE_COLLECTOR_PANEL_20260928.md).

The amendment is now frozen. `benchmark_tools.derive_threadripper_resources`
checks its externally supplied digest and source identities, replays all raw
measurement evidence and derives the three endpoints. Its output preserves
CPU/memory scope details and failed native outcomes, and cannot grant timing
admission. [Actual retained-fixture replay](threadripper_resource_endpoint_fixture_replay_20260929.json)
uses the three existing measurements; no fixture or inference was rerun.
Bind this amendment into the final collector/runtime/readiness review before
production. Existing execution guards remain intact. Continuous containment,
full-scale overhead and the environmental handoff still require validation.

```sh
python -B -m benchmark_tools.derive_threadripper_resources \
  --directory /local/measurement --job JOB_ID \
  --protocol benchmark_tools/results/threadripper_resource_endpoints_20260929.json \
  --protocol-sha256 dba9602c365184f7e802425b9116b8021175b225d1f2266c44077da924bcfcdc \
  --output /fresh/scoped-resources.json
```

## Independent Isolation Blocker

[Fresh read-only evidence](threadripper_resource_scope_load_20260929.json)
finds 104.586267 observed persistent foreign CPU-core equivalents across two
snapshots separated by a requested five seconds. The Slurm queue is empty.
Three sampling errors and eight unmatched identities are retained; this is
not a complete whole-host or whole-run measurement. The largest observed
processes are IQ-TREE instances. No command lines, environments or signals
were collected or sent. The sentinel exclusion scope was absent in both
snapshots, so only the observer PID was excluded from CPU totals.

Production still requires workload coordination, a reviewed service
policy, whole-run monitoring, full-scale observer validation and the tested
environmental handoff. No timing job or DGX operation was started.
