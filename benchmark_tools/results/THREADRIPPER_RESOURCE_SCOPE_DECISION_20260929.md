# Timing Resource Scope Decision

Status: proposed, awaiting user response. Not a production amendment, execution
permit or waiver of any existing admission check. The 27 identities, scientific
settings, inputs, ordering and repeats remain unchanged.

## Decision Needed

The goal asks for matched hardware, actual CPU use and consistent peak-memory
accounting, with preparation, inference, conversion and scoring separated.
The implementation additionally pursued terminal whole-job accounting,
including controller/reporting teardown. These are different endpoints.

Recommendation: prospectively designate native-inference measurements as the
primary resource endpoints, with their exact instrumentation scopes disclosed.
Retain preparation/reporting measurements separately and final whole-job values
as unavailable. This could avoid requiring scheduler changes solely to obtain
teardown counters. It must not relabel partial job counters as complete job use.
An asynchronous question requests this scope decision; no response is assumed.

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

If this proposal is selected, write and freeze a prospective resource-endpoint
amendment before any production outcome is inspected. Bind it to the collector,
replayer, runtime and readiness review; test missing/decreasing counters,
scope mismatches, incomplete descendants and nonzero exits. Reuse validated
fixture evidence where its bindings remain unchanged. Do not manufacture a
positive readiness receipt or silently weaken existing guards.

## Independent Isolation Blocker

[Fresh read-only evidence](threadripper_resource_scope_load_20260929.json)
finds 104.586267 observed persistent foreign CPU-core equivalents across two
snapshots separated by a requested five seconds. The Slurm queue is empty.
Three sampling errors and eight unmatched identities are retained; this is
not a complete whole-host or whole-run measurement. The largest observed
processes are IQ-TREE instances. No command lines, environments or signals
were collected or sent. The sentinel exclusion scope was absent in both
snapshots, so only the observer PID was excluded from CPU totals.

Either scope decision still requires workload coordination, a reviewed service
policy, whole-run monitoring, full-scale observer validation and the tested
environmental handoff. No timing job or DGX operation was started.
