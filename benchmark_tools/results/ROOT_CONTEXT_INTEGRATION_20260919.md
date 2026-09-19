# Root Context Collector Integration

The existing lineage measurement loop now accepts an optional point reader.
Its default resolves to the original reader; worker launch, native timer,
completion/release handling, cadence, memory collection and original/narrow
screen calculations are otherwise unchanged. Frozen DGX recipes were not
modified or redeployed.

`measure_native_root_context.py` supplies the new reader: it completes each
ordinary lineage point, then collects a separately timed root context.
Validation binds boot identity and clock tick rate and requires supplementary
reads to follow the native observation and precede the next one. Context
windows are not substituted for the original native windows. Partial root
read failures are retained in `failed_root_context.json`.

The ordinary lineage report remains present, and a separate
`root_context_report.json` records context comparisons with a byte/hash
binding to that report. The binding uses a relative filename so relocated
archives can replay. `replay_native_root_context.py` first invokes the full
lineage replay, then independently validates and reconstructs the context
comparisons. Original flags, negative complements and unknown coverage are
not changed by the added observations.

All 133 focused tests pass across the new wrapper/replay, existing lineage
measurement/replay, boundary replay, root probe and complete-panel auditor.
Checks include unchanged original screens, timer/resource evidence, explicit
default-reader equivalence, malformed/missing context, identity/read-window
conflicts, partial-failure retention, modified reports and archive relocation.
The measurement integration uses simulated process fixtures; these tests do
not prove DGX execution, live overhead, workload validity or specificity.

Next implement the fixed 12-trial coordinator and independent workload audit
from `ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md`, then freeze/export
committed sources before deployment. No native control has run with this
supplementary collector, and no scientific timing is admitted. Added observer
work may affect native execution even though the timer definition is unchanged.
