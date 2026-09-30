# Calibration Reporting Failure

Job 22378 is terminal FAILED (1:0); its native step is COMPLETED (0:0).
All 32 witness processes finished, but reporting raised
`Missing or irregular observation interval`. No complete lineage/resource
report was written. The deferred independent audit in job 22379 therefore
also failed, reporting the missing lineage report. Neither attempt is admitted.

The [raw-evidence diagnosis](publication_terminal_readback_20260930.json)
binds all 32 raw points, process stream, native receipt, command, complete
witnesses and original failure log. All interval midpoints except the final
pair are within the unchanged 0.5-1.5-second range. Points 30 to 31 have a
3.6452974795-second gap. A synchronous host scan occupies that gap:

| Host scan | Duration | Process rows | Retained errors |
| --- | ---: | ---: | ---: |
| Initial | 3.324138 seconds | 1,988 | 2 |
| Periodic | 3.507657 seconds | 2,026 | 38 |
| Final | 3.247997 seconds | 1,993 | 3 |

The periodic scan starts at monotonic 368351.437196425 and ends at
368354.94485372. Point 30 starts 29.986209111 seconds after native launch;
point 31 starts 33.632393083 seconds after launch. Native execution lasts
31.41512004 seconds. Thus the delayed gap overlaps native execution; it
cannot be discarded as purely post-command reporting or waived by increasing
the cadence limit. Process-race errors also remain, not quiet-host evidence.

Before editing the local collector, all 823 source pins from the original
deferred protocol were checked and copied exactly to the failed run's local
`frozen_sources/` directory with an index. These are retained source bodies,
not a rerun or an executable portable study archive. Original receipts and
raw samples remain unmodified. Historical validators require their recorded
source revision, not the changed current collector.

## Prospective Change

The local collector now keeps periodic host scans on one owned observation
thread. Initial inventory still precedes native release; the thread starts
before `go.json`, and any release-guard freshness check includes that startup.
All periodic scans are serialized. It joins the thread before making a final
post-command host sample and summarizing, with a 15-second bounded join.
Unexpected thread failure or an unfinished join raises, not a quiet result.
Missed host periods are not followed by a burst of catch-up scans.

Only the local collector and its focused tests change; historical DGX helpers,
scientific inference settings, resource scopes and cadence thresholds remain
unchanged. The serialized report schema remains v5 because its fields and
interpretation are unchanged; source pins distinguish the scheduling revision.
173 focused tests pass across observer/collector/replay, calibration/resource,
executor/lifecycle, scientific admission and deferred-audit components.

A separate prospective diagnostic is justified by this identified scheduling
defect, not an automatic retry or an attempt to hide 22378. It must use fresh
output, current source pins and a new resource-protocol source binding with
unchanged endpoints. Passing it would still not prove causal observer overhead,
native environmental handoff, continuous containment or host isolation. No
production scaling identity is launched by this fix.
