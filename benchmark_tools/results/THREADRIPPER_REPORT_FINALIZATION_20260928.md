# Reporting-Stage Accounting

The v4 Threadripper collector retains an additional
`report_finalization.json` after lineage/context report construction and
serialization. It records monotonic start/end times, elapsed reporting time
and a subsequent job-cgroup memory observation. Native timer boundaries,
commands, resources, release guard, samples and scientific settings are
unchanged. Historical v1/v2/v3 measurements are not rewritten.

This closes the previously missing reporting-stage observation, not all
end-of-job accounting: the new read excludes later runtime/output validation,
the receipt's own serialization and job teardown. It records cumulative job
peak, not reporting-only peak. Native-step and job peaks overlap; do not add
them or subtract a baseline. Elapsed reporting time excludes the final memory
read and receipt serialization. A context-manager failure still attempts the
memory read and preserves the original exception; missing evidence after an
abrupt kill remains incomplete, never successful reporting.

Independent replay requires the new receipt for v4 and validates job/scope,
native-to-report phase ordering, elapsed arithmetic, raw memory gauges/events,
the 128-GiB limit and nondecreasing peaks/events. Historical schemas retain
their existing requirements and return no finalization record. Filename
presence alone is insufficient. This adds no timing admission authority.

## Live Check

[Diagnostic job 22366](threadripper_reporting_probe_22366.json) completed
0:0 in five scheduler seconds. It ran `/usr/bin/sleep 1`, not OrthoHMM or a
production timing identity. The fixed allocation and release-budget guard
were used. Independent replay and terminal controller validation passed.

- Native wall: 1.018798277 seconds.
- Report construction/serialization: 0.077414250 seconds.
- Job peak at both post-native and post-report reads: 38,035,456 bytes.
- Post-report memory event counters report no OOM or OOM-kill events.
- Maximum observed foreign interval-average CPU use: 72.5394 cores.

The unchanged peak is a valid observation for this tiny workload, not evidence
that large reports cost no memory. The earlier full-day report-size projection
remains relevant. The host was not quiet; no measured time is admitted to the
controlled comparison. Raw observations, reports and submission provenance
remain under `benchmarks/work/threadripper_reporting_probe_20260928` and the
adjacent submission receipt. No unrelated process was signalled.

## Tests And Remaining Integration

73 focused tests passed after correcting an order-dependent memory-test
measurement. Two combined runs initially failed the existing disk-observation
test's retained-allocation bound (~2.65 MB); its isolated run passed. A third
diagnostic reproduced the failure and attributed a single 2,560-KiB allocation
to `pathlib.py:74`, `sys.intern`, rather than retained observation payloads.
That allocation depends on the shared interpreter's earlier intern-table
population. The memory test now uses a fresh isolated Python process with
the original 200 x 100,000-byte payload workload and unchanged 1/2-MB bounds.
The new reporting tests cover success, reporting failure, memory-read failure,
simultaneous failures, timing/identity corruption and v4 disk-backed replay.

The collector and replay source changes invalidate their old runtime hashes.
The v3 runtime binding and v2 lookup receipt remain historical artifacts;
do not bypass their expected drift rejection. Regenerate and review the
binding, include the new reporting helper, and rerun checked native fixtures
before any production execution. Production orchestration, quiet-window
eligibility, large-workload observer validation and complete-job resource
accounting remain open. No production identity was launched.
