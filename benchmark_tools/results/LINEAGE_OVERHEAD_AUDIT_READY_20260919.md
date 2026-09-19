# Lineage Overhead Audit Ready

The existing complete-panel auditor now has a separate pinned `lineage_21999`
registry entry. It binds array 21999, the exact plan/recipe/authorization,
actual local submission-script path, native worker/collector, original inputs,
runtime records and receipts. Historical panel identities and the default
audit selection remain unchanged. This is local postprocessing only: no
running DGX recipe was modified or redeployed.

## Evidence Handling

Periodic tasks read `lineage_report.json` and invoke the lineage raw replay.
Boundary tasks read `lineage_boundary_report.json` and invoke its separate
raw replay. Both retain complete lineage screening. Periodic tasks retain
original and narrow interval flags; boundary flags remain null with explicit
unavailable interval coverage. Periodic pressure summaries are labeled as
the observation window, not exact native-command boundaries. Boot identity
uses the lineage evidence, not the absent historical frontier inventory.

The audit requires all 18 accounting rows and detailed terminal controller
records, with matching states/exit codes, before native/context inspection.
It independently invokes native-output validation and canonical fingerprints,
checks archived recipe identity, preserves failed/missing/invalid tasks and
partial lineage reports, and applies unchanged signed nine-pair arithmetic
and numerical budgets. Missing pairs yield an incomplete overall budget
conclusion. Numerical budget results never establish environmental validity,
and boundary-only evidence cannot exclude transient competing activity.

## Post-Run Invocation

Do not collect or inspect DGX native outputs until all 18 tasks are terminal.
After recorder 22000 completes, independently replay the full controller
capture. Collect the complete archive preserving paths relative to the DGX
project root, including deployed recipe, native runs, original scaling inputs
and detailed records copied as `scheduler_0.txt` through `scheduler_17.txt`.
Keep the original controller capture and raw polls separately.

Accounting must use the existing seven-column format
`JobID|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList`, with array task IDs
such as `21999_0`, not numeric `JobIDRaw` alone. Then use a fresh output path:

```sh
python -B -m benchmark_tools.audit_frontier_overhead \
  --panel lineage_21999 \
  --archive benchmarks/work/lineage_overhead_archive_21999 \
  --results benchmark_tools/results \
  --accounting benchmarks/work/lineage_overhead_array_accounting_21999.txt \
  --output benchmarks/work/lineage_overhead_audit_21999.json
```

## Verification And Limits

All 355 focused tests pass across four-panel provenance, archive orchestration,
paired arithmetic, lineage periodic/boundary replay and historical replay.
Coverage includes every assigned task, changed source/authorization/input/
runtime/worker/scheduler evidence, cross-panel substitution, real context
checksums, missing terminal records, contradictory accounting, correct lineage
report/replay dispatch, retained flags, unavailable boundary coverage and
failed-pair handling. Orchestration fixtures mock native products and raw
replay; these tests do not prove that the live outputs pass the audit.

Latest controller check: task 21999_0 RUNNING at 6:00, tasks 1-17 pending,
recorder 22000 RUNNING at 6:41; BLAST 21713 RUNNING at 6:12:34. No native
output was inspected. Runtime validity, actual overhead, remaining CPU
discrepancies, during-read lifecycle behavior and scientific inclusion remain
unresolved. The complete publication goal is not achieved.
