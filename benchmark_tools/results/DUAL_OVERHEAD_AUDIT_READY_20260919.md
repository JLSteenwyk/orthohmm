# Dual-Collector Overhead Audit Ready

Extended the existing overhead audit with a separate checksum-pinned
`dual_21920` registry entry. The historical frontier/pressure panel identities
and default audit choice remain unchanged. This code is local postprocessing;
the live recipe exported from `6e222a4` is not edited or redeployed.

The new entry binds the exact18-task plan, recipe, authorization, submitted
batch path, native worker/collector launch, inputs/order, runtime evidence
and the additional recipe hash in task receipts. Periodic tasks use the
dual raw replay; boundary tasks retain the pressure-enabled boundary replay.
Native output validation and canonical within-pair fingerprints remain the
existing independent checks.

Periodic results retain original interval flags and complete dual screening,
including narrow flags. The arithmetic summary continues to use signed
overhead and frozen5% median/10% pair budgets; it does not turn narrow-screen
passage into environmental or scientific admission. Failed-task evidence now
also retains any partial dual report.

Before reading native files for this panel, require all18 accounting rows
terminal and all18 detailed controller records. Their terminal states and
exit codes must agree. Missing/contradictory records block the audit rather
than silently treating unobserved tasks as successful or failed. Boundary
arms still have no interval screening.

## Invocation After Complete Collection

Archive layout preserves remote paths relative to
`/home/jlsteenwyk/projects/orthohmm-publication`, with
`scheduler_0.txt` through `scheduler_17.txt` copied from the completed
controller capture. Include the recipe tree, native runs, original scaling
inputs and every result/log needed by output validation. Do not inspect the
live DGX archive during the quiet panel.

```bash
python -m benchmark_tools.audit_frontier_overhead \
  --panel dual_21920 \
  --archive benchmarks/work/dual_overhead_archive_21920 \
  --results benchmark_tools/results \
  --accounting benchmarks/work/dual_overhead_array_accounting_21920.txt \
  --output benchmarks/work/dual_overhead_audit_21920.json
```

Accounting rows use the existing seven-column format:
`JobID|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList`.
Use `JobID` so rows carry array task identifiers such as `21920_3`;
`JobIDRaw` alone provides numeric allocation identifiers and cannot match
the task-index gate. This corrects the earlier preparation-note typo;
the auditor and original accounting evidence are unchanged.
Retain the controller capture and raw polls separately and independently
replay the first-terminal observations before interpreting the final audit.

## Verification

247 focused tests passed across provenance, audit orchestration, paired
arithmetic and both raw replay implementations. Tests cover all18 new task
identities, legacy panels, dual versus boundary dispatch, separate original
and narrow flags, missing terminal records, accounting disagreement and
retention of failed assigned tasks. Synthetic orchestration tests do not
establish production validity; real audit must follow complete execution.

At preparation, task21920_0 remained RUNNING at5:11 and recorder21922 at4:41.
Corrected full QfO OrthoFinder21706_1 remained RUNNING at10:53:28. No result
or environmental-validity claim is made by this prepared audit.
