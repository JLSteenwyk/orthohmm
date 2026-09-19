# Prospective Native Root-Context Integration

## Scope

Run the three existing frozen four-proteome native commands, in order:
OrthoHMM high sensitivity, OrthoHMM satellite_v2, and full OrthoFinder 3.1.5.
Derive them from the same pressure-panel v2 tasks 1, 3 and 8 used for lineage
diagnostics. Preserve the 73,266 proteins, original enumeration, exact method
settings, inference source and runtime manifests. Relocate only output/cache
paths and diagnostic metadata into a fresh `root_context_native_v1` tree.

The question is whether separately timestamped root/system/user/init context
can be collected and replayed during actual native tools, while retaining
unchanged output semantics, original/narrow CPU screens, native timer,
pressure and native-step cgroup memory accounting. The completed synthetic
controls 22020 do not establish this native integration or its overhead.

## Execution And Evidence

Use a single exclusive DGX20CPU/96GiB allocation with a one-hour limit and
no requeue, containing the three serial native steps. Keep each native
timeout at 900 seconds and observation cadence at one second. Retain all
outcomes and stop subsequent execution on a failed preparation, command,
measurement or cleanup; unrun entries must remain explicit. No selective
retry or substitution. A failure in the scientific output audit must not be
reinterpreted as an executable success for comparative timing.

Freeze the derived plan and complete deployment recipe before submission.
Verify an empty DGX scheduler queue, and use a bounded submitting login
session lasting until the whole job terminates, following the demonstrated
held-session mechanism. The submission/client bound must exceed the one-hour
job limit and be recorded before launch. No extra SSH inspection or transfers
during the panel. Retain the session/client and unrelated user-manager
services as parts of the environment; do not change persistent login settings
or stop unrelated services. Retrieve the bounded manager journal afterward.

Use the existing `measure_run` preparation and before/after input/runtime
checks, with `measure_native_root_context.measure_native_run` as adapter.
Keep `lineage_report.json` as the workflow measurement and separately retain
`root_context_report.json`; never collapse the supplementary windows into
the earlier lineage read windows. Reuse independent replay for both reports.

After termination, audit every scheduler, source, input, runtime, native
command and output witness. Validate native output semantics and compare
canonical work fingerprints to the existing lineage diagnostics, retaining
differences rather than assuming identity. Check within- and between-run
boot/scope identity and observation order; record missing coverage explicitly.

## Reporting And Limits

Report all three statuses, native wall/CPU/memory/pressure, canonical output
comparison, and complete original/narrow flag inventories. Describe named
scope CPU, signed residuals, host categories, observed membership changes
and read durations for all intervals and flagged/unflagged subsets. Do not
sum overlapping guest/user host categories, treat dependent intervals as
independent replicates, infer causality from co-occurrence, or remove flags.

These are one-off native integration diagnostics, not matched-resource
scaling estimates or a method-speed ranking. Native-tool context-collector
overhead still needs a separate prespecified comparison. Non-CPU isolation,
background-service limitations and prospective scientific timing inclusion
remain unresolved. No method default, threshold or historical result changes
under this protocol, and success cannot establish publication readiness.
