# Full Integrated OrthoBench Reproduction Protocol

Frozen before submission. The [machine-readable protocol](integrated_full_ob_protocol_20260927.json)
pins one integrated workflow execution on all 251,378 genes, 12 input species,
70 reference families and 11 low-certainty files. This validates the combined
controller at dataset scale after its successful installation-fixture run.
It is a selective reproducibility rerun, not tuning or a replacement for the
earlier admitted full run 22326.

Submission update: protocol commit `8ae66d1` preceded
[job 22337](integrated_full_ob_submission_20260927.json). Its latest observed
state is RUNNING; installations completed and built-in HMM search started.
This update does not change the frozen endpoints or imply completion/admission.

## Fixed Execution

- Local `bizon` node, 32 CPUs, 128 GiB requested RAM and 24-hour job limit.
- Shared-host reproduction only; no controlled timing claim and no DGX access.
- Exactly one attempt, fresh inference/reader installations, no native checkpoint reuse.
- Controller bytes from `4d037b1`, frozen scientific package and canonical
  candidate ordering as in the validated recovery workflow.
- Separate inference and reader hash locks. No reference files passed to inference.
- Persistent private copies of wheels, entrypoint, readers and native tools.
  Two MAFFT convenience links are relative in the new copy; original assets remain intact.
- Full native inference, four independent scientific readers and frozen scoring.
- Stop on failure, preserve logs and outputs, and do not automatically retry.

The local run root is
`benchmarks/work/publication_integrated_full_ob_20260927`.
Its 64,737-byte plan pins 215 file records and has SHA256
`5cab1cd592e6f6d69276c84d0ab5e1b9dc91e1e12a44ace757db87d08df2d899`.
The scheduler launcher checks the plan, its own source and every recorded file
before execution and again after successful completion. Native outputs are
written under a fresh `run/` directory; existing attempts are refused.
The trusted installer is the separately validated patched reader environment.

## Prespecified Admission And Comparison

Admission requires successful terminal scheduler accounting; exact command,
plan, input, source, tool and lock identities; all eight completed stages;
successful structure, sequence, event and hierarchy readback; and a full
partition covering each input gene exactly once. Independently check installed
inventories and wheel payloads after completion.

Compare the complete root partition without group-label dependence and the
entire 70-family score object against admitted run 22326. The expected retained
F1 is 74.10607351873405%, with precision 81.77045380181866% and recall
67.75533630827974%. Exact equality is a reproduction endpoint, not an excuse
to discard changed results. Report any changed groups, genes, family records
or metrics; preserve all prior outputs and scores. No result-driven ordering,
threshold, dependency or default selection is allowed.

GNU time covers the entire combined workflow, including installations and
readback/scoring. Its maximum process RSS is not simultaneous process-tree
memory. These descriptive measurements must not be mixed with the pending
dedicated-host scaling evidence or labeled inference-only timings.

## Remaining Requirements

Completion of this run alone does not establish publication readiness,
cross-host portability, new biological validation, QfO uncertainty estimates,
controlled timing, complete redistribution clearance or public deposition.
The independent post-completion admission must be recorded separately from
the controller's success marker. Nineteen focused launcher/controller tests
pass before submission.
